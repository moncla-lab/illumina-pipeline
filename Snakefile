from mlip import *

import os
import sys
import csv
import json
from itertools import product

import pandas as pd


configfile: "config.yml"

# Validate that analysis directory is specified
if 'analysis' not in config or not config['analysis']:
    print("ERROR: The 'analysis' key in 'config.yml' cannot be empty.")
    print("Please specify a name for your analysis (e.g., 'h5n1_cattle').")
    sys.exit(1)

ANALYSIS_DIR = config['analysis']

def data(path):
    """Construct paths within the analysis directory"""
    return f"{ANALYSIS_DIR}/{path}"

# Check for file manifest in the analysis directory
if not os.path.exists(data("file_manifest.json")):
    print(f"ERROR: '{data('file_manifest.json')}' not found.")
    print("Please configure the pipeline. You can see your status by running:")
    print("  python mlip/dataflow.py check")
    sys.exit(1)


wildcard_constraints:
  segment="[^/]+",
  sample="[^/]+",
  replicate="[^/]+",
  mapping_stage="[^/]+"

REFERENCE = config['reference']
USING_ZIP = REFERENCE.endswith(".zip")
reference_dictionary = load_reference_dictionary(REFERENCE, USING_ZIP)
metadata_dictionary = load_metadata_dictionary(ANALYSIS_DIR)
SAMPLES = samples_to_analyze(ANALYSIS_DIR)
DUPLICATE_SAMPLES = get_duplicate_samples(metadata_dictionary)
SEGMENTS = load_segments(config, ANALYSIS_DIR)
NUMBER_OF_REMAPPINGS = config['number_of_remappings']

if not USING_ZIP:
    rule fetch_reference_data:
        message:
            'Fetching reference data for segment {wildcards.segment}...'
        output:
            fasta=data('reference/{segment}/sequence.fasta'),
            genbank=data('reference/{segment}/metadata.gb')
        resources:
            ncbi_fetches=1
        params:
            genbank_accession=(
                lambda wildcards:
                reference_dictionary[wildcards.segment]['genbank_accession']
            ),
        shell:
            '''
                efetch -db nuccore \
                    -id {params.genbank_accession} \
                    -format genbank \
                    > {output.genbank}

                efetch -db nuccore \
                    -id {params.genbank_accession} \
                    -format fasta \
                | seqkit replace -p "^(.+)" -r "{wildcards.segment} genbank"\
                    > {output.fasta}
            '''

rule build_full_reference:
    message:
        'Concatenating reference data into single FASTA...'
    input:
        expand(data("reference/{segment}/sequence.fasta"), segment=SEGMENTS)
    output:
        data('reference/sequences.fasta'),
    shell:
        'cat {input} > {output}'


def get_genbank_input(wildcards):
    if USING_ZIP:
        return data(f"reference/{wildcards.segment}/metadata.gb")
    else:
        return rules.fetch_reference_data.output.genbank


rule genbank_to_gtf:
    message:
        'Converting Genbank data to GTF...'
    input: get_genbank_input
    output:
        data('reference/{segment}/metadata.gtf')
    run:
        genbank_to_gtf(input[0], output[0], wildcards.segment)

rule full_gtf:
    input:
        expand(data('reference/{segment}/metadata.gtf'), segment=SEGMENTS)
    output:
        data('reference/metadata.gtf')
    shell:
        'cat {input} > {output}'

rule gene_list:
    input:
        rules.full_gtf.output[0]
    output:
        data('reference/gene_list.txt')
    run:
        extract_genes(input[0], output[0])

rule sample_list:
    input:
        data('metadata.tsv')
    output:
        data('sample_list.txt')
    shell:
        'csvcut -t -c SampleId {input} | sort | uniq | grep -v SampleId > {output}'

def forward_fastq_merge_inputs(wildcards):
    experiments = metadata_dictionary[wildcards.sample][wildcards.replicate]
    forward_path = data('%s/sequencing-{sequencing}/forward.fastq.gz') % wildcards.sample
    result = expand(
        forward_path,
        sequencing=metadata_dictionary[wildcards.sample][wildcards.replicate]
    )
    return result


def reverse_fastq_merge_inputs(wildcards):
    experiments = metadata_dictionary[wildcards.sample][wildcards.replicate]
    reverse_path = data('%s/sequencing-{sequencing}/reverse.fastq.gz') % wildcards.sample
    return expand(
        reverse_path,
        sequencing=metadata_dictionary[wildcards.sample][wildcards.replicate]
    )


# Helper function to determine input files for the consensus summary comparison
def get_consensus_summary_inputs(wildcards):
    final_remapping_num = NUMBER_OF_REMAPPINGS
    
    if final_remapping_num < 1:
        return []

    final_path = data(f"{wildcards.sample}/replicate-{wildcards.replicate}/remapping-{final_remapping_num}/consensus.fasta")

    if final_remapping_num == 1:
        penultimate_path = data(f"{wildcards.sample}/replicate-{wildcards.replicate}/initial/consensus.fasta")
    else:
        penultimate_num = final_remapping_num - 1
        penultimate_path = data(f"{wildcards.sample}/replicate-{wildcards.replicate}/remapping-{penultimate_num}/consensus.fasta")
        
    return [penultimate_path, final_path]


rule concatenate_replicates_from_manifest:
    input:
        manifest=data("file_manifest.json")
    output:
        forward=temp(data("{sample}/replicate-{replicate}/forward.fastq")),
        reverse_=temp(data("{sample}/replicate-{replicate}/reverse.fastq"))
    run:
        concatenate_replicates_from_manifest_py(
            manifest_filepath=input.manifest,
            sample_id=wildcards.sample,
            replicate_num_str=wildcards.replicate,
            output_forward_fastq_path=output.forward,
            output_reverse_fastq_path=output.reverse_
        )

rule trimmomatic:
    message:
        '''
            Trimming replicate {wildcards.replicate} of sample {wildcards.sample}...
            Parameters:
                Window size: {params.trimming_window_size}
                Q-score: {params.minimum_quality_score}
                Minimum length: {params.trimming_minimum_length}
        '''
    input:
        forward=rules.concatenate_replicates_from_manifest.output.forward,
        reverse_=rules.concatenate_replicates_from_manifest.output.reverse_
    output:
        forward_paired=temp(data('{sample}/replicate-{replicate}/forward_paired.fastq')),
        reverse_paired=temp(data('{sample}/replicate-{replicate}/reverse_paired.fastq')),
        forward_unpaired=temp(data('{sample}/replicate-{replicate}/forward_unpaired.fastq')),
        reverse_unpaired=temp(data('{sample}/replicate-{replicate}/reverse_unpaired.fastq')),
        stdout=data('{sample}/replicate-{replicate}/trimmomatic-stdout.txt'),
        log=data('{sample}/replicate-{replicate}/trimmomatic.log'),
    params: **config
    priority: 1
    shell:
        '''
            trimmomatic PE \
                {input.forward} {input.reverse_} \
                {output.forward_paired} {output.forward_unpaired} \
                {output.reverse_paired} {output.reverse_unpaired} \
                SLIDINGWINDOW:{params.trimming_window_size}:{params.minimum_quality_score} \
                MINLEN:{params.trimming_minimum_length} \
                > {output.stdout} 2> {output.log}
        '''

#rule vapor_segment:
#    input:
#        fastq=rules.trimmomatic.output.concat,
#        reference_db=data('reference/{segment}/all.fasta'),
#        mlip_reference=data('reference/{segment}/sequence.fasta')
#    output:
#        vapor_reference=data('{sample}/replicate-{replicate}/initial_reference_{segment}.fa'),
#        vapor_id=data('{sample}/replicate-{replicate}/initial_id_{segment}.txt'),
#        unaligned=data('{sample}/replicate-{replicate}/initial_unaligned_{segment}.fasta'),
#        aligned=data('{sample}/replicate-{replicate}/initial_aligned_{segment}.fasta')
#    params:
#        data('{sample}/replicate-{replicate}/initial_reference_{segment}')
#    shell:
#        '''
#            vapor.py -m .01 -fq {input.fastq} -fa {input.reference_db} -o {params}
#            head -n 1 {output.vapor_reference} | cut -c 2- > {output.vapor_id}
#            sed -i '1s/.*/>{wildcards.segment} vapor/' {output.vapor_reference}
#            cat {input.mlip_reference} {output.vapor_reference} > {output.unaligned}
#            mafft --preservecase {output.unaligned} > {output.aligned}
#        '''
#
#rule hybrid_segment:
#    input:
#        rules.vapor_segment.output.aligned
#    output:
#        data('{sample}/replicate-{replicate}/hybrid_{segment}.fasta'),
#    run:
#        fill(input[0], output[0])
#
#rule hybrid_reference:
#    input:
#        expand(
#            data('{{sample}}/replicate-{{replicate}}/hybrid_{segment}.fasta'),
#            segment=SEGMENTS
#        )
#    output:
#        data('{sample}/replicate-{replicate}/reference.fasta')
#    shell:
#        'cat {input} > {output}'

def situate_reference_input(wildcards):
    if wildcards.mapping_stage == 'initial':
        if config['use_vapor']:
            return data(f'{wildcards.sample}/replicate-{wildcards.replicate}/reference.fasta')
        else:
            return data('reference/sequences.fasta'),
    elif wildcards.mapping_stage == 'remapping-1':
        return data(f'{wildcards.sample}/replicate-{wildcards.replicate}/initial/filler.fasta')
    mapping_stage_int = int(wildcards.mapping_stage.split('-')[1]) - 1
    return data(f'{wildcards.sample}/replicate-{wildcards.replicate}/remapping-{mapping_stage_int}/filler.fasta')


rule situate_reference:
    input:
        situate_reference_input
    output:
        data('{sample}/replicate-{replicate}/{mapping_stage}/reference/sequences.fasta')
    shell:
        'cp {input} {output}'

rule index:
    message:
        'Indexing reference sequence...'
    input:
        rules.situate_reference.output[0]
    params:
        data('{sample}/replicate-{replicate}/{mapping_stage}/reference/index')
    output:
        index1=data('{sample}/replicate-{replicate}/{mapping_stage}/reference/index.1.bt2'),
        index2=data('{sample}/replicate-{replicate}/{mapping_stage}/reference/index.2.bt2'),
        index3=data('{sample}/replicate-{replicate}/{mapping_stage}/reference/index.3.bt2'),
        index4=data('{sample}/replicate-{replicate}/{mapping_stage}/reference/index.4.bt2'),
        indexrev1=data('{sample}/replicate-{replicate}/{mapping_stage}/reference/index.rev.1.bt2'),
        indexrev2=data('{sample}/replicate-{replicate}/{mapping_stage}/reference/index.rev.2.bt2'),
        stdout=data('{sample}/replicate-{replicate}/{mapping_stage}/reference/bowtie2-stdout.txt'),
        stderr=data('{sample}/replicate-{replicate}/{mapping_stage}/reference/bowtie2-stderr.txt')
    shell:
        'bowtie2-build {input} {params} > {output.stdout} 2> {output.stderr}'


rule mapping:
    message:
        '''
            Mapping replicate {wildcards.replicate} of sample {wildcards.sample} 
            at stage {wildcards.mapping_stage} to reference...
        '''
    input:
        forward_paired=temp(rules.trimmomatic.output.forward_paired),
        reverse_paired=temp(rules.trimmomatic.output.reverse_paired),
        forward_unpaired=temp(rules.trimmomatic.output.forward_unpaired),
        reverse_unpaired=temp(rules.trimmomatic.output.reverse_unpaired),
        index=rules.index.output.index1
    params:
        data('{sample}/replicate-{replicate}/{mapping_stage}/reference/index')
    output:
        sam=temp(data('{sample}/replicate-{replicate}/{mapping_stage}/mapped.sam')),
        stdout=data('{sample}/replicate-{replicate}/{mapping_stage}/bowtie2-stdout.txt'),
        stderr=data('{sample}/replicate-{replicate}/{mapping_stage}/bowtie2-stderr.txt')
    priority: 2
    shell:
        '''
            bowtie2 --local --very-sensitive-local -x {params} \
                -1 {input.forward_paired} -2 {input.reverse_paired} \
                -U {input.forward_unpaired},{input.reverse_unpaired} \
                -S {output.sam} \
                > {output.stdout} 2> {output.stderr}
        '''

rule samtools:
    message:
        'Running various samtools modules on {wildcards.replicate} of sample {wildcards.sample}...'
    input:
        sam=rules.mapping.output.sam,
        reference=rules.situate_reference.output[0]
    output:
        mapped=temp(data('{sample}/replicate-{replicate}/{mapping_stage}/mapped.bam')),
        sorted_=temp(data('{sample}/replicate-{replicate}/{mapping_stage}/sorted.bam')),
        index=data('{sample}/replicate-{replicate}/{mapping_stage}/sorted.bam.bai'),
        stdout=data('{sample}/replicate-{replicate}/{mapping_stage}/samtools-stdout.txt'),
        pileup=temp(data('{sample}/replicate-{replicate}/{mapping_stage}/samtools.pileup')),
        stderr=data('{sample}/replicate-{replicate}/{mapping_stage}/samtools-stderr.txt')
    priority: 3
    shell:
        '''
            samtools view -S -b {input.sam} > {output.mapped} 2> {output.stderr}
            samtools sort {output.mapped} -o {output.sorted_} > {output.stdout} 2>> {output.stderr}
            samtools index {output.sorted_} >> {output.stdout} 2>> {output.stderr}
            samtools mpileup -a -A -d 0 -B -Q 0 \
                -f {input.reference} {output.sorted_} > {output.pileup} 2>> {output.stderr}
        '''

rule call_variants:
    message:
        '''
            Calling variants on replicate {wildcards.replicate} of sample {wildcards.sample}...
            Parameters:
                Mapping stage: {wildcards.mapping_stage}
                SNP frequency: {params.variants_minimum_frequency}
                Minimum coverage for variant calling: {params.variants_minimum_coverage}
                Strand filter: {params.strand_filter}
                SNP quality threshold: {params.minimum_quality_score}
        '''
    input:
        pileup=rules.samtools.output.pileup,
        stderr=rules.samtools.output.stderr,
        reference=situate_reference_input
    output:
        vcf=    data('{sample}/replicate-{replicate}/{mapping_stage}/varscan.vcf'),
        tsv=    data('{sample}/replicate-{replicate}/{mapping_stage}/varscan.tsv'),
        vcf_zip=data('{sample}/replicate-{replicate}/{mapping_stage}/varscan.vcf.gz'),
        index=  data('{sample}/replicate-{replicate}/{mapping_stage}/varscan.vcf.gz.tbi'),
        stderr= data('{sample}/replicate-{replicate}/{mapping_stage}/varscan-stderr.txt')
    params:
        **config
    priority: 4
    shell:
        '''
        (
            varscan mpileup2snp {input.pileup} \
                --min-coverage {params.variants_minimum_coverage} \
                --min-avg-qual {params.minimum_quality_score} \
                --min-var-freq {params.variants_minimum_frequency} \
                --strand-filter {params.strand_filter} \
                --output-vcf 1 > {output.vcf} 2> {output.stderr}
            grep -v '^##' {output.vcf} > {output.tsv}
            bgzip -c {output.vcf} > {output.vcf_zip}
            tabix -p vcf {output.vcf_zip}
        ) || true
        '''

rule coverage:
    message:
        'Computing coverage of replicate {wildcards.replicate} for sample {wildcards.sample}...'
    input:
        rules.samtools.output.sorted_
    output:
        bg= data('{sample}/replicate-{replicate}/{mapping_stage}/coverage.bedGraph'),
        tsv=data('{sample}/replicate-{replicate}/{mapping_stage}/coverage.tsv')
    priority: 4
    shell:
        '''
            echo "segment\tstart\tend\tcoverage" > {output.tsv}
            bedtools genomecov -ibam {input} -bga > {output.bg}
            cat {output.bg} >> {output.tsv}
        '''

rule coverage_summary:
    message:
        'Computing coverage summary of replicate {wildcards.replicate} for sample {wildcards.sample}...'
    input:
        rules.coverage.output.tsv
    output:
        data('{sample}/replicate-{replicate}/{mapping_stage}/coverage-report.tsv')
    run:    
        compute_coverage_categories_io(input[0], output[0])

rule call_segment_consensus:
    input:
        bam=rules.samtools.output.sorted_,
        pileup=rules.samtools.output.pileup,
        reference=situate_reference_input,
        original_reference=data('{sample}/replicate-{replicate}/initial/reference/sequences.fasta')
    output:
        ivar_fasta=data('{sample}/replicate-{replicate}/{mapping_stage}/segments/{segment}/ivar.fa'),
        fasta=data('{sample}/replicate-{replicate}/{mapping_stage}/segments/{segment}/consensus.fasta'),
        reference=data('{sample}/replicate-{replicate}/{mapping_stage}/segments/{segment}/reference.fasta'),
        bam=temp(data('{sample}/replicate-{replicate}/{mapping_stage}/segments/{segment}/segment.bam')),
        bai=data('{sample}/replicate-{replicate}/{mapping_stage}/segments/{segment}/segment.bam.bai'),
        samtools=data('{sample}/replicate-{replicate}/{mapping_stage}/segments/{segment}/samtools.fasta'),
        unaligned=data('{sample}/replicate-{replicate}/{mapping_stage}/segments/{segment}/unaligned.fasta'),
        aligned=data('{sample}/replicate-{replicate}/{mapping_stage}/segments/{segment}/aligned.fasta')
    params: ** { \
        **config, \
        'ivar': data('{sample}/replicate-{replicate}/{mapping_stage}/segments/{segment}/ivar') \
    }
    priority: 5
    shell:
        '''
        (
            seqkit grep -p {wildcards.segment} {input.reference} > {output.reference}
            samtools faidx {output.reference}
            if [ ! -s {output.reference} ]; then
                echo "WARNING: empty reference, this is just so the pipeline runs batches to completion"
                cp {input.original_reference} {output.reference}
            fi
            samtools view -b -h {input.bam} {wildcards.segment} > {output.bam}
            samtools index {output.bam}
            grep {wildcards.segment} {input.pileup} | ivar consensus -p {params.ivar} \
                -m {params.consensus_minimum_coverage} \
                -q 0 \
                -t {params.consensus_minimum_frequency} \
                -c {params.consensus_minimum_frequency}
            echo ">{wildcards.segment} {wildcards.mapping_stage}" > {output.fasta}
            tail -n +2 {output.ivar_fasta} >> {output.fasta}
            seqkit grep -p {wildcards.segment} {input.original_reference} > {output.unaligned}
            cat {output.fasta} >> {output.unaligned}
            samtools consensus --mode simple -d {params.consensus_minimum_coverage} --call-fract 0 {output.bam} > {output.samtools}
            mafft --preservecase {output.unaligned} > {output.aligned}
        ) || true
        '''

rule full_consensus:
    input:
        expand(
            data('{{sample}}/replicate-{{replicate}}/{{mapping_stage}}/segments/{segment}/consensus.fasta'),
            segment=SEGMENTS
        )
    output:
        data('{sample}/replicate-{replicate}/{mapping_stage}/consensus.fasta')
    priority: 6
    shell:
        'cat {input} > {output}'

rule fill_consensus:
    input:
        rules.call_segment_consensus.output.aligned
    output:
        data('{sample}/replicate-{replicate}/{mapping_stage}/segments/{segment}/filler.fasta'),
    run:
        fill(input[0], output[0])


rule full_filler:
    input:
        expand(
            data('{{sample}}/replicate-{{replicate}}/{{mapping_stage}}/segments/{segment}/filler.fasta'),
            segment=SEGMENTS
        )
    output:
        data('{sample}/replicate-{replicate}/{mapping_stage}/filler.fasta')
    shell:
        'cat {input} > {output}'


def call_sample_consensus_input(wildcards):
    remapping_string = f'remapping-{NUMBER_OF_REMAPPINGS}'
    return expand(
        data('{{sample}}/replicate-{replicate}/%s/consensus.fasta') % remapping_string,
        replicate=metadata_dictionary[wildcards.sample].keys()
    )


rule call_sample_consensus:
    input: call_sample_consensus_input
    output:
        data('{sample}/consensus.fasta')
    run:
        call_sample_consensus(input, output[0])

rule call_sample_proteins:
    input:
        rules.call_sample_consensus.output[0],
        expand(data('reference/{segment}/metadata.gb'), segment=SEGMENTS)
    output:
        directory(data('{sample}/protein'))
    run:
        translate_consensus_genes(input[0], output[0], wildcards.sample, ANALYSIS_DIR)

#rule multiqc:
#    message:
#        'Running Multi QC on {wildcards.replicate} of sample {wildcards.sample}...'
#    input:
#        rules.trimmomatic.output.log,
#        rules.samtools.output.stats,
#        rules.samtools.output.flagstat,
#        rules.samtools.output.depth
#    output:
#        data('{sample}/replicate-{replicate}/{mapping_stage}/multiqc_report.html')
#    params:
#        data('{sample}/replicate-{replicate}/{mapping_stage}')
#    shell:
#        'multiqc -f {params} --outdir {params}'

rule coding_regions:
    input:
        annotated_references=expand(data('reference/{segment}/metadata.gb'), segment=SEGMENTS),
        replicate_consensus=rules.full_consensus.output[0]
    output:
        data('{sample}/replicate-{replicate}/{mapping_stage}/coding_regions.json')
    run:
        extract_coding_regions_io(SEGMENTS, input.replicate_consensus, output[0], ANALYSIS_DIR)

rule annotate_varscan:
    input:
        coding_regions=rules.coding_regions.output[0],
        varscan=rules.call_variants.output.vcf
    output:
        data('{sample}/replicate-{replicate}/{mapping_stage}/varscan-annotated.tsv')
    run:
        with open(input.coding_regions) as json_file:
            coding_regions = json.load(json_file)
        annotate_amino_acid_changes(
            coding_regions, input.varscan, output[0]
        )

rule clean_varscan:
    message:
        'Cleaning varscan VCF from replicate {wildcards.replicate} of sample {wildcards.sample}...'
    input:
        rules.call_variants.output.tsv
    output:
        data('{sample}/replicate-{replicate}/{mapping_stage}/ml.tsv')
    run:
        df = pd.read_csv(input[0], sep='\t')
        clean_varscan(df).to_csv(output[0], sep='\t', index=False)

def merge_varscan_inputs(wildcards):
    return expand(
        data('{{sample}}/replicate-{replicate}/remapping-%s/varscan-annotated.tsv') % NUMBER_OF_REMAPPINGS,
        replicate=range(1, len(metadata_dictionary[wildcards.sample])+1)
    )

rule merge_varscan_across_replicates:
    message:
        'Merging variant calls of sample {wildcards.sample}...'
    input: merge_varscan_inputs
    output:
        data('{sample}/ml.tsv')
    run:
        merge_varscan_io(input, output[0])

rule visualize_replicate_calls:
    message:
        'Visualizing replicate variant calls of sample {wildcards.sample}...'
    input:
        rules.merge_varscan_across_replicates.output[0]
    output:
        data('{sample}/ml.html')
    run:
        replicate_variant_plot(input[0], output[0])


def full_coverage_summary_input(wildcards):
    coverage_filepaths = []
    for sample, replicates in metadata_dictionary.items():
        for replicate in replicates.keys():
            coverage_filepaths.append(
                data(f'{sample}/replicate-{replicate}/remapping-{NUMBER_OF_REMAPPINGS}/coverage-report.tsv')
            )
    return coverage_filepaths


rule full_coverage_summary:
    input: full_coverage_summary_input
    output:
        data('coverage-report.tsv'),
    run:
        coverage_summary(input, output[0])

rule full_genome:
    input:
        expand(data('{sample}/consensus.fasta'), sample=SAMPLES)
    output:
        data('{segment}.fasta')
    params:
        samples=' '.join(SAMPLES)
    shell:
        '''
          for sample in {params.samples}; do
            consensus_file={ANALYSIS_DIR}/$sample/consensus.fasta
            seqkit grep -p {wildcards.segment} $consensus_file | \
              seqkit replace -p {wildcards.segment} -r "$sample {wildcards.segment}" >> \
              {output}
          done
        '''

#rule check_replicate_consensus:
#    input:
#        fasta=rules.call_segment_consensus.output.fasta,
#        pileup=rules.call_segment_consensus.output.pileup
#    output:
#        data('{sample}/replicate-{replicate}/{mapping_stage}/segments/{segment}/consensus-report.tsv')
#    run:
#        check_consensus_io(
#            input.fasta, input.pileup, output[0],
#            wildcards.sample, wildcards.replicate
#        )
#
#
#def full_consensus_summary_input(wildcards):
#    consensus_filepaths = []
#    replicates = metadata_dictionary[wildcards.sample]
#    for replicate in replicates.keys():
#        for segment in SEGMENTS:
#            consensus_filepaths.append(
#                data(f'{wildcards.sample}/replicate-{replicate}/remapping-{NUMBER_OF_REMAPPINGS}/segments/{segment}/consensus-report.tsv')
#            )
#    return consensus_filepaths
#
#
#rule check_sample_consensus:
#    input:
#        full_consensus_summary_input
#    output:
#        data('{sample}/consensus-report.tsv')
#    shell:
#        'csvstack {input} > {output}'

rule check_consensus_summary:
    message:
        "Comparing penultimate vs. final remapping for {wildcards.sample}, replicate {wildcards.replicate}..."
    input:
        get_consensus_summary_inputs
    output:
        data("{sample}/replicate-{replicate}/consensus_summary.tsv")
    run:
        penultimate_fasta, final_fasta = input
        compare_remappings_io(
            penultimate_fasta,
            final_fasta,
            output[0],
            wildcards.sample,
            wildcards.replicate
        )

def get_sample_consensus_summary_inputs(wildcards):
    return expand(
        data("{{sample}}/replicate-{replicate}/consensus_summary.tsv"),
        replicate=metadata_dictionary[wildcards.sample].keys()
    )

rule aggregate_sample_consensus_summary:
    message:
        "Aggregating consensus summary reports for sample {wildcards.sample}..."
    input:
        get_sample_consensus_summary_inputs
    output:
        data("{sample}/consensus_summary_report.tsv")
    run:
        aggregate_consensus_summaries_io(input, output[0])

rule aggregate_all_consensus_summary:
    message:
        "Aggregating all sample consensus summary reports into a final project summary..."
    input:
        expand(data("{sample}/consensus_summary_report.tsv"), sample=SAMPLES)
    output:
        data("consensus_summary_report.tsv")
    run:
        aggregate_consensus_summaries_io(input, output[0])

rule all_variants:
    input:
        tsv=expand(data('{sample}/ml.tsv'), sample=DUPLICATE_SAMPLES),
        html=expand(data('{sample}/ml.html'), sample=DUPLICATE_SAMPLES)
    output: data('variants.tsv')
    run:
        merge_variant_calls(input.tsv, output[0])

#rule full_consensus_summary:
#    input:
#        expand(data('{sample}/consensus-report.tsv'), sample=SAMPLES)
#    output:
#        data('consensus-report.tsv'),
#    shell:
#        'csvstack {input} > {output}'

rule all_full_segments:
    input:
        expand(
            data('{segment}.fasta'),
            segment=SEGMENTS
        )
    output:
        data('all.fasta')
    shell:
        'cat {input} > {output}'

rule all_preliminary:
    input:
        rules.full_coverage_summary.output[0]

rule all_consensus:
    input:
        #rules.full_consensus_summary.output[0],
        rules.all_full_segments.output[0]

rule all_protein:
    input:
        expand(data('{sample}/protein'), sample=SAMPLES),
        genes=rules.gene_list.output[0]
    output:
        data("protein/.done")
    shell:
        '''
        mkdir -p {ANALYSIS_DIR}/protein
        for gene in $(cat {input.genes}); do
            cat {ANALYSIS_DIR}/*/protein/$gene.fasta > {ANALYSIS_DIR}/protein/$gene.fasta
        done
        touch {output}
        '''

rule zip:
    input:
        rules.full_coverage_summary.output[0]
        #rules.full_consensus_summary.output[0]
    output:
        data('project.zip')
    shell:
        'zip -r {output} data -x "*.fastq" "*.bam" "*.sam" "*.pileup"'

def preserved_bam_input(wildcards):
    bam_filepaths = []
    for sample, replicates in metadata_dictionary.items():
        for replicate in replicates.keys():
            bam_filepaths.append(
                data(f'{sample}/replicate-{replicate}/final.bam')
            )
    return bam_filepaths

rule preserve_final_bam:
    message:
        'Preserving final BAM files for replicate {wildcards.replicate} of sample {wildcards.sample}...'
    input:
        bam=data('{sample}/replicate-{replicate}/remapping-%d/sorted.bam') % NUMBER_OF_REMAPPINGS,
        bai=data('{sample}/replicate-{replicate}/remapping-%d/sorted.bam.bai') % NUMBER_OF_REMAPPINGS
    output:
        bam=data('{sample}/replicate-{replicate}/final.bam'),
        bai=data('{sample}/replicate-{replicate}/final.bam.bai')
    shell:
        '''
        cp {input.bam} {output.bam}
        cp {input.bai} {output.bai}
        '''

rule copy_config:
    input:
        "config.yml"
    output:
        data("config.yml")
    shell:
        "cp {input} {output}"

rule all:
    input:
        rules.copy_config.output,
        rules.all_preliminary.input,
        rules.all_consensus.input,
        rules.all_protein.output,
        rules.all_variants.output,
        rules.zip.output,
        rules.aggregate_all_consensus_summary.output,
        preserved_bam_input
