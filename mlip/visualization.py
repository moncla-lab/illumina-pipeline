import os
import http.server
import socketserver
import webbrowser
import socket

import pandas as pd
import altair as alt


def replicate_variant_plot(input_tsv_filepath, output_html_filepath):
    df = pd.read_csv(input_tsv_filepath, sep="\t")
    chart = (
        alt.Chart(df, width=400, height=400)
        .mark_point(size=100, filled=True, opacity=0.6)
        .encode(
            x=alt.X(
                "Frequency_1:Q",
                title="Sample 1 variant frequency",
                scale=alt.Scale(domain=[0, 0.5]),
            ),
            y=alt.Y(
                "Frequency_2:Q",
                title="Sample 2 variant frequency",
                scale=alt.Scale(domain=[0, 0.5]),
            ),
            color=alt.Color("segment:N", title="Segment"),
            tooltip=[
                alt.Tooltip("segment:N", title="Segment"),
                alt.Tooltip("position:Q", title="Position"),
                alt.Tooltip("allele:N", title="Allele"),
                alt.Tooltip("coding_region_change:N", title="Coding region change"),
                alt.Tooltip("gene:N", title="Gene"),
                alt.Tooltip("Frequency_1:Q", title="Frequency 1"),
                alt.Tooltip("Frequency_2:Q", title="Frequency 2"),
            ],
        )
        .interactive()
        .properties(
            title=alt.TitleParams(
                text="Replicate variant calls", fontSize=24, anchor="middle"
            ),
            padding={"left": 20, "right": 20, "top": 10, "bottom": 10},
        )
        .configure_axis(labelFontSize=20, titleFontSize=20)
        .configure_legend(titleFontSize=20, labelFontSize=20)
    )
    chart.save(output_html_filepath)


def find_free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("", 0))
        return s.getsockname()[1]


def start_server():
    os.chdir("data")
    PORT = find_free_port()
    handler = http.server.SimpleHTTPRequestHandler

    with socketserver.TCPServer(("", PORT), handler) as httpd:
        print(f"Serving ./data on port {PORT}")
        webbrowser.open(f"http://localhost:{PORT}")
        httpd.serve_forever()


if __name__ == "__main__":
    start_server()
