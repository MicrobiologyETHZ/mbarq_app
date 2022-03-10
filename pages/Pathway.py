import escher
from escher import Builder
import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
import plotly.express as px
from pathlib import Path
import escher
from escher import Builder
import cobra
import streamlit.components.v1 as components
escher.rc['never_ask_before_quit'] = True


def app():
    st.markdown('## Metabolic Maps')
    clrs = px.colors.qualitative.D3

    rfile = st.file_uploader('Load the analysis file')
    rfile = "/Users/ansintsova/git_repos/mbarq_app/data/SL1344_test/summarized_mageck_results.csv"
    fdf = pd.read_csv(rfile)

    c1, c2 = st.columns(2)
    contrast_to_show = c1.selectbox("Select a contrast", fdf.contrast.sort_values().unique())
    maps = {file.stem: file for file in Path("/Users/ansintsova/git_repos/mbarq_app/data/SL1344_test/maps").glob("*json")}
    map_to_show = c2.selectbox('Select a metabolic map', [""] + list(maps.keys()))
    if not map_to_show:
        st.stop()
    ci_to_show = fdf[fdf.contrast == contrast_to_show][['gene', 'logCI']].set_index('gene')
    builder = Builder(map_json=str(maps[map_to_show]))
    builder.gene_data = ci_to_show
    builder.reaction_scale = [{'type': 'min', 'color': 'red', 'size': 12},
                              {'type': 'max', 'color': 'blue', 'size': 12},
                              {'type': 'value', 'value': 0, 'color': 'white', 'size': 8}]
    builder.show_gene_reaction_rules = True
    builder.save_html('./data/example_map.html')
    HtmlFile = open("./data/example_map.html", 'r', encoding='utf-8')
    source_code = HtmlFile.read()
    components.html(source_code, height=800)
