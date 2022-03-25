import streamlit as st
import pandas as pd
import plotly.express as px
from pathlib import Path
import escher
from escher import Builder
import streamlit.components.v1 as components
escher.rc['never_ask_before_quit'] = True


def app():
    st.markdown('## Metabolic Maps')

    with st.container():
        st.subheader('Browse an example data set')
        data_type = st.radio('Choose dataset to show', ['Look at an example'], index=0)
        if data_type == 'Load my data':
            result_file = st.file_uploader('Upload library map file', accept_multiple_files=True)
        else:
            result_file = Path("/Users/ansintsova/git_repos/mbarq_app/data/SL1344_test/library_10_1-unfiltered-results.kegg.csv")
            st.subheader('Example results file')
            st.write(pd.read_csv(result_file, index_col=0).sample(5))
        if not result_file:
            st.stop()

        clrs = px.colors.qualitative.D3

    with st.container():
        fdf = pd.read_csv(result_file)
        c1, c2 = st.columns(2)
        contrast_to_show = c1.selectbox("Select a contrast", fdf.contrast.sort_values().unique())
        maps = {file.stem: file for file in Path("/Users/ansintsova/git_repos/mbarq_app/data/SL1344_test/maps").glob("*json")}
        map_to_show = c2.selectbox('Select a metabolic map', [""] + list(maps.keys()))
        if not map_to_show:
            st.stop()
        ci_to_show = fdf[fdf.contrast == contrast_to_show][['Name', 'lfc']].set_index('Name')
        builder = Builder(map_json=str(maps[map_to_show]))
        builder.gene_data = ci_to_show.to_dict()['lfc']
        builder.reaction_scale = [{'type': 'min', 'color': 'red', 'size': 12},
                                  {'type': 'max', 'color': 'blue', 'size': 12},
                                  {'type': 'value', 'value': 0, 'color': 'white', 'size': 8}]
        builder.show_gene_reaction_rules = True
        builder.save_html('./data/example_map.html')
        HtmlFile = open("./data/example_map.html", 'r', encoding='utf-8')
        source_code = HtmlFile.read()
        components.html(source_code, height=800)
