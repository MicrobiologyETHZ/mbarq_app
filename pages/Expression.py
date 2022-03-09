import escher
from escher import Builder
import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
import plotly.express as px
from itertools import cycle


def app():
    st.write('# Barcode Abundance')
    clrs = px.colors.qualitative.Plotly
    cfile = st.file_uploader('Load merged count file')
    mfile = st.file_uploader('Load metadata')

    cfile = "/Users/ansintsova/git_repos/mbarq_app/data/SL1344_test/library_10_1_count_batchcorrected.txt"
    mfile = "/Users/ansintsova/git_repos/mbarq_app/data/SL1344_test/library_10_1_sampleData.csv"


    df = pd.read_table(cfile).set_index(['gene', 'barcode'])
    df = np.log2(df/df.sum()*1000000 +0.5)
    df = df.reset_index()
    sampleData = pd.read_csv(mfile)
    c1, c2 = st.columns(2)
    compare_by = c1.selectbox('Compare by', sampleData.columns)
    color_by = c2.selectbox('Color by',  ['barcode']+ list(sampleData.columns) )
    genes = st.multiselect("Choose gene(s) of interest", df.gene.unique())

    if not genes:
        st.stop()
    c3, c4 = st.columns(2)
    for col, gene in zip(cycle([c3, c4]), genes):
       gene_df = df[df.gene == gene]
       gene_df = (gene_df.melt(id_vars=['barcode', 'gene'], value_name='log2CPM', var_name='sampleID')
                  .merge(sampleData, how='left', on='sampleID'))
       fig = px.strip(gene_df, title = gene, x=compare_by, y='log2CPM', color=color_by,
                        hover_data=['barcode'] + list(sampleData.columns), stripmode='overlay')
       fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                         font=dict(size=16))
       fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='LightGrey')
       col.plotly_chart(fig, use_container_width=True)