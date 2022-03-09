import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
import plotly.express as px


def app():
    st.write('# Library Map')
    mfile = st.file_uploader('Upload library map file')
    mfile = '/Users/ansintsova/git_repos/mbarq_app/data/SL1344_test/library_10_1_closest.annotated.csv'
    df = pd.read_csv(mfile)
    df['library'] = 'library_10_1'
    #st.write(f"Number of unique insertions: {df.groupby(['library', 'chr', 'insertion_site']).barcode.nunique()}")
   # st.write(f"Number of CDS disrupted: ")
    fixed_col_names = ['barcode', 'number_of_reads', 'insertion_site', 'chr', 'strand',
                       'multimap', 'distance_to_feature', 'library']
    attr_names = [c for c in df.columns if c not in fixed_col_names]
    c1, c2 = st.columns(2)
    seqid = c1.selectbox('Choose sequence to display', df.chr.unique())

    df['in CDS'] = df.distance_to_feature == 0
    colorby = c2.selectbox('Color by', ['in CDS', 'library', 'multimap'])
    df = df.fillna('NaN')
    df_to_show = df[df.chr == seqid].sort_values("insertion_site")

    fig = px.scatter(df_to_show, x='insertion_site', y='number_of_reads', color= colorby, log_y=True, height=600,
                     hover_data=attr_names, labels={'insertion_site': 'Position, bp', 'number_of_reads': 'Read Counts'})
    fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                      font=dict(size=18))
    fig.update_traces(marker=dict(size=8,
                                  line=dict(width=2,
                                            color='DarkSlateGrey')),
                      selector=dict(mode='markers'))
    st.plotly_chart(fig, use_container_width=True)


    table1 = (df.groupby('library')
              .agg({'barcode': ['nunique'], 'Name': ['nunique'], 'distance_to_feature':[lambda x: sum(x != 0)],
                    'multimap': ['sum']})
              .reset_index())
    table1.columns = ["Library", '# of insertions', '# of genes with insertion',
                      '# of insertions outside of CDS', '# of barcodes mapped to multiple locations']

    table2 = (df.groupby(['library', 'Name'])
              .barcode.count().reset_index().groupby('library')
              .agg({'barcode': ['median', 'max']})
              .reset_index())
    table2.columns = ['Library', 'Median insertions per gene', 'Max insertions per gene']

    st.markdown("### Table 1: Insertion Summary")
    st.table(table1)
    st.markdown("### Table 2: Insertions per gene ")
    st.dataframe(table2)