import streamlit as st
import pandas as pd
from pathlib import Path
import plotly.express as px


def process_library_map(uploaded_map):
    df = pd.read_csv(uploaded_map)
    if 'library' not in df.columns:
        library_name = st.text_input("Please provide library name", value=uploaded_map.name)
        df['library'] = library_name
    fixed_col_names = ['barcode', 'number_of_reads', 'insertion_site', 'chr',
                       'multimap', 'distance_to_feature', 'strand', 'library']
    missing_cols = [c for c in fixed_col_names if c not in df.columns]
    attr_names = [c for c in df.columns if c not in fixed_col_names]
    if len(missing_cols) > 0:
        st.markdown(
            f"""The following columns are missing from the map files: {', '.join(missing_cols)}. 
            Please rename the columns/rerun mBARq and try again. Skipping {uploaded_map.name}""")
        return pd.DataFrame()
    df['in CDS'] = df.distance_to_feature == 0
    df = df.fillna('NaN')  # todo check if this does anything, shouldn't
    return df, attr_names


def app():
    st.markdown(""" # Library Map
    
    Write a little explanation of what this page shows.
    
    - Required Inputs
    - How to use the graph
    
    """)

    with st.expander("Upload library map file"):
        map_files = st.file_uploader('Upload library map file', accept_multiple_files=True)
        # map_files = [Path('/Users/ansintsova/git_repos/mbarq_app/data/SL1344_test/library_10_1_closest.annotated.csv'),
        #              Path('/Users/ansintsova/git_repos/mbarq_app/data/SL1344_test/library_10_1_copy_closest.annotated.csv')]

        if len(map_files) < 1:
            st.stop()
        map_dfs = []
        attr_names = []
        for uploaded_map in map_files:
            st.write(f"Processing {uploaded_map.name}")
            df, attrs = process_library_map(uploaded_map)
            map_dfs.append(df)
            attr_names.extend(attrs)
        map_df = pd.concat(map_dfs)

    with st.container():
        chr_col = 'chr'
        color_by_cols = ['in CDS', 'multimap', 'library']
        attr_names = set(attr_names)
        c1, c2 = st.columns(2)
        seqid = c1.selectbox('Choose sequence to display', map_df[chr_col].unique())
        colorby = c2.selectbox('Color by', color_by_cols)
        df_to_show = map_df[map_df[chr_col] == seqid].sort_values("insertion_site")
        fig = px.scatter(df_to_show, x='insertion_site', y='number_of_reads', color=colorby, log_y=True, height=600,
                         hover_data=attr_names,
                         labels={'insertion_site': 'Position, bp', 'number_of_reads': 'Read Counts'})
        fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                          font=dict(size=18))
        fig.update_traces(marker=dict(size=6,
                                      line=dict(width=2,
                                                color='DarkSlateGrey')),
                          selector=dict(mode='markers'))
        st.plotly_chart(fig, use_container_width=True)

    with st.container():
        name_col = st.selectbox("Choose gene name column", attr_names)
        table1 = (map_df.groupby('library')
                  .agg({'barcode': ['nunique'], name_col: ['nunique'], 'distance_to_feature': [lambda x: sum(x != 0)],
                        'multimap': ['sum']})
                  .reset_index())
        table1.columns = ["Library", '# of insertions', '# of genes with insertion',
                          '# of insertions outside of CDS', '# of barcodes mapped to multiple locations']
        table2 = (map_df.groupby(['library', name_col])
                  .barcode.count().reset_index().groupby('library')
                  .agg({'barcode': ['median', 'max']})
                  .reset_index())
        table2.columns = ['Library', 'Median insertions per gene', 'Max insertions per gene']
        st.markdown("### Table 1: Insertion Summary")
        st.table(table1)
        st.markdown("### Table 2: Insertions per gene ")
        st.dataframe(table2)
