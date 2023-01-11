import streamlit as st
import pandas as pd
from pathlib import Path
import plotly.express as px
from scripts.graphs import graph_library_map

def load_library_map(uploaded_map):  # todo refactor so can be cached; don't think can do this...
    df = pd.read_csv(uploaded_map)
    df['in CDS'] = df.distance_to_feature == 0
    df = df.fillna('NaN')  # todo check if this does anything, shouldn't
    return df


def process_library_maps(map_files):
    map_dfs = []
    attr_names = []
    fixed_col_names = ['barcode', 'abundance_in_mapping_library', 'insertion_site',
                       'chr', 'distance_to_feature', 'strand', 'library']   # todo make sure this is up to date
    for uploaded_map in map_files:
        st.write(f"Processing {uploaded_map.name}")
        df = load_library_map(uploaded_map)
        if 'library' not in df.columns:
            library_name = st.text_input("Change library name?", value=uploaded_map.name)
            df['library'] = library_name
        if 'number_of_reads' in df.columns:
            df = df.rename({'number_of_reads': 'abundance_in_mapping_library'}, axis=1)
        missing_cols = [c for c in fixed_col_names if c not in df.columns]
        if len(missing_cols) > 0:
            st.markdown(
                f"""⚠️ The following columns are missing from the map files: {', '.join(missing_cols)}. 
                Please rename the columns/rerun mBARq and try again. Skipping {uploaded_map.name} ⚠️""")
        attrs = [c for c in df.columns if c not in fixed_col_names]
        map_dfs.append(df)
        attr_names.extend(attrs)
    return pd.concat(map_dfs), attr_names


@st.cache
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv(index=False).encode('utf-8')


def app():
    st.markdown(""" # Library Map """)

    with st.expander('How this works: '):
        st.markdown("""
    
        ### Visualize insertion position along the genome.
    
        - Takes in library map **CSV** file (`*.annotated.csv`) produced by `mbarq map`. Expects to find the following columns in the file: `barcode`, `abundance_in_mapping_library`, `insertion_site`, `chr`, `distance_to_feature`, `strand`. 
        - Can load more than one library file at the same time to compare! 
        - You can select which sequence (e.g. chromosome or plasmids) to display, and color the insertions by distance to feature, multimapping or library.
        - You can click on the figure legend to only show a specific subset of data (i.e. if looking at multiple libraries, double clicking on the specific library name will show data for that library only)
    
        """)

    with st.container():
        st.subheader('Load your own data or browse the example data set')
        data_type = st.radio('Choose dataset to show', ['Look at an example', 'Load my data'], index=1, key='lib')
        if data_type == 'Load my data':
            map_files = st.file_uploader('Upload library map file', accept_multiple_files=True)
        else:
            map_files = [Path('examples/example_library_mapping_file.annotated.csv')]
            st.subheader('Example mapping file')  # todo update example mapping file not to include multimap column
            example_df = pd.read_csv(map_files[0])
            st.write(example_df.sample(5))
            st.download_button(
                label="Download example data as CSV",
                data=convert_df(example_df),
                file_name='example_library_mapping_file.csv',
                mime='text/csv',
            )

    if len(map_files) > 0:
        map_df, attr_names = process_library_maps(map_files)
        with st.container():
            graph_library_map(map_df, attr_names)

        with st.container():
            name_col = st.selectbox("Choose attribute column", attr_names)
            table1 = (map_df.groupby('library')
                      .agg({'barcode': ['nunique'], 'distance_to_feature': [lambda x: sum(x != 0)]})
                      .reset_index())
            table1.columns = ["Library", '# of insertions', '# of insertions outside of CDS']
            table1 = table1.set_index('Library')
            table1[f'# of {name_col}s with insertion'] = (map_df[map_df.distance_to_feature == 0]
                                                  .groupby('library')[name_col].nunique())
            # # table1['Library'] = table1.Library.str.replace("_", '-')
            table2 = (map_df.groupby(['library', name_col])
                      .barcode.count().reset_index().groupby('library')
                      .agg({'barcode': ['median', 'max']})
                      .reset_index())
            table2.columns = ['Library', 'Median insertions per gene', 'Max insertions per gene']
            table2['Median insertions per gene'] = table2["Median insertions per gene"].astype(int)
            table2 = table2.set_index('Library')

            st.markdown("#### Insertion Summary")
            st.table(table1.merge(table2, left_index=True, right_index=True))




