import streamlit as st
import pandas as pd
from scripts.datasets import LibraryMap, ResultDataSet

def app():
    st.markdown("## Upload library map")
    map_files = st.file_uploader('Upload a mutant library map file',
                                 accept_multiple_files=True, key='map')

    if st.button('Clear', key='map_clear_button') and 'map_file' in st.session_state.keys():
        del st.session_state['map_file']
        del st.session_state['lib_map']

    if map_files:
        lm = LibraryMap(map_files=map_files)
        lm.load_map()
        lm.validate_lib_map()
        st.session_state['map_file'] = ", ".join([m.name for m in map_files])
        st.session_state['lib_map'] = lm

    if 'map_file' in st.session_state.keys():
        st.write(f'Currently loaded libary map: {st.session_state["map_file"]}')

    st.markdown("## Upload count file and sample data file")
    count_file = st.file_uploader('Upload a file containing earthquake data', key='cnt')
    sample_file = st.file_uploader('Upload a file containing earthquake data', key='sample')

    if count_file is not None:
        count_df = pd.read_csv(count_file)
        st.session_state['count_file'] = count_file
        st.session_state['count_df'] = count_df

    if sample_file is not None:
        sample_df = pd.read_csv(sample_file)
        st.session_state['sample_file'] = sample_file
        st.session_state['sample_df'] = sample_df

    if 'count_file' in st.session_state.keys():
        st.write(f'Currently loaded count file: {st.session_state["count_file"].name}')

    st.markdown("## Upload results file")
    results_files = st.file_uploader('Upload a file final results table',
                                    accept_multiple_files=True,
                                    key='res')
    gene_id = st.text_input('Unique gene identifier used in the result files', value='Name')

    if st.button('Clear', key='res_button') and 'results_files' in st.session_state.keys():
        del st.session_state['results_files']
        del st.session_state['results_ds']
        del st.session_state['results_gene_id']

    if results_files:
        rds = ResultDataSet(results_files, gene_id)
        rds.load_results()
        rds.validate_results_df()
        st.session_state['results_files'] = ", ".join([m.name for m in results_files])
        st.session_state['results_ds'] = rds
        st.session_state['results_gene_id'] = gene_id

    if 'results_file' in st.session_state.keys():
        st.write(f'Currently loaded count file: {st.session_state["results_files"].name}')


app()
