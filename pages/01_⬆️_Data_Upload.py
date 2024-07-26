import streamlit as st
from scripts.datasets import LibraryMap, CountDataSet, ResultDataSet
st.set_page_config(layout='wide')

def app():

    st.markdown("# Data Upload")
    info_url = "https://docs.streamlit.io/knowledge-base/using-streamlit/where-file-uploader-store-when-deleted"
    st.info(f""" mBARq app does not store your data. You can read more on this [here]({info_url}). """)

    # Uploading mapping data
    st.markdown("## Upload library map\n"
                "Required for **Library Map** page. You can upload multiple files for comparison.")
    url_map = 'https://mbarq.readthedocs.io/en/latest/mapping.html'
    upload_msg = f'Upload a RB-TnSeq [library map]({url_map}) file (ex. `library.annotated.csv`)'
    with st.form("map-form", clear_on_submit=True):
        map_files = st.file_uploader(upload_msg, accept_multiple_files=True, key='map_file_key')
        submitted = st.form_submit_button("Submit")

    if map_files:
        lm = LibraryMap(map_files=map_files)
        lm.load_map()
        lm.validate_lib_map()
        if not lm.lib_map.empty:
            st.session_state['lib_map'] = lm
            st.session_state['annotations'] = lm.lib_map[lm.attributes].drop_duplicates()

    # If clear button is pressed removed the file and lib_map from session state
    if st.button('Clear loaded map', key='map_clear_button') and 'lib_map' in st.session_state.keys():
        del st.session_state['lib_map']
        del st.session_state['annotations']

    if 'lib_map' in st.session_state.keys():
        st.info(f'**Currently loaded libary map**: {", ".join([m.name for m in st.session_state["lib_map"].map_files])}')
    else:
        st.info("**No library map is currently loaded**")

    # Uploading count data
    count_url = "https://mbarq.readthedocs.io/en/latest/counting.html"
    st.markdown("## Upload barcode **count data**  file and **sample data** file\n"
                f"Learn more about file format [here]({count_url}). "
                "Required for **Exploratory Analysis** page.")
    with st.form("count-form", clear_on_submit=True):
        count_file = st.file_uploader('Upload a file containing merged count data', key='count_file_key')
        sample_file = st.file_uploader('Upload a file containing sample data', key='sample_data_key')
        st.form_submit_button("Submit")

    if count_file is not None and sample_file is not None:
        cds = CountDataSet(count_file, sample_file)
        # todo add validation step?
        st.session_state['count_ds'] = cds

    if st.button('Clear loaded count and sample data', key='count_clear_button') and 'count_ds' in st.session_state.keys():
        del st.session_state['count_ds']
    if 'count_ds' in st.session_state.keys():
        st.info(f'**Currently loaded  files**: `{st.session_state["count_ds"].count_file.name}` '
                f'and `{st.session_state["count_ds"].sample_data_file.name}`')
    else:
        st.info("**No count/sample files are currently loaded**")

    # Upload the analysis results
    analysis_url = "https://mbarq.readthedocs.io/en/latest/analysis.html"
    st.markdown("## Upload **fitness data** file\n"
                "Required for **Differential Abundance**, **STRING**, and **KEGG** pages. "
                f"Learn more about file format [here]({analysis_url}). "
                "You can upload multiple results files (i.e. results from different mutant libraries).")
    with st.form("fitness-form", clear_on_submit=True):
        results_files = st.file_uploader('Upload the final results table',
                                     accept_multiple_files=True,
                                     key='results_file_key')
        st.form_submit_button("Submit")

    if results_files:
        gene_id = st.text_input('Unique gene identifier used in the result files', value='Name')
        rds = ResultDataSet(results_files, gene_id=gene_id)
        rds.load_results()
        rds.validate_results_df()
        if 'annotations' in st.session_state.keys():
            if rds.gene_id in st.session_state['annotations'].columns:
                rds.results_df = (rds.results_df.merge(st.session_state['annotations'],
                                                       how='left', on=rds.gene_id))
        st.session_state['results_ds'] = rds

    if st.button('Clear loaded fitness data', key='res_button') and 'results_ds' in st.session_state.keys():
        del st.session_state['results_ds']

    if 'results_ds' in st.session_state.keys():
        st.info(f'**Currently loaded count file**: `{", ".join([m.name for m in st.session_state["results_ds"].result_files])}`')
    else:
        st.info("**No result file is currently loaded**")


if __name__ == "__main__":
    app()
