import streamlit as st
from scripts.datasets import LibraryMap, CountDataSet, ResultDataSet
st.set_page_config(layout='wide')


def app():
    st.markdown("# Data Upload")
    st.markdown("## Upload library map\n"
                "Required for **Library Map** page. You can upload multiple files for comparison.")
    url_map = 'https://mbarq.readthedocs.io/en/latest/mapping.html'
    st.markdown(f'Upload a mutant [library map]({url_map}) file (ex. `library.annotated.csv`)')
    map_files = st.file_uploader(f'Upload', accept_multiple_files=True, key='map')

    if st.button('Clear', key='map_clear_button') and 'map_file' in st.session_state.keys():
        del st.session_state['map_file']
        del st.session_state['lib_map']
        del st.session_state['attributes']
        del st.session_state['annotations']

    if map_files:
        lm = LibraryMap(map_files=map_files)
        lm.load_map()
        lm.validate_lib_map()
        st.session_state['map_file'] = ", ".join([m.name for m in map_files])
        st.session_state['lib_map'] = lm
        st.session_state['attributes'] = lm.attributes
        st.session_state['annotations'] = lm.lib_map[lm.attributes].drop_duplicates()

    if 'map_file' in st.session_state.keys():
        st.write(f'**Currently loaded libary map**: {st.session_state["map_file"]}')

    st.markdown("## Upload barcode **count data**  file and **sample data** file\n"
                "Required for **Exploratory Analysis** page.")
    count_file = st.file_uploader('Upload a file containing merged count data', key='cnt')
    sample_file = st.file_uploader('Upload a file containing sample data', key='sample')

    if st.button('Clear', key='count_clear_button') and 'count_file' in st.session_state.keys():
        del st.session_state['count_file']
        del st.session_state['count_ds']

    if count_file is not None and sample_file is not None:
        cds = CountDataSet(count_file, sample_file)
        st.session_state['count_ds'] = cds
        st.session_state['count_file'] = count_file

    if 'count_file' in st.session_state.keys():
        st.write(f'**Currently loaded count file**: {st.session_state["count_file"].name}')

    st.markdown("## Upload **fitness data** file\n"
                "Required for **Differential Abundance**, **STRING**, and **KEGG** pages. "
                "You can upload multiple results files (i.e. results from different mutant libraries).")
    results_files = st.file_uploader('Upload the final results table',
                                     accept_multiple_files=True,
                                     key='res')
    if st.button('Clear', key='res_button') and 'results_files' in st.session_state.keys():
        del st.session_state['results_files']
        del st.session_state['results_ds']

    if results_files:
        gene_id = st.text_input('Unique gene identifier used in the result files', value='Name')
        rds = ResultDataSet(results_files, gene_id=gene_id)
        rds.load_results()
        rds.validate_results_df()
        if 'annotations' in st.session_state.keys():
            if rds.gene_id in st.session_state['annotations'].columns:
                rds.results_df = (rds.results_df.merge(st.session_state['annotations'],
                                                       how='left', on=rds.gene_id))
        st.session_state['results_files'] = ", ".join([m.name for m in results_files])
        st.session_state['results_ds'] = rds

    if 'results_files' in st.session_state.keys():
        st.write(f'**Currently loaded count file**: {st.session_state["results_files"]}')


if __name__ == "__main__":
    app()
