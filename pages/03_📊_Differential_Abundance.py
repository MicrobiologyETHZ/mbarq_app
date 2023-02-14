import streamlit as st
from scripts.datasets import ResultDataSet
from pathlib import Path


def app():
    st.markdown(""" # Fitness Results """)
    with st.expander('How this works: '):
        st.markdown("""TBC""")

    with st.container():
        # Get the data
        if 'results_ds' in st.session_state.keys():
            rds = st.session_state['results_ds']
            gene_id = st.selectbox['results_gene_id']
        else:
            st.info('Browse example results file or upload your data in **Data Upload**')
            result_files = [Path("examples/example_rra_results.csv")]
            gene_id = 'Name'
            rds = ResultDataSet(result_files=result_files, gene_id=gene_id)
            rds.load_results()
            rds.validate_results_df()

        if st.checkbox('Show sample of the dataset'):
            try:
                st.write(rds.results_df.sample(5))
            except ValueError:
                st.write('Result table is empty')

    if not rds.results_df.empty:
        contrasts = rds.results_df[rds.contrast_col].sort_values().unique()
        libraries = rds.results_df[rds.library_col].sort_values().unique()
        libraries = ['All'] + list(libraries) if len(libraries) > 1 else libraries

        st.subheader("Fitness Results")
        # SUBSET DATAFRAME TO SPECIFIC CONTRAST AND LIBRARY
        contrast_col, lfc_col, fdr_col, lfc_lib_col = st.columns(4)
        contrast_to_show = contrast_col.multiselect('Select a contrast', contrasts, default=contrasts[0])
        library_to_show = lfc_lib_col.selectbox('Select library to show', libraries)
        fdr_th = fdr_col.number_input('FDR cutoff', value=0.05)
        type_lfc_th = lfc_col.radio('Absolute LFC cutoff or define range', ['Absolute', 'Range'])
        if type_lfc_th == 'Absolute':
            lfc_low = lfc_col.number_input('Log FC cutoff (absolute)', min_value=0.0, step=0.5, value=1.0)
            lfc_hi = None
        else:
            lfc_low = lfc_col.number_input('Min Log FC',  step=0.5, value=-5.0)
            lfc_hi = lfc_col.number_input('Max Log FC',  step=0.5, value=-1.0)

        rds.identify_hits(library_to_show,  lfc_low, lfc_hi, fdr_th)
        st.write(rds.results_df.head())


app()
