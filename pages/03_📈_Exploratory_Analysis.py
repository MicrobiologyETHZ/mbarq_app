import streamlit as st
from scripts.datasets import CountDataSet, convert_df, define_color_scheme
import pandas as pd
from scripts.layouts import pca_layout, barcode_abundance_layout
from pathlib import Path
st.set_page_config(layout='wide')


def app():
    hide_dataframe_row_index = """
                <style>
                .row_heading.level0 {display:none}
                .blank {display:none}
                </style>
                """
    st.markdown(hide_dataframe_row_index, unsafe_allow_html=True)
    st.markdown(""" # Exploratory Analysis """)
    with st.expander("How this works: "):

        count_url = 'https://mbarq.readthedocs.io/en/latest/counting.html'
        st.markdown(f"""
        #### Count data:
        
        - For this page, you need to upload a **csv** file of merged counts produced by `mbarq count` + `mbarq merge`. For instructions on how to generate this file, please see [here]({count_url}).
        - The first column must contain the barcodes, and the second column must contain the gene identifier (ex. locus tag). 
        - All other columns must be sample names. 
        #
        Example structure:
        """)
        test = pd.DataFrame([['ACACACGT', 'abcD', '450', '700'],
                             ['GACCCCAC', 'efgH', '100', '0']], columns=['barcode', 'geneName', 'sample1', 'sample2'])
        st.table(test)
        st.markdown("""
        #### Sample data:
        - A **csv** file containing sample data. 
        - The first column must contain sample names that correspond to sample names in the count file.  
        - All other columns will be read in as metadata.
        #
        
        Example structure:
        """)
        test = pd.DataFrame([['sample1', 'treatment'], ['sample2', 'control']],
                            columns=['sampleID', 'treatment'])
        st.table(test)
        st.markdown("""

        - The merged count table produced by `mbarq` will contain barcodes found in the mapping file, as well as unannotated barcodes (e.g. control spike-ins, artifacts). Only annotated barcodes are used for the exploratory analysis.
        - Simple TSS normalization and log2 transformation are performed.
        - For PCA plot, you can choose how many barcodes are used for the analysis, as well as which components to visualize. The scree plot shows the % of variance explained by each of the PCs. 
        - For Barcode Abundance, normalized barcode counts can be visualized for any gene of interest and compared across different sample data variables. 
        """)

    # LOAD THE DATA
    with st.container():
        st.subheader('Load your own data or browse the example data set')
        if 'count_ds' in st.session_state.keys():
            cds = st.session_state['count_ds']
        else:
            st.info('Browse the example data set below or load your own data on **⬆️ Data Upload** page')
            cfile = "examples/example_mbarq_merged_counts.csv"
            mfile = "examples/example_sample_data.csv"
            c1, c2 = st.columns(2)
            c1.subheader('Example count file (sample)')
            ex_df = pd.read_csv(cfile)
            ex_sample_df = pd.read_csv(mfile)
            samples_to_show = ['dnaid1315_10', 'dnaid1315_107']
            c1.write(ex_df[['barcode', 'Name'] + samples_to_show].dropna().head())
            c2.subheader('Example metadata file (sample)')
            c2.write(pd.read_csv(mfile, index_col=0).loc[samples_to_show].reset_index())
            c1.download_button(
                label="Download example count data as CSV",
                data=convert_df(ex_df),
                file_name='example_counts_file.csv',
                mime='text/csv',
            )

            c2.download_button(
                label="Download example sample data as CSV",
                data=convert_df(ex_sample_df),
                file_name='example_sample_data_file.csv',
                mime='text/csv',
            )
            cds = CountDataSet(cfile, mfile)
        # IF DATA IS LOADED VISUALIZE
        if cds.valid:
            cds.normalize_counts()
            st.write('## PCA plot')
            # PCA GRAPH
            pca_layout(cds)
            # BARCODE ABUNDANCE
            barcode_abundance_layout(cds)
app()