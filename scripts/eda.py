import streamlit as st
import pandas as pd
from sklearn.decomposition import PCA
import numpy as np
import plotly.express as px
from itertools import cycle
from pathlib import Path

from scripts.graphs import pca_figure, barcode_abundance

@st.cache
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv(index=False).encode('utf-8')


def app():

    st.markdown(""" # Exploratory Analysis """)

    with st.expander("How this works: "):

        st.markdown(""" ### Visualizing barcode count data. """)

        st.markdown("""
        - Takes in a **CSV** file of merged counts produced by `mbarq count` + `mbarq merge`. The first column must contain the barcodes, the second column must contain barcode annotation (e.g. gene name or locus tag). All other columns must be sample names. Example structure:
        #
        """)

        test = pd.DataFrame([['ACACACGT', 'abcD', '450', '700']], columns=['barcode', 'geneName', 'sample1', 'sample2'])
        st.dataframe(test)

        st.markdown("""

        - Takes in a **CSV** file containing sample data. First column must contain sample names that correspond to sample names in the count file. Example structure:
        #
        """)
        test = pd.DataFrame([['sample1', 'batch1', 'treatment'], ['sample2', 'batch2', 'control']],
                        columns=['sampleID', 'batch', 'treatment'])
        st.dataframe(test)
        st.markdown("""

        - Merged count table produced by `mbarq` will contain barcodes found in the mapping file, as well as unannotated barcodes (e.g. control spike-ins, artifacts). Only annotated barcodes are used for the exploratory analysis.

        - Simple TSS normalisation and log2 transformation is performed
        - For PCA plot, you can choose how many barcodes are used for the analysis, as well as which components to visualise. Scree plot shows the % of variance explained by each of the PCs. 
        - For Barcode Abundance, normalised barcode counts can be visualised for any gene of interest and compared across different sample data variables. 

        """)

    with st.container():
        st.subheader('Load your own data or browse the example data set')
        data_type = st.radio('Choose dataset to show', ['Look at an example', 'Load my data'], index=1, key='exp')
        if data_type == 'Load my data':
            cfile = st.file_uploader('Load merged count file')
            mfile = st.file_uploader('Load metadata')
        else:
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

        if cfile and mfile:
            # Requirements: first column has barcodes, second column has attributes in the countData, rest need to be sampleIDs.
            # first column are sampleIDs column in the sampleData
            # Will only look at barcodes that were mapped to a feature / alternatively filter out low counts?
            df = pd.read_csv(cfile)
            countData = (df.rename({df.columns[0]: 'barcode'}, axis=1)
                         .dropna()
                         .drop([df.columns[1]], axis=1)
                         .drop_duplicates()
                         .set_index('barcode'))
            sampleData = pd.read_csv(mfile, index_col=0).fillna('N/A')
            sampleData.index.name = 'sampleID'
            to_keep = list(sampleData.index)
            countData = countData[to_keep]
            to_log = st.checkbox('Apply log2 transform (recommended)', value=True)
            if to_log:
                countData = np.log2((countData / countData.sum()) * 1000000 + 0.5)
            else:
                countData = countData / countData.sum() * 1000000
            st.write('## PCA plot')
            with st.expander('Show PCA'):
                pca_figure(countData, sampleData)

            st.write('## Barcode Abundance')
            with st.expander('Show Barcode Abundance'):
                barcode_abundance(df, sampleData)




