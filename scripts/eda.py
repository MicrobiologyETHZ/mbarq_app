import streamlit as st
import pandas as pd
from sklearn.decomposition import PCA
import numpy as np
import plotly.express as px
from itertools import cycle
from pathlib import Path

from scripts.graphs import pca_figure, barcode_abundance, define_color_scheme, find_PCs

@st.cache
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv(index=False).encode('utf-8')


def app():
    hide_dataframe_row_index = """
                <style>
                .row_heading.level0 {display:none}
                .blank {display:none}
                </style>
                """
    st.markdown(hide_dataframe_row_index, unsafe_allow_html=True)
    st.markdown(""" # Exploratory Analysis """)
    # EXPLAIN WHAT HAPPENS ON THIS PAGE
    with st.expander("How this works: "):
        st.markdown(""" ### Visualizing barcode count data. """)
        c1, c2 = st.columns(2)
        c1.markdown("""
        - Takes in a **CSV** file of merged counts produced by `mbarq count` + `mbarq merge`. 
        - The first column must contain the barcodes, the second column must contain barcode annotation. 
        - All other columns must be sample names. 
        Example structure:
        """)
        test = pd.DataFrame([['ACACACGT', 'abcD', '450', '700'],
                             ['GACCCCAC', 'efgH', '100', '0']], columns=['barcode', 'geneName', 'sample1', 'sample2'])
        c1.table(test)
        c2.markdown("""

        - Takes in a **CSV** file containing sample data. 
        - First column must contain sample names that correspond to sample names in the count file.  
        - All other columns will be read in as metadata
        Example structure:
        """)
        test = pd.DataFrame([['sample1',  'treatment'], ['sample2',  'control']],
                        columns=['sampleID',  'treatment'])
        c2.table(test)
        st.markdown("""

        - Merged count table produced by `mbarq` will contain barcodes found in the mapping file, as well as unannotated barcodes (e.g. control spike-ins, artifacts). Only annotated barcodes are used for the exploratory analysis.
        - Simple TSS normalisation and log2 transformation is performed
        - For PCA plot, you can choose how many barcodes are used for the analysis, as well as which components to visualise. Scree plot shows the % of variance explained by each of the PCs. 
        - For Barcode Abundance, normalised barcode counts can be visualised for any gene of interest and compared across different sample data variables. 

        """)
    # LOAD THE DATA
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
        # IF DATA IS LOADED VISUALIZE
        if cfile and mfile:
            """
            Requirements: first column has barcodes, second column has attributes in the countData, rest need to be sampleIDs.
            Sample Data: first column are sampleIDs 
            Will only look at barcodes that were mapped to a feature 
            """
            df = pd.read_csv(cfile)
            countData = (df.rename({df.columns[0]: 'barcode', df.columns[1]: 'Name'}, axis=1)
                         .dropna()  # todo make sure that does the correct thing
                         .drop_duplicates()
                         .set_index(['barcode', 'Name']))
            countData = np.log2((countData / countData.sum()) * 1000000 + 0.5).reset_index()
            sampleData = pd.read_csv(mfile).fillna('N/A')
            sampleData = sampleData.rename({sampleData.columns[0]: 'sampleID'})
            to_keep = list(set(sampleData.sampleID.unique()).intersection(countData.columns))
            # CHECK THAT METADATA MATCHES COLUMNS IN THE COUNT FILE
            if len(to_keep) == 0:
                st.write("Sample IDs do not match any columns in the count file")
                st.stop()
            pcaDf = countData[['barcode'] + to_keep].copy()
            sampleData = sampleData[sampleData.sampleID.isin(to_keep)]
            # Stay consistent -> no index
            st.write('## PCA plot')
            # PCA GRAPH
            with st.expander('Show PCA'):
                _, aC, all_clrs = define_color_scheme()
                c1, c2 = st.columns((3, 1))
                c2.write('### PCA Options')
                numPCs = c2.number_input("Select number of Principal Components", min_value=2, max_value=50, value=10)
                numGenes = c2.number_input("Number of genes to use", min_value=numPCs, value=min(250, pcaDf.shape[0]),
                                           max_value=pcaDf.shape[0])
                chooseBy = 'variance'
                numGenes = int(numGenes)
                numPCs = int(numPCs)
                pcDf, pcVar = find_PCs(pcaDf, sampleData, numPCs, numGenes, chooseBy)
                missingMeta = " ,".join(list(pcDf[pcDf.isna().any(axis=1)].index))
                if missingMeta:
                    st.write(f"The following samples have missing_metadata and will not be shown: {missingMeta}")
                pcDf = pcDf[~pcDf.isna().any(axis=1)]
                pcxLabels = [f'PC{i}' for i in range(1, numPCs + 1)]
                expVars = [c for c in pcDf.columns if c not in pcxLabels]
                pcX = c2.selectbox('X-axis component', pcxLabels)
                pcY = c2.selectbox('Y-axis component', [pc for pc in pcxLabels if pc != pcX])
                pcVarHi = c2.radio('Variable to highlight', expVars)
                pcSym = c2.radio('Variable to show as symbol', [None] + expVars)
                fig1, fig2, fig3 = pca_figure(pcDf, pcX, pcY, pcVarHi, pcVar, pcSym, expVars, all_clrs)
                c1.write(f'### {pcX} vs {pcY}, highlighting {pcVarHi}')
                c1.plotly_chart(fig1, use_container_width=True)
                c3, c4 = st.columns(2)
                c3.write('### Scree Plot')
                c3.plotly_chart(fig2)
                c4.write(f'### PCs summarized by {pcVarHi}')
                c4.plotly_chart(fig3, use_container_width=True)

            # BARCODE ABUNDANCE
            st.write('## Barcode Abundance')
            with st.expander('Show Barcode Abundance'):
                # Process the dataframe
                abDf = countData.dropna()
                barcode = abDf.columns[0]
                gene_name = abDf.columns[1]
                # Get user input
                c1, c2 = st.columns(2)
                compareCondition = c1.selectbox('Which conditions to compare?', sampleData.columns)
                conditionCategories = c1.multiselect(f'Categories of {compareCondition} to display',
                                            ['All'] + list(sampleData[compareCondition].unique()))

                filterCondition = c2.selectbox("Filter by", ['No filter'] + list(sampleData.columns))
                if filterCondition == 'No filter':
                    filterCategories = []
                else:
                    filterCategories = c2.multiselect(f'Which category(ies) of {filterCondition} to keep?',
                                                      list(sampleData[filterCondition].unique()))
                if 'All' in conditionCategories:
                    conditionCategories = list(sampleData[compareCondition].unique())
                genes = st.multiselect("Choose gene(s) of interest", abDf[gene_name].unique())
                if len(genes) * len(conditionCategories) > 40:
                    st.write('Too many genes/categories to display, consider choosing fewer genes')
                else:
                    geneDf = abDf[abDf[gene_name].isin(genes)]
                    abSampleDf = sampleData[sampleData[compareCondition].isin(conditionCategories)]
                    if filterCategories:
                        abSampleDf = abSampleDf[abSampleDf[filterCondition].isin(filterCategories)]
                    geneDf = (geneDf.melt(id_vars=[barcode, gene_name], value_name='log2CPM', var_name='sampleID')
                                    .merge(abSampleDf, how='inner', on='sampleID')
                                    .sort_values(compareCondition))
                    groupBy = st.radio('Group by', [gene_name, compareCondition])
                    colorBy = [c for c in [gene_name, compareCondition] if c != groupBy][0]
                    fig = barcode_abundance(geneDf, groupBy, colorBy, all_clrs)
                    st.plotly_chart(fig, use_container_width=True)



