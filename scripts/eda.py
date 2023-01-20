import streamlit as st
import pandas as pd
from sklearn.decomposition import PCA
import numpy as np
import pandera as pa
from pandera.typing import Index, DataFrame, Series
from pandera.errors import SchemaError
from scripts.graphs import pca_figure, barcode_abundance_box, barcode_abundance_violin, define_color_scheme
import yaml

with open('scripts/config.yaml', 'r') as cf:
    config = yaml.load(cf, Loader=yaml.SafeLoader)['eda']

# Load column naming schema
col_name_config = config['fixed_column_names']
FIXED_COLUMN_NAMES = list(col_name_config.values())
BARCODE_COL = col_name_config['barcode_col']
GENENAME_COL = col_name_config['genename_col']
SAMPLEID_COL = col_name_config['sampleID_col']
NAME_COL = col_name_config['name_col']

@st.cache
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv(index=False).encode('utf-8')


class CountsSchema(pa.SchemaModel):
    barcode: Series[str] = pa.Field(coerce=True)


class CountDataSet:
    def __init__(self, countFile, sampleDataFile):
        self.countData = pd.read_csv(countFile)
        self.sampleData = pd.read_csv(sampleDataFile).fillna('N/A')
        self.valid = self._validate()

    def _validate(self):
        """
        First column of sampleData should be sampleIDs
        """
        self.sampleData = self.sampleData.rename({self.sampleData.columns[0]: 'sampleID'})
        samplesFound = list(set(self.sampleData.sampleID.unique()).intersection(self.countData.columns))
        if not samplesFound or any([x in samplesFound for x in [self.countData.columns[0], self.countData.columns[1]]]):
            return False
        self.sampleData = self.sampleData[self.sampleData.sampleID.isin(samplesFound)]
        self.countData = (self.countData.rename({self.countData.columns[0]: 'barcode',
                                                 self.countData.columns[1]: 'Gene Name'}, axis=1)
                          .dropna(subset=['Gene Name'])
                          .drop_duplicates())
        self.countData = self.countData[['barcode', 'Gene Name'] + samplesFound]
        if self.countData.empty:
            return False
        return True

    def normalize_counts(self):
        self.countData = self.countData.set_index(['barcode', 'Gene Name'])
        self.countData = self.countData.loc[:, self.countData.sum() > 0]
        self.countData = np.log2((self.countData / self.countData.sum()) * 1000000 + 0.5).reset_index()

    def get_principal_components(self, numPCs, numGenes, chooseBy):
        """
        :param numPCs:
        :param numGenes:
        :param chooseBy:
        :return:
        """
        pcaDf = self.countData.set_index('barcode').copy()
        pcaDf = pcaDf.drop('Gene Name', axis=1)
        pcaSd = self.sampleData.set_index('sampleID').apply(lambda x: x.astype('category'))
        if numGenes:
            # calculate var for each, pick numGenes top var across samples -> df
            if chooseBy == 'variance':
                genes = pcaDf.var(axis=1).sort_values(ascending=False).head(int(numGenes)).index
                pcaDf = pcaDf.loc[genes].T
            else:
                pass
                # todo implement log2fc selection
        else:
            pcaDf = pcaDf.T
        pca = PCA(n_components=numPCs)
        principalComponents = pca.fit_transform(pcaDf)
        pcs = [f'PC{i}' for i in range(1, numPCs + 1)]
        pcDf = (pd.DataFrame(data=principalComponents, columns=pcs)
                .set_index(pcaDf.index))
        pcVar = {pcs[i]: round(pca.explained_variance_ratio_[i] * 100, 2) for i in range(0, numPCs)}
        pcDf = pcDf.merge(pcaSd, how="left", left_index=True, right_index=True)
        return pcDf, pcVar


#########
#  APP  #
#########


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
        test = pd.DataFrame([['sample1', 'treatment'], ['sample2', 'control']],
                            columns=['sampleID', 'treatment'])
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
            # """
            cds = CountDataSet(cfile, mfile)
            if not cds.valid:
                st.write("Sample IDs do not match any columns in the count file")
                st.stop()
            cds.normalize_counts()
            st.write('## PCA plot')
            # PCA GRAPH
            with st.expander('Show PCA'):
                _, aC, all_clrs = define_color_scheme()
                st.write('### PCA Options')
                c1, c2, c3, c4 = st.columns(4)
                numPCs = c1.number_input("Select number of Principal Components", min_value=2, max_value=50, value=10)
                numGenes = c2.number_input("Number of genes to use", min_value=int(numPCs),
                                           value=int(min(250, cds.countData.shape[0])),
                                           max_value=int(cds.countData.shape[0]))
                chooseBy = 'variance'
                numGenes = int(numGenes)
                numPCs = int(numPCs)
                # pcDf, pcVar = find_PCs(pcaDf, sampleData, numPCs, numGenes, chooseBy)
                pcDf, pcVar = cds.get_principal_components(numPCs, numGenes, chooseBy)
                missingMeta = " ,".join(list(pcDf[pcDf.isna().any(axis=1)].index))
                if missingMeta:
                    st.write(f"The following samples have missing_metadata and will not be shown: {missingMeta}")
                pcDf = pcDf[~pcDf.isna().any(axis=1)]  # todo this should be included in the function
                pcxLabels = [f'PC{i}' for i in range(1, numPCs + 1)]
                expVars = [c for c in pcDf.columns if c not in pcxLabels]
                pcX = c1.selectbox('X-axis component', pcxLabels)
                pcY = c2.selectbox('Y-axis component', [pc for pc in pcxLabels if pc != pcX])
                pcVarHi = c3.radio('Variable to highlight', expVars)
                pcSym = c4.radio('Variable to show as symbol', [None] + expVars)
                pcDf = pcDf.sort_values(pcVarHi)
                fig1, fig2, fig3 = pca_figure(pcDf, pcX, pcY, pcVarHi, pcVar, pcSym, expVars, all_clrs)
                c1.write(f'### {pcX} vs {pcY}, highlighting {pcVarHi}')
                st.plotly_chart(fig1, use_container_width=True)
                c5, c6 = st.columns(2)
                c5.write('### Scree Plot')
                c5.plotly_chart(fig2)
                c6.write(f'### PCs summarized by {pcVarHi}')
                c6.plotly_chart(fig3, use_container_width=True)

            # BARCODE ABUNDANCE
            st.write('## Barcode Abundance')
            with st.expander('Show Barcode Abundance'):
                # Process the dataframe
                # abDf = countData.dropna()
                barcode = cds.countData.columns[0]
                gene_name = cds.countData.columns[1]

                # Get user input
                c1, c2 = st.columns(2)
                compare_condition = c1.selectbox('Which conditions to compare?', cds.sampleData.columns)
                condition_categories = c1.multiselect(f'Categories of {compare_condition} to display',
                                                      ['All'] + list(cds.sampleData[compare_condition].unique()))
                filter_condition = c2.selectbox("Filter by", ['No filter'] + list(cds.sampleData.columns))
                if filter_condition == 'No filter':
                    filter_categories = []
                else:
                    filter_categories = c2.multiselect(f'Which category(ies) of {filter_condition} to keep?',
                                                       list(cds.sampleData[filter_condition].unique()))
                if 'All' in condition_categories:
                    condition_categories = list(cds.sampleData[compare_condition].unique())
                genes = st.multiselect("Choose gene(s) of interest", cds.countData[gene_name].unique())
                if len(genes) * len(condition_categories) > 40:
                    st.write('Too many genes/categories to display, consider choosing fewer genes')
                else:
                    gene_df = cds.countData[cds.countData[gene_name].isin(genes)]
                    ab_sample_df = cds.sampleData[cds.sampleData[compare_condition].isin(condition_categories)]
                    if filter_categories:
                        ab_sample_df = ab_sample_df[ab_sample_df[filter_condition].isin(filter_categories)]
                    gene_df = (gene_df.melt(id_vars=[barcode, gene_name], value_name='log2CPM', var_name='sampleID')
                               .merge(ab_sample_df, how='inner', on='sampleID')
                               .sort_values(compare_condition))

                    col1, col2 = st.columns(2)
                    groupBy = col1.radio('Group by', [gene_name, compare_condition])
                    colorBy = [c for c in [gene_name, compare_condition] if c != groupBy][0]

                    # Choose Plot Type
                    box = "Box"
                    violin = "Violin"
                    plotType = col2.radio('Plot Type', (box, violin))
                    if plotType == box:
                        fig = barcode_abundance_box(gene_df, groupBy, colorBy, all_clrs)
                    if plotType == violin:
                        fig = barcode_abundance_violin(gene_df, groupBy, colorBy, all_clrs)

                    st.plotly_chart(fig, use_container_width=True)

