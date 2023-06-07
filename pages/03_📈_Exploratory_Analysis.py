import streamlit as st
from scripts.datasets import CountDataSet, convert_df, define_color_scheme
import pandas as pd
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
        
        - A **csv** file of merged counts produced by `mbarq count` + `mbarq merge`. For instruction on how to generate this file, please see [here]({count_url}).
        - The first column must contain the barcodes, the second column must contain gene identifier (ex. locus tag). 
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
        - First column must contain sample names that correspond to sample names in the count file.  
        - All other columns will be read in as metadata.
        #
        
        Example structure:
        """)
        test = pd.DataFrame([['sample1', 'treatment'], ['sample2', 'control']],
                            columns=['sampleID', 'treatment'])
        st.table(test)
        st.markdown("""

        - Merged count table produced by `mbarq` will contain barcodes found in the mapping file, as well as unannotated barcodes (e.g. control spike-ins, artifacts). Only annotated barcodes are used for the exploratory analysis.
        - Simple TSS normalisation and log2 transformation is performed.
        - For PCA plot, you can choose how many barcodes are used for the analysis, as well as which components to visualise. Scree plot shows the % of variance explained by each of the PCs. 
        - For Barcode Abundance, normalised barcode counts can be visualised for any gene of interest and compared across different sample data variables. 
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
            with st.expander('Show PCA'):
                aC, sushi_clrs, all_clrs = define_color_scheme()
                st.write('### PCA Options')
                c1, c2, c3, c4 = st.columns(4)

                numPCs = c1.number_input("Number of PCs", min_value=2, max_value=50, value=10)
                numGenes = c2.number_input("Number of genes to use", min_value=int(numPCs),
                                           value=int(min(250, cds.count_data.shape[0])),
                                           max_value=int(cds.count_data.shape[0]))
                chooseBy = 'variance'
                numGenes = int(numGenes)
                numPCs = int(numPCs)
                pcDf, pcVar = cds.get_principal_components(numPCs, numGenes, chooseBy)

                missingMeta = " ,".join(list(pcDf[pcDf.isna().any(axis=1)].index))
                if missingMeta:
                    st.write(f"The following samples have missing_metadata and will not be shown: {missingMeta}")
                pcDf = pcDf[~pcDf.isna().any(axis=1)]  # todo this should be included in the function
                pcxLabels = [f'PC{i}' for i in range(1, numPCs + 1)]
                expVars = [c for c in pcDf.columns if c not in pcxLabels]
                w, h, font_size = None, None, 24
                if st.checkbox('Modify PCA plot dimensions'):
                    wc, hc, font = st.columns(3)
                    w = wc.number_input('Width', min_value=200, max_value=1000, value=800)
                    h = hc.number_input('Height', min_value=200, max_value=1000, value=400)
                    font_size = font.number_input("Font size", min_value=8, max_value=40, value=24)
                pcX = c1.selectbox('X-axis component', pcxLabels)
                pcY = c2.selectbox('Y-axis component', [pc for pc in pcxLabels if pc != pcX])
                pcVarHi = c3.radio('Variable to highlight', expVars)
                pcSym = c4.radio('Variable to show as symbol', [None] + expVars)
                pcDf = pcDf.sort_values(pcVarHi)
                fig1, fig2, fig3 = cds.pca_figure(pcDf, pcX, pcY, pcVarHi, pcVar, pcSym, expVars,
                                                  all_clrs, w, h, font_size=font_size)
                st.plotly_chart(fig1, use_container_width=False)
                c5, c6 = st.columns(2)
                c5.write('### Scree Plot')
                c5.plotly_chart(fig2)
                c6.write(f'### PCs summarized by {pcVarHi}')
                c6.plotly_chart(fig3, use_container_width=True)
            # BARCODE ABUNDANCE
            st.write('## Barcode Abundance')
            with st.expander('Show Barcode Abundance'):
                # Process the dataframe
                # abDf = count_data.dropna()
                # Get user input
                c1, c2 = st.columns(2)
                compare_condition = c1.selectbox('Which conditions to compare?', cds.sample_data.columns)
                condition_categories = c1.multiselect(f'Categories of {compare_condition} to display',
                                                      ['All'] + list(cds.sample_data[compare_condition].unique()),
                                                      default='All')
                filter_condition = c2.selectbox("Filter by", ['No filter'] + list(cds.sample_data.columns),
                                                index=0)
                if filter_condition == 'No filter':
                    filter_categories = []
                else:
                    filter_categories = c2.multiselect(f'Which category(ies) of {filter_condition} to keep?',
                                                       list(cds.sample_data[filter_condition].unique()))
                if 'All' in condition_categories:
                    condition_categories = list(cds.sample_data[compare_condition].unique())
                genes = st.multiselect("Choose gene(s) of interest", cds.count_data[cds.gene_name_col].unique())
                if len(genes) * len(condition_categories) > 40:
                    st.write('Too many genes/categories to display, consider choosing fewer genes.')
                else:
                    gene_df = cds.count_data[cds.count_data[cds.gene_name_col].isin(genes)]
                    ab_sample_df = cds.sample_data[cds.sample_data[compare_condition].isin(condition_categories)]
                    if filter_categories:
                        ab_sample_df = ab_sample_df[ab_sample_df[filter_condition].isin(filter_categories)]
                    gene_df = (gene_df.melt(id_vars=[cds.barcode_col, cds.gene_name_col],
                                            value_name='log2CPM', var_name=cds.sample_id_col)
                               .merge(ab_sample_df, how='inner', on=cds.sample_id_col)
                               .sort_values(compare_condition))

                    col1, col2 = st.columns(2)
                    groupBy = col1.radio('Group by', [cds.gene_name_col, compare_condition])
                    colorBy = [c for c in [cds.gene_name_col, compare_condition] if c != groupBy][0]

                    # Choose Plot Type
                    box = "Box"
                    violin = "Violin"
                    plotType = col2.radio('Plot Type', (box, violin))
                    if plotType == box:
                        fig = cds.barcode_abundance_plot(gene_df, groupBy, colorBy, all_clrs)
                    if plotType == violin:
                        fig = cds.barcode_abundance_plot(gene_df, groupBy, colorBy, all_clrs, box=False)
                    st.plotly_chart(fig, use_container_width=True)
app()