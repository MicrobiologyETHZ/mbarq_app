import streamlit as st
import pandas as pd
from sklearn.decomposition import PCA
import numpy as np
import plotly.express as px
from itertools import cycle
from pathlib import Path


@st.cache
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv(index=False).encode('utf-8')


def find_PCs(countData, sampleData, numPCs=2, numGenes=None, choose_by='variance'):
    """
    :param countData: each column is a sampleID, index is featureID
    :param sampleData:
    :param numPCs:
    :param numGenes:
    :return:
    """
    if numGenes:
        # calculate var for each, pick numGenes top var across samples -> df
        if choose_by == 'variance':
            genes = countData.var(axis=1).sort_values(ascending=False).head(numGenes).index
            df = countData.loc[genes].T
        else:
            pass
            # todo implement log2fc selection
    else:
        df = countData.T

    pca = PCA(n_components=numPCs)
    principalComponents = pca.fit_transform(df)
    pcs = [f'PC{i}' for i in range(1, numPCs + 1)]
    pDf = (pd.DataFrame(data=principalComponents, columns=pcs)
           .set_index(df.index))
    pc_var = {pcs[i]: round(pca.explained_variance_ratio_[i] * 100, 2) for i in range(0, numPCs)}
    pDf2 = pDf.merge(sampleData, how="left", left_index=True, right_index=True)
    return pDf2, pc_var


def app():
    # CSS to inject contained in a string
    hide_dataframe_row_index = """
                <style>
                .row_heading.level0 {display:none}
                .blank {display:none}
                </style>
                """

    # Inject CSS with Markdown
    st.markdown(hide_dataframe_row_index, unsafe_allow_html=True)

    st.markdown(""" # Exploratory Analysis
        
    ### Visualizing barcode count data.
        
    - Takes in a **CSV** file of merged counts produced by `mbarq count` + `mbarq merge`. The first column must contain the barcodes, the second column must contain barcode annotation (e.g. gene name or locus tag). All other columns must be sample names. Example structure:
    """)

    test = pd.DataFrame([['ACACACGT', 'abcD', '450', '700']], columns=['barcode', 'geneName', 'sample1', 'sample2'])
    st.dataframe(test)

    st.markdown("""
    - Takes in a **CSV** file containing sample data. First column must contain sample names that correspond to sample names in the count file. Example structure:
    """)
    test = pd.DataFrame([['sample1', 'batch1', 'treatment'], ['sample2', 'batch2', 'control']],
                        columns=['sampleID', 'batch', 'treatment'])
    st.dataframe(test)
    st.markdown("""
    
    - Merged count table produced by `mbarq` will contain barcodes found in the mapping file, as well as unannotated barcodes (e.g. control spike-ins, artifacts). Only annotated barcodes are used for the exploratory analysis.
    
    - Simple TSS normalisation and log2 transformation is performed
    - For PCA plot, you can how many barcodes are used for the analysis, as well as which components to visualise. Scree plot shows the % of variance explained by each of the PCs. 
    - For Barcode Abundance, normalised barcode counts can be visualised for any gene of interest and compared across different sample data variables. 

        """)

    with st.container():
        st.subheader('Load your own data or browse the example data set')
        data_type = st.radio('Choose dataset to show', ['Look at an example', 'Load my data'], index=1)
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
        if not cfile or not mfile:
            st.stop()

        # Requirements: first column has barcodes, second column has attributes in the countData, rest need to be sampleIDs.
        # first column are sampleIDs column in the sampleData
        # Will only look at barcodes that were mapped to a feature / alternatively filter out low counts?
        df = pd.read_csv(cfile)
        countData = (df.rename({df.columns[0]: 'barcode'}, axis=1)
                     .dropna()
                     .drop([df.columns[1]], axis=1)
                     .drop_duplicates()
                     .set_index('barcode'))
        to_log = st.checkbox('Apply log2 transform (recommended)', value=True)
        if to_log:
            countData = np.log2(countData / countData.sum() * 1000000 + 0.5)
        else:
            countData = countData / countData.sum() * 1000000
        sampleData = pd.read_csv(mfile, index_col=0)
        sampleData.index.name = 'sampleID'

    st.write('## PCA plot')
    with st.expander('Show PCA'):
        c1, c2 = st.columns((4, 1))
        c2.write('### PCA Options')
        numPCs = c2.slider("Select number of Principal Components", min_value=2, max_value=50, value=10)
        numGenes = c2.slider("Number of genes to use", value=250, max_value=countData.shape[0])
        # choose_by = c2.selectbox('Choose genes based on highest', ['variance', 'log2FoldChange (not implemented)'])
        pDf, pc_var = find_PCs(countData, sampleData, numPCs, numGenes, 'variance')
        pcX_labels = [f'PC{i}' for i in range(1, numPCs + 1)]
        expVars = [c for c in pDf.columns if c not in pcX_labels]
        pcX = c2.selectbox('X-axis component', pcX_labels)
        pcY = c2.selectbox('Y-axis component', [pc for pc in pcX_labels if pc != pcX])
        pcVar = c2.radio('Variable to highlight', expVars)
        fig = px.scatter(pDf, x=pcX, y=pcY, color=pcVar,
                         labels={pcX: f'{pcX}, {pc_var[pcX]} % Variance',
                                 pcY: f'{pcY}, {pc_var[pcY]} % Variance'},
                         height=700, hover_data=expVars, hover_name=pDf.index)
        fig.update_layout(autosize=True, font=dict(size=18), paper_bgcolor='rgba(0,0,0,0)',
                          )
        fig.update_traces(marker=dict(size=12,
                                      line=dict(width=2,
                                                color='DarkSlateGrey')),
                          selector=dict(mode='markers'))
        c1.write(f'### {pcX} vs {pcY}, highlighting {pcVar}')
        c1.plotly_chart(fig, use_container_width=True)

        c3, c4 = st.columns(2)
        pDf_sum = pDf.groupby(pcVar).median()
        varDf = pd.DataFrame.from_dict(pc_var, orient='index').reset_index()
        varDf.columns = ['PC', '% Variance']
        fig2 = px.line(varDf, x='PC', y='% Variance', markers=True,
                       labels={'PC': ''})
        fig2.update_traces(marker=dict(size=12,
                                       line=dict(width=2,
                                                 color='DarkSlateGrey')))
        c3.write('### Scree Plot')
        c3.plotly_chart(fig2)
        c4.write(f'### PCs summarized by {pcVar}')
        c4.plotly_chart(px.imshow(pDf_sum), use_container_width=True)

    st.write('## Barcode Abundance')
    with st.expander('Show Barcode Abundance'):
        df = df.dropna()
        barcode = df.columns[0]
        gene_name = df.columns[1]
        df = df.set_index([barcode, gene_name])
        df = np.log2(df / df.sum() * 1000000 + 0.5)
        sampleDataAb = sampleData.reset_index()
        df = df.reset_index()
        c1, c2 = st.columns(2)
        compare_by = c1.selectbox('Compare by', sampleDataAb.columns)
        color_by = c2.selectbox('Color by', [barcode] + list(sampleDataAb.columns))
        genes = st.multiselect("Choose gene(s) of interest", df[gene_name].unique())

        if not genes:
            st.stop()
        c3, c4 = st.columns(2)
        for col, gene in zip(cycle([c3, c4]), genes):
            gene_df = df[df[gene_name] == gene]
            gene_df = (gene_df.melt(id_vars=[barcode, gene_name], value_name='log2CPM', var_name='sampleID')
                       .merge(sampleData, how='left', on='sampleID'))
            gene_df = gene_df.sort_values(compare_by)
            fig = px.strip(gene_df, title=gene, x=compare_by, y='log2CPM', color=color_by,
                           hover_data=[barcode] + list(sampleData.columns), stripmode='overlay')
            fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                              font=dict(size=16))
            fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='LightGrey')
            col.plotly_chart(fig, use_container_width=True)
