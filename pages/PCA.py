import streamlit as st
import pandas as pd
from sklearn.decomposition import PCA
import numpy as np
import plotly.express as px

from pathlib import Path

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
    pcs = [f'PC{i}' for i in range(1, numPCs+1)]
    pDf = (pd.DataFrame(data=principalComponents, columns=pcs)
           .set_index(df.index))
    pc_var = {pcs[i]: round(pca.explained_variance_ratio_[i] * 100, 2) for i in range(0, numPCs)}
    pDf2 = pDf.merge(sampleData, left_index=True, right_index=True)
    return pDf2, pc_var


def app():
    st.write('# PCA')
    cfile = st.file_uploader('Load merged count file')
    mfile = st.file_uploader('Load metadata')
    cfile = "/Users/ansintsova/git_repos/mbarq_app/data/SL1344_test/library_10_1_count_batchcorrected.txt"
    mfile = "/Users/ansintsova/git_repos/mbarq_app/data/SL1344_test/library_10_1_sampleData.csv"

    # Requirements: barcode and gene columns in the countData, rest need to be sampleIDs.
    # sampleID column in the sampleData

    df = (pd.read_table(cfile).drop(['gene'], axis=1)
          .drop_duplicates()
          .set_index('barcode'))
    countData = np.log2(df/df.sum()*1000000 +0.5)

    sampleData = pd.read_csv(mfile).set_index('sampleID')
    c1, c2 = st.columns((4, 1))
    c2.write('### PCA Options')
    numPCs = c2.slider("Select number of Principal Components", min_value=2, max_value=50, value=10)
    numGenes = c2.slider("Number of genes to use", value=500, max_value=countData.shape[0])
    #choose_by = c2.selectbox('Choose genes based on highest', ['variance', 'log2FoldChange (not implemented)'])
    pDf, pc_var = find_PCs(countData, sampleData, numPCs, numGenes, 'variance')
    pcX_labels = [f'PC{i}' for i in range(1, numPCs+1)]
    expVars = [c for c in pDf.columns if c not in pcX_labels]
    pcX = c2.selectbox('X-axis component', pcX_labels)
    pcY = c2.selectbox('Y-axis component', [pc for pc in pcX_labels if pc != pcX])
    pcVar = c2.radio('Variable to highlight', expVars)
    fig = px.scatter(pDf, x=pcX, y=pcY, color=pcVar,
                     labels ={pcX: f'{pcX}, {pc_var[pcX]} % Variance',
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


