import streamlit as st
import pandas as pd
from sklearn.decomposition import PCA
import numpy as np
import plotly.express as px
from itertools import cycle
from pathlib import Path
import requests
from time import sleep


def define_color_scheme():
    alphabetClrs = px.colors.qualitative.Alphabet
    clrs = ["#f7ba65", "#bf4713", "#9c002f", "#d73d00", "#008080", "#004c4c"]
    colors = {'grey': alphabetClrs[8],
              'light_yellow': clrs[0],
              'darko': clrs[1],
              'maroon': clrs[2],
              'brighto': clrs[3],
              'teal': clrs[4],
              'darkteal': clrs[5]
              }
    all_clrs = [colors['brighto'], colors['teal'], colors['maroon']] + alphabetClrs[13:]
    return colors, alphabetClrs, all_clrs


def find_PCs(countData, sampleData, numPCs=2, numGenes=None, choose_by='variance'):
    """
    :param countData: each column is a sampleID, index is featureID
    :param sampleData:
    :param numPCs:
    :param numGenes:
    :return:
    """
    if countData.columns[0] != 'barcode':
        st.write('Barcode column not found.')
        return ()
    if sampleData.columns[0] != 'sampleID':
        st.write('sampleID column not found')
        return ()
    df = countData.set_index('barcode')
    sampleData = sampleData.set_index('sampleID').apply(lambda x: x.astype('category'))
    if numGenes:
        # calculate var for each, pick numGenes top var across samples -> df
        if choose_by == 'variance':
            genes = df.var(axis=1).sort_values(ascending=False).head(int(numGenes)).index
            df = df.loc[genes].T
        else:
            pass
            # todo implement log2fc selection
    else:
        df = countData.T
    pca = PCA(n_components=numPCs)
    principalComponents = pca.fit_transform(df)
    pcs = [f'PC{i}' for i in range(1, numPCs + 1)]
    pcDf = (pd.DataFrame(data=principalComponents, columns=pcs)
              .set_index(df.index))
    pcVar = {pcs[i]: round(pca.explained_variance_ratio_[i] * 100, 2) for i in range(0, numPCs)}
    pcDf = pcDf.merge(sampleData, how="left", left_index=True, right_index=True)
    return pcDf, pcVar


def graph_library_map(map_df, attr_names, chr_col='chr', color_by_cols=('in CDS', 'library')):
    attr_names = set(attr_names)
    colors, alphabetClrs, all_clrs = define_color_scheme()
    graph_type = st.radio("Choose graph", ['Coverage Histogram', 'Individual Insertions'])
    c1, c2 = st.columns(2)
    seqid = c1.selectbox('Choose sequence to display', map_df[chr_col].unique())
    df_to_show = map_df[map_df[chr_col] == seqid].sort_values("insertion_site")
    if graph_type == 'Individual Insertions':
        colorby = c2.selectbox('Color by', color_by_cols)
        fig = px.scatter(df_to_show, x='insertion_site', y='abundance_in_mapping_library', color=colorby, log_y=True,
                         height=600,
                         color_discrete_sequence=all_clrs,
                         hover_data=attr_names, template='plotly_white',
                         labels={'insertion_site': 'Position, bp', 'abundance_in_mapping_library': 'Read Counts'})
        fig.update_traces(marker=dict(size=8, line=dict(width=2, color='DarkSlateGrey')),
                          selector=dict(mode='markers'), opacity=0.9)

    elif graph_type == 'Coverage Histogram':
        num_bins = c2.number_input('Number of bins', value=100)
        st.write(int(num_bins))
        fig = px.histogram(df_to_show, x='insertion_site', nbins=int(num_bins),
                           labels={'insertion_site': 'Position, bp'}, color_discrete_sequence=[colors['teal']])
        fig.update_layout(bargap=0.1)

    fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'},
                      autosize=True,
                      font=dict(size=24))
    st.plotly_chart(fig, use_container_width=True)


def barcode_abundance(geneDf, groupBy, colorBy, colorSeq):
        fig = px.box(geneDf, x=groupBy, y='log2CPM', color=colorBy,
                     hover_data=geneDf.columns, points='all',
                     color_discrete_sequence=colorSeq,)
        fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                          font=dict(size=16))
        fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='LightGrey')
        return fig


def pca_figure(pcDf, pcX, pcY, pcVarHi, pcVar, pcSym, expVars, colorSeq):
    fig = px.scatter(pcDf, x=pcX, y=pcY, color=pcVarHi, symbol=pcSym,
                     labels={pcX: f'{pcX}, {pcVar[pcX]} % Variance',
                             pcY: f'{pcY}, {pcVar[pcY]} % Variance'},
                     color_discrete_sequence=colorSeq,
                     template='plotly_white',
                     height=800, hover_data=expVars, hover_name=pcDf.index)
    fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'},
                      autosize=True,
                      font=dict(size=20))
    fig.update_traces(marker=dict(size=24,
                                  line=dict(width=2,
                                            color='DarkSlateGrey'), opacity=0.9),
                      selector=dict(mode='markers'))
    varDf = pd.DataFrame.from_dict(pcVar, orient='index').reset_index()
    varDf.columns = ['PC', '% Variance']
    fig2 = px.line(varDf, x='PC', y='% Variance', markers=True,
                   labels={'PC': ''})
    fig2.update_traces(marker=dict(size=12,
                                   line=dict(width=2,
                                             color='DarkSlateGrey')))
    pcSum = pcDf.groupby(pcVarHi).median()
    fig3 = px.imshow(pcSum)
    return fig, fig2, fig3


def show_lfc_ranking(fdf, contrasts, libraries):
    colors, alphabetClrs, all_clrs = define_color_scheme()
    contrast_col, lfc_col, fdr_col, lfc_lib_col = st.columns(4)
    contrast_to_show = contrast_col.selectbox('Select a contrast', contrasts)
    library_to_show = lfc_lib_col.selectbox('Select library to show', libraries)
    fdr = fdr_col.number_input('FDR cutoff', value=0.05)
    lfc_th = lfc_col.number_input('Log FC cutoff (absolute)', min_value=0.0, step=0.5, value=1.0)
    df = fdf[(fdf.contrast == contrast_to_show) & (fdf.library == library_to_show)].copy()
    df['hit'] = ((abs(df['LFC']) > lfc_th) & (df['fdr'] < fdr))
    show_kegg = st.selectbox('Show KEGG Pathway', ['All'] + list(df.KEGG_Pathway.unique()))
    if show_kegg != 'All':
        df = df[df.KEGG_Pathway == show_kegg]
    df = df.sort_values('LFC').reset_index().reset_index().rename({'level_0': 'ranking'}, axis=1)
    fig = px.scatter(df, x='ranking', y='LFC', color='hit',
                     height=800,
                     color_discrete_map={
                         True: colors['teal'],
                         False: colors['grey']},
                     hover_name='Name',
                     title=f"{contrast_to_show} - {show_kegg}",
                     hover_data={'LFC': True,
                                 'log10FDR': False,
                                 'ranking': False,
                                 'fdr': True,
                                 'KEGG_Pathway': True},
                     labels={"ranking": '', 'LFC': 'Log2 FC'}
                     )
    fig.add_hline(y=0, line_width=2, line_dash="dash", line_color="grey")
    fig.update_xaxes(showticklabels=False)
    fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                      font=dict(size=18))
    fig.update_traces(marker=dict(size=20,
                                  line=dict(width=0.2,
                                            color='DarkSlateGrey'), opacity=0.8),
                      selector=dict(mode='markers'))
    #st.plotly_chart(fig, use_container_width=True)
    return fig, df[df.hit == True]


def link_to_string(hits_df, st_col, lfc_col='LFC', gene_name='Name'):
    up = st_col.radio('Up or Down?', ('Upregulated Only', 'Downregulated Only', 'Both'))
    if up == 'Upregulated Only':
        hits_df = hits_df[hits_df[lfc_col] > 0]
    elif up == 'Downregulated Only':
        hits_df = hits_df[hits_df[lfc_col] < 0]

    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = 'tsv-no-header'
    method = 'get_link'
    if gene_name:
        my_genes = set(hits_df[gene_name].values)
    else:
        my_genes = list(hits_df.index)
    request_url = "/".join([string_api_url, output_format, method])
    species = st_col.number_input("NCBI species taxid", value=99287, help='Salmonella Typhimurium: 99287')
    params = {
        "identifiers": "\r".join(my_genes),  # your protein
        "species": species,  # species NCBI identifier
        "network_flavor": "confidence",  # show confidence links
        "caller_identity": "explodata"  # your app name
    }
#
    if st_col.button('Get STRING network'):
        network = requests.post(request_url, data=params)
        network_url = network.text.strip()
        st_col.markdown(f"[Link to STRING network]({network_url})")
        sleep(1)


def plot_rank(rank_df, colors, hover_dict, gene_id):
    fig = px.scatter(rank_df, x='ranking', y='LFC_median', color='hit', symbol='contrast',
                     height=800,
                     color_discrete_map={
                         True: colors['teal'],
                         False: colors['grey']},
                     hover_name=gene_id,
                     hover_data=hover_dict,
                     labels={"ranking": '', 'LFC_median': 'Log2 FC'}
                     )
    fig.add_hline(y=0, line_width=2, line_dash="dash", line_color="grey")
    fig.update_xaxes(showticklabels=False)
    fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                      font=dict(size=18))
    fig.update_traces(marker=dict(size=14,
                                  line=dict(width=1,
                                            color='DarkSlateGrey'), opacity=0.8),
                      selector=dict(mode='markers'))
    return fig


def plot_position(position_df, hover_dict, gene_id):
    fig = px.scatter(position_df, x='Start', y='LFC_median', color='contrast', size='library_nunique',
                     height=800,
                     color_discrete_sequence=px.colors.qualitative.D3,
                     hover_name=gene_id,
                     hover_data=hover_dict,
                     labels={'LFC_median': 'Log2 FC'}
                     )
    fig.add_hline(y=0, line_width=2, line_dash="dash", line_color="grey")

    fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                      font=dict(size=18))
    fig.update_traces(marker=dict(
        line=dict(width=1,
                  color='DarkSlateGrey'), opacity=0.8),
        selector=dict(mode='markers'))
    return fig


def plot_heatmap(heatDf):
    fig = px.imshow(heatDf, color_continuous_scale=px.colors.diverging.Geyser,
                     color_continuous_midpoint=0,

                     width=1000, height=900)
    fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                       font=dict(size=10))
    return fig