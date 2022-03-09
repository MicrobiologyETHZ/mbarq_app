import streamlit as st
from pathlib import Path
import pandas as pd
import numpy as np
import plotly.express as px
from itertools import cycle
import dash_bio as dashbio


def load_data(design_file, count_file):
    sampleData = pd.read_csv(design_file, index_col=0)
    countData = pd.read_csv((count_file), index_col=0)
    return sampleData, countData


def app():
    st.write('# Differentical Abundance Results')
    rfile = st.file_uploader('Load the analysis file')
    rfile = "/Users/ansintsova/git_repos/mbarq_app/data/SL1344_test/library_10_1-unfiltered-results.kegg.csv"

    clrs = px.colors.qualitative.D3
    fdf = pd.read_csv(rfile, index_col=0)

    # there are fdr, lfc, contrast and KEGG_Pathway columns
    fdf['log10FDR'] = -10*np.log10(fdf['fdr'])
    fdf["KEGG_Pathway"] = fdf.KEGG_Pathway.str.split("- Salmonella", expand=True)[0]
    fdf['Library'] = 'library_10_1'
    contrasts = fdf['contrast'].sort_values().unique()
    contrast_col, lfc_col, fdr_col = st.columns(3)
    contrast_to_show = contrast_col.selectbox('Select a contrast', contrasts)
    fdr = fdr_col.number_input('FDR cutoff', value=0.05)
    lfc_th = lfc_col.number_input('Log FC cutoff (absolute)', value=1)

    df = fdf[fdf.contrast == contrast_to_show].copy()
    df['hit'] = ((abs(df['lfc']) > lfc_th) & (df['fdr'] < fdr))

    show_kegg =st.selectbox('Show KEGG Pathway', ['All'] + list(df.KEGG_Pathway.unique()))

    if show_kegg != 'All':
        df = df[df.KEGG_Pathway == show_kegg]
    df = df.sort_values('lfc').reset_index().reset_index().rename({'index':'ranking'}, axis=1)

    fig = px.scatter(df, x='ranking', y='lfc', color='hit', size='log10FDR',
                     height=700,
                     color_discrete_map={
                         True: clrs[1],
                         False: clrs[0]},
                     hover_name = 'Name',
                     title = f"{contrast_to_show} - {show_kegg}",
                     hover_data={'lfc':True,
                                 'log10FDR': False,
                                'ranking': False,
                                 'fdr': True},
                     labels = {"ranking": '', 'lfc': 'Log2 FC'}
                     )
    fig.add_hline(y=0, line_width=2, line_dash="dash", line_color="grey")
    fig.update_xaxes(showticklabels=False)
    fig.update_layout( {'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                      font=dict(size=18))
    fig.update_traces(marker=dict(
                                  line=dict(width=1,
                                            color='DarkSlateGrey')),
                      selector=dict(mode='markers'))
    st.plotly_chart(fig, use_container_width=True)
    #

    c1, c2 = st.columns(2)
    show_kegg_heat1 = c1.selectbox('Show KEGG Pathway',  list(fdf.KEGG_Pathway.dropna().unique()), key='leftKegg')
    heatDfLeft = fdf[fdf.KEGG_Pathway == show_kegg_heat1][['Name', 'lfc', 'contrast']].pivot(index="Name",
                                                                                             columns='contrast',
                                                                                             values='lfc')
    show_kegg_heat2 = c2.selectbox('Show KEGG Pathway', list(fdf.KEGG_Pathway.dropna().unique()), key='rightKegg')


    heatDfRight = fdf[fdf.KEGG_Pathway == show_kegg_heat2][['Name', 'lfc', 'contrast']].pivot(index="Name",
                                                                                             columns='contrast',
                                                                                            values='lfc')

    fig2 = px.imshow(heatDfLeft, color_continuous_scale=px.colors.sequential.RdBu_r, color_continuous_midpoint=0,
                         width=600, height=900)
    fig3 = px.imshow(heatDfRight, color_continuous_scale=px.colors.sequential.RdBu_r, color_continuous_midpoint=0,
                         width=600, height=900)
    c1.plotly_chart(fig2, use_container_width=False)
    c2.plotly_chart(fig3, use_container_width=False)

    # compContrasts = st.multiselect('Select contrasts to compare', contrasts)
    # if not compContrasts:
    #     st.stop()
    # c1, c2 = st.columns(2)
    # filters = {}
    # for col, contrast in zip(cycle([c1, c2]), compContrasts):
    #    col.write(contrast)
    #    l = col.number_input('LFC cutoff', value=2, key=f'{contrast}_lfc')
    #    f = col.number_input('FDR cutoff', value=0.05, key=f'{contrast}_fdr')
    #    filters[contrast] = (l, f)
    #
    # compDfs = []
    # for key, value in filters.items():
    #     if value[0] > 0:
    #         df = fdf[(fdf.contrast == key) & (fdf.lfc > value[0])&(fdf.fdr < value[1])]
    #     else:
    #         df = fdf[(fdf.contrast == key) & (fdf.lfc < value[0]) & (fdf.fdr < value[1])]
    #     compDfs.append(df)
    #
    # vennDf = pd.concat(compDfs)
    # vennDf = vennDf.reset_index().pivot(index='index', columns='contrast', values='lfc').dropna()
    #
    # if vennDf.empty:
    #     st.write('No genes matching the filtering criteria found')
    #     st.stop()
    # st.write(vennDf.shape)
    # columns = list(vennDf.columns)
    # rows = list(vennDf.index)
    # fix error if only one gene or 1 sample

    #
    # cluster = 'all' if len(columns) > 1 else 'row'
    # clustergram = dashbio.Clustergram(
    #     data=vennDf.loc[rows].values,
    #     row_labels=rows,
    #     column_labels=columns,
    #     color_threshold={
    #         'row': 250,
    #         'col': 700
    #     },
    #     height=800,
    #     width=700,
    #     center_values=False,
    #     cluster=cluster
    # )
    #
    # if vennDf.shape[0] < 100:
    #     st.plotly_chart(clustergram, use_container_width=True)
    # else:
    #     st.write('Too many genes to display. Download table?')
