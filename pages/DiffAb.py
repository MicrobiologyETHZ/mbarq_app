import streamlit as st
from pathlib import Path
import pandas as pd
import numpy as np
import plotly.express as px
from itertools import cycle
import dash_bio as dashbio


def process_library_results(result_file, file_type='mageck'):
    df = pd.read_csv(result_file, index_col=0)
    if 'library' not in df.columns:
        library_name = st.text_input("Add library name?", value='')
        df['library'] = library_name
    fixed_col_names = {'mageck': ['library', 'lfc', 'fdr', 'contrast', 'KEGG_Pathway']}
    missing_cols = [c for c in fixed_col_names[file_type] if c not in df.columns]
    if len(missing_cols) > 0:
        st.markdown(
            f"""The following columns are missing from the map files: {', '.join(missing_cols)}. 
            Please rename the columns/rerun mBARq and try again. Skipping {result_file.name}""")
        return pd.DataFrame()
    df = df.rename({df.columns[0]: 'Name'}, axis=1)
    return df

def app():
    st.markdown(""" # Library Map

    ### Write a little explanation of what this page shows.

    - Required Inputs
    - How to use the graph

        """)

    with st.container():
        st.subheader('Browse an example data set')
        data_type = st.radio('Choose dataset to show', ['Look at an example'], index=0)
        if data_type == 'Load my data':
            result_files = st.file_uploader('Upload library map file', accept_multiple_files=True)
        else:
            result_files = [Path("/Users/ansintsova/git_repos/mbarq_app/data/SL1344_test/library_10_1-unfiltered-results.kegg.csv")]
            st.subheader('Example mapping file')
            st.write(pd.read_csv(result_files[0], index_col=0).sample(5))
        if len(result_files) < 1:
            st.stop()

    with st.container():
        clrs = px.colors.qualitative.D3
        result_dfs = []
        for uploaded_result in result_files:
            st.write(f"Processing {uploaded_result.name}")
            df = process_library_results(uploaded_result)
            result_dfs.append(df)
        fdf = pd.concat(result_dfs)
        fdf['log10FDR'] = -10*np.log10(fdf['fdr'])
        fdf["KEGG_Pathway"] = fdf.KEGG_Pathway.str.split("- Salmonella", expand=True)[0] #todo needs to be done ahead of the app
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
                         hover_name='Name',
                         title=f"{contrast_to_show} - {show_kegg}",
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
