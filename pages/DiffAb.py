import streamlit as st
from pathlib import Path
import pandas as pd
import numpy as np
import plotly.express as px
from itertools import cycle


def process_library_results(result_file, file_type='mageck'):
    df = pd.read_csv(result_file, index_col=0)
    if 'library' not in df.columns:
        library_name = st.text_input("Add library name?", value='')
        df['library'] = library_name
    fixed_col_names = {'mageck': ['library', 'LFC', 'neg_selection_fdr', 'pos_selection_fdr', 'contrast']}
    missing_cols = [c for c in fixed_col_names[file_type] if c not in df.columns]
    if len(missing_cols) > 0:
        st.markdown(
            f"""The following columns are missing from the results files: {', '.join(missing_cols)}. 
            Please rename the columns/rerun mBARq and try again. Skipping {result_file.name}""")
        return pd.DataFrame()
    df = df.rename({df.columns[0]: 'Name'}, axis=1)
    return df


def process_gmt(gmt_file, identifier='Name', path_name=1):
    genes = []
    paths = []
    with open(gmt_file, 'r', errors='ignore') as fh:
        for line in fh.readlines():
            for gene in line.split('\t')[path_name:]:
                genes.append(gene.strip("\n"))
                paths.append(":".join(line.split('\t')[0:path_name]))
    return pd.DataFrame([genes, paths], index=[identifier, 'KEGG_Pathway']).T


def app():
    url = "https://doi.org/10.1186/s13059-014-0554-4"
    st.markdown(f""" # Differential Gene Fitness

    ### Visualize gene log fold changes (LFCs) between two conditions of interest

    - Statistical analysis is performed using [MAGeCK]({url}). MAGeCK uses a modified version of Robust Ranking Aggregation algorithm to identify negativly and positively selected genes.
    
    - :hourglass: To add KEGG annotation to the results, `emapper` annotation file must be provided 

    - The central graph displays LFC for each gene ranked from lowest to highest, the size of the circle is inversely proportional to the FDR
    
    - The heatmaps show LFC for all the genes for a specific KEGG pathway for all the contrasts (comparisions)
    
        """)

    with st.container():
        st.subheader('Browse an example data set')
        data_type = st.radio('Choose dataset to show', ['Look at an example'], index=0)
        if data_type == 'Load my data':
            result_files = st.file_uploader('Upload results file', accept_multiple_files=True)
            gmt_file = st.file_uploader('Upload gmt file')
        else:
            result_files = [Path("examples/example_rra_results.csv")]
            gmt_file = Path("examples/04-03-2022-SL1344-KEGG-API.gmt")
            identifier = 'Name'
            st.subheader('Example results file')
            ex_df = (pd.read_csv(result_files[0], index_col=0))
            ex_df = ex_df[ex_df.contrast == 'd1'].sort_values('number_of_barcodes', ascending=False)
            st.write(ex_df.head(5))
        if len(result_files) < 1:
            st.stop()
        clrs = px.colors.qualitative.D3
        result_dfs = []
        for uploaded_result in result_files:
            st.write(f"Processing {uploaded_result.name}")
            df = process_library_results(uploaded_result)
            result_dfs.append(df)
        fdf = pd.concat(result_dfs)
        if gmt_file:
            gmt_df = process_gmt(gmt_file, identifier, path_name=2)
            fdf = fdf.merge(gmt_df, how='left', on=identifier)
        else:
            if "KEGG_Pathway" not in list(fdf.columns):
                fdf['KEGG_Pathway'] = np.nan


    st.subheader("Results by log fold change")
    with st.expander('LFC rankings'):
        fdf['fdr'] = np.where(fdf['LFC'] < 0, fdf['neg_selection_fdr'], fdf['pos_selection_fdr'])
        fdf['log10FDR'] = -10*np.log10(fdf['fdr'])
        contrasts = fdf['contrast'].sort_values().unique()
        contrast_col, lfc_col, fdr_col = st.columns(3)
        contrast_to_show = contrast_col.selectbox('Select a contrast', contrasts)
        fdr = fdr_col.number_input('FDR cutoff', value=0.05)
        lfc_th = lfc_col.number_input('Log FC cutoff (absolute)', min_value=0.0, step=0.5, value=1.0)
        df = fdf[fdf.contrast == contrast_to_show].copy()
        df['hit'] = ((abs(df['LFC']) > lfc_th) & (df['fdr'] < fdr))
        show_kegg = st.selectbox('Show KEGG Pathway', ['All'] + list(df.KEGG_Pathway.unique()))
        if show_kegg != 'All':
            df = df[df.KEGG_Pathway == show_kegg]
        df = df.sort_values('LFC').reset_index().reset_index().rename({'level_0': 'ranking'}, axis=1)
        fig = px.scatter(df, x='ranking', y='LFC', color='hit', size='log10FDR',
                         height=700,
                         color_discrete_map={
                             True: clrs[1],
                             False: clrs[0]},
                         hover_name='Name',
                         title=f"{contrast_to_show} - {show_kegg}",
                         hover_data={'LFC': True,
                                     'log10FDR': False,
                                    'ranking': False,
                                     'fdr': True},
                         labels={"ranking": '', 'LFC': 'Log2 FC'}
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

    st.subheader("Pathway Heatmaps")
    with st.expander('Show LFC Heatmaps'):
        c1, c2 = st.columns(2)
        show_kegg_heat1 = c1.selectbox('Show KEGG Pathway',  list(fdf.KEGG_Pathway.dropna().unique()), key='leftKegg')
        show_kegg_heat2 = c2.selectbox('Show KEGG Pathway', list(fdf.KEGG_Pathway.dropna().unique()), key='rightKegg')
        if not show_kegg_heat1:
            st.markdown('No KEGG Annotation found :thinking_face:')
        else:
            heatDfLeft = fdf[fdf.KEGG_Pathway == show_kegg_heat1][['Name', 'LFC', 'contrast']].pivot(index="Name",
                                                                                                     columns='contrast',
                                                                                                     values='LFC')
            heatDfRight = fdf[fdf.KEGG_Pathway == show_kegg_heat2][['Name', 'LFC', 'contrast']].pivot(index="Name",
                                                                                                     columns='contrast',
                                                                                                    values='LFC')

            fig2 = px.imshow(heatDfLeft, color_continuous_scale=px.colors.sequential.RdBu_r, color_continuous_midpoint=0,
                                 width=600, height=900)
            fig3 = px.imshow(heatDfRight, color_continuous_scale=px.colors.sequential.RdBu_r, color_continuous_midpoint=0,
                                 width=600, height=900)
            c1.plotly_chart(fig2, use_container_width=False)
            c2.plotly_chart(fig3, use_container_width=False)

    st.subheader("Browse results as a table")
    with st.expander('Tabular Results'):
        compContrasts = st.multiselect('Select contrasts to display', contrasts)
        if not compContrasts:
            st.write("No contrast selected")
        else:
            c1, c2 = st.columns(2)
            filters = {}
            for col, contrast in zip(cycle([c1, c2]), compContrasts):
               col.write(contrast)
               l = col.number_input('LFC cutoff', value=-1.0, step=0.5, key=f'{contrast}_lfc')
               f = col.number_input('FDR cutoff', value=0.05, step=0.01,  key=f'{contrast}_fdr')
               filters[contrast] = (l, f)

            comp_genes = []
            for key, value in filters.items():
                if value[0] > 0:
                    genes = set(fdf[(fdf.contrast == key) & (fdf['LFC'] > value[0]) & (fdf.fdr < value[1])][identifier].values)
                else:
                    genes = set(fdf[(fdf.contrast == key) & (fdf['LFC'] < value[0]) & (fdf.fdr < value[1])][identifier].values)
                comp_genes.append(genes)
            intersect_genes = set.intersection(*comp_genes)
            vennDf = fdf[fdf[identifier].isin(intersect_genes)].copy()
            vennDf = vennDf[[identifier, 'library', 'LFC', 'fdr', 'contrast']].drop_duplicates()
            vennDf2 = vennDf.pivot(index='Name', columns='contrast', values=['LFC', 'fdr'])
            st.write(vennDf2)


    st.subheader("Gene Selector")
    with st.expander('Gene Selector'):
        genes = st.multiselect("Choose gene(s) of interest", fdf[identifier].unique())
        if not genes:
            st.stop()
        c3, c4 = st.columns(2)

        for col, gene in zip(cycle([c3, c4]), genes):
            gene_df = fdf[fdf[identifier] == gene].copy()
            gene_df = gene_df[[identifier, 'library', 'contrast', 'LFC', 'fdr']].drop_duplicates()
            gene_df = gene_df.sort_values('contrast')
            fig = px.strip(gene_df, title=gene, x="contrast", y='LFC',
                           hover_data=[identifier,"library", "contrast", "fdr"], stripmode='overlay')
            fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                              font=dict(size=16))
            fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='LightGrey')
            col.plotly_chart(fig, use_container_width=True)


