import streamlit as st
import pathlib
from pathlib import Path
import pandas as pd
import numpy as np
import plotly.express as px
from itertools import cycle
import requests
from time import sleep
from io import StringIO

from scripts.graphs import show_lfc_ranking, link_to_string

def process_library_results(result_file, file_type='mageck'):
    df = pd.read_csv(result_file)
    if 'library' not in df.columns:
        library_name = st.text_input("Add library name?", value=result_file.name)
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


@st.cache
def convert_df(df):
    return df.to_csv(index=False).encode('utf-8')


@st.cache
def process_gmt(gmt_file, identifier='Name', path_name=1):
    genes = []
    paths = []
    gmt_lines = StringIO(gmt_file.getvalue().decode("utf-8")).read().split("\n")
    for line in gmt_lines:
        for gene in line.split('\t')[path_name:]:
            genes.append(gene.strip("\n"))
            paths.append(":".join(line.split('\t')[0:path_name]))
    return pd.DataFrame([genes, paths], index=[identifier, 'KEGG_Pathway']).T


def process_gmt_file(gmt_file, path_name=1):
    genes = []
    paths = []
    with open(gmt_file, 'r', errors='ignore') as fh:
        for line in fh.readlines():
            for gene in line.split('\t')[path_name:]:
                genes.append(gene.strip("\n"))
                paths.append(":".join(line.split('\t')[0:path_name]))
    return pd.DataFrame([genes, paths], index=['Name', 'KEGG_Pathway']).T


def app():
    url = "https://doi.org/10.1186/s13059-014-0554-4"
    st.markdown(f""" # Differential Gene Fitness""")

    with st.expander("How this works: "):
        st.markdown(f""" ### Visualize gene log fold changes (LFCs) between two conditions of interest

    - Statistical analysis is performed using [MAGeCK]({url}). MAGeCK uses a modified version of Robust Ranking Aggregation algorithm to identify negativly and positively selected genes.

    - `gmt` file is needed to add pathway visualization

    - The central graph displays LFC for each gene ranked from lowest to highest, the size of the circle is inversely proportional to the FDR

    - The heatmaps show LFC for all the genes for a specific KEGG pathway for all the contrasts (comparisions)

        """)

    with st.container():
        st.subheader('Browse an example data set')
        data_type = st.radio('Choose dataset to show', ['Load my data', 'Look at an example'], index=0)
        if data_type == 'Load my data':
            result_files = st.file_uploader('Upload results file', accept_multiple_files=True)
            gmt_file = st.file_uploader('Upload gmt file')
            identifier = st.text_input('Gene attribute used in the result files', value='Name')
        else:
            result_files = [Path("examples/example_rra_results.csv")]
            gmt_file = Path("examples/04-03-2022-SL1344-KEGG-API.gmt")
            identifier = 'Name'
            st.subheader('Example results file')
            ex_df = (pd.read_csv(result_files[0]))
            ex_df = ex_df[ex_df.contrast == 'd1'].sort_values('number_of_barcodes', ascending=False)
            st.write(ex_df.head(5))
        if len(result_files) < 1:
            st.stop()
        result_dfs = []
        for uploaded_result in result_files:
            st.write(f"Processing {uploaded_result.name}")
            df = process_library_results(uploaded_result)
            result_dfs.append(df)
        fdf = pd.concat(result_dfs)
        if gmt_file:
            if type(gmt_file) is pathlib.PosixPath:
                gmt_df = process_gmt_file(gmt_file, path_name=2)
            else:
                gmt_df = process_gmt(gmt_file, identifier, path_name=2)
            fdf = fdf.merge(gmt_df, how='left', on=identifier)
        else:
            if "KEGG_Pathway" not in list(fdf.columns):
                fdf['KEGG_Pathway'] = np.nan

        fdf['fdr'] = np.where(fdf['LFC'] < 0, fdf['neg_selection_fdr'], fdf['pos_selection_fdr'])
        fdf['log10FDR'] = -10 * np.log10(fdf['fdr'])
        contrasts = fdf['contrast'].sort_values().unique()
        libraries = fdf['library'].unique()

    st.subheader("Results by log fold change")
    with st.expander('LFC rankings'):
        fig, df = show_lfc_ranking(fdf, contrasts, libraries)
        c1, c2 = st.columns((2,1))
        c1.plotly_chart(fig, use_container_width=True)
        link_to_string(df, c2, lfc_col='LFC', gene_name='Name')


    st.subheader("Pathway Heatmaps")
    with st.expander('Show LFC Heatmaps'):
        c1, c2 = st.columns(2)
        heat_by = st.radio('Show heatmap by:', ('Genes of Interest', 'KEGG Pathway'))
        if heat_by == 'KEGG Pathway':
            show_kegg_heat1 = st.selectbox('Show KEGG Pathway', list(fdf.KEGG_Pathway.dropna().unique()), key='leftKegg')
            if not show_kegg_heat1:
                st.markdown('No KEGG Annotation found :thinking_face:')
            else:
                heatDfLeft = (fdf[fdf.KEGG_Pathway == show_kegg_heat1][['Name', 'LFC', 'contrast', 'library']]
                              .groupby(['contrast', 'Name']).LFC.median().reset_index())
                heatDfLeft.columns = ['contrast', 'Name', 'LFC']
                heatDfLeft = heatDfLeft.pivot(index="Name", columns='contrast', values='LFC')

                fig2 = px.imshow(heatDfLeft, color_continuous_scale=px.colors.diverging.Geyser,
                                 color_continuous_midpoint=0,
                                 width=600, height=900)
                fig2.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                                   font=dict(size=22))
                st.plotly_chart(fig2, use_container_width=True)
        elif heat_by == 'Genes of Interest':
            show_genes = st.multiselect("Choose gene(s) of interest", fdf[identifier].unique())
            if not show_genes:
                st.markdown("No genes selected :thinking_face:")
            else:
                heatDfRight = (fdf[fdf[identifier].isin(show_genes)][['Name', 'LFC', 'contrast', 'library']]
                               .groupby(['contrast', 'Name']).LFC.mean().reset_index())
                heatDfRight.columns = ['contrast', 'Name', 'LFC']
                heatDfRight = heatDfRight.pivot(index="Name", columns='contrast', values='LFC')
                fig2 = px.imshow(heatDfRight, color_continuous_scale=px.colors.diverging.Geyser,
                                 color_continuous_midpoint=0,
                                 width=600, height=900)
                fig2.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                                   font=dict(size=22))
                st.plotly_chart(fig2, use_container_width=True)


    # st.subheader("Protein-protein interactions")
    # with st.expander('Tabular Results'):
    #     contCol, libCol = st.columns(2)
    #     compContrasts = contCol.multiselect('Select contrasts to display', contrasts)
    #     libs = libCol.multiselect('Select libraries to display', ['All'] + list(libraries), default='All')
    #     if 'All' in libs:
    #         vennDf = fdf.copy()
    #     else:
    #         vennDf = fdf[fdf.library.isin(libs)].copy()
    #     if not compContrasts:
    #         st.write("No contrast selected")
    #     else:
    #         c1, c2 = st.columns(2)
    #         filters = {}
    #         for col, contrast in zip(cycle([c1, c2]), compContrasts):
    #             col.write(contrast)
    #             l = col.number_input('LFC cutoff', value=-5.0, step=0.5, key=f'{contrast}_lfc')
    #             f = col.number_input('FDR cutoff', value=0.05, step=0.01, key=f'{contrast}_fdr')
    #             filters[contrast] = (l, f)
    #         comp_genes = []
    #         for key, value in filters.items():
    #             if value[0] > 0:
    #                 genes = set(vennDf[(vennDf.contrast == key) & (vennDf['LFC'] > value[0]) & (vennDf.fdr < value[1])][
    #                                 identifier].values)
    #             else:
    #                 genes = set(vennDf[(vennDf.contrast == key) & (vennDf['LFC'] < value[0]) & (vennDf.fdr < value[1])][
    #                                 identifier].values)
    #             comp_genes.append(genes)
    #         intersect_genes = set.intersection(*comp_genes)
    #         vennDf = vennDf[vennDf[identifier].isin(intersect_genes)].copy()
    #         vennDf = vennDf[[identifier, 'library', 'LFC', 'fdr', 'contrast']].drop_duplicates()
    #         # vennDf2 = vennDf.pivot(index='Name', columns='contrast', values=['LFC', 'fdr'])
    #         st.write(vennDf)
    #
    #         string_col, download_col = st.columns(2)
    #
    #         string_api_url = "https://version-11-5.string-db.org/api"
    #         output_format = 'tsv-no-header'
    #         method = 'get_link'
    #         my_genes = set(vennDf[identifier].values)
    #         request_url = "/".join([string_api_url, output_format, method])
    #         string_col.markdown("### STRING Interaction Network")
    #         species = string_col.number_input("NCBI species taxid", value=99287, help='Salmonella Typhimurium: 99287')
    #
    #         params = {
    #             "identifiers": "\r".join(my_genes),  # your protein
    #             "species": species,  # species NCBI identifier
    #             "network_flavor": "confidence",  # show confidence links
    #             "caller_identity": "mbarq"  # your app name
    #         }
    #
    #         if string_col.button('Get STRING network'):
    #             results = requests.post(request_url, data=params)
    #             network_url = results.text.strip()
    #             st.markdown(f"[Link to STRING network]({network_url})")
    #             sleep(1)
    #         download_col.markdown("### Download Files as csv")
    #         fname = download_col.text_input("File name", value="mbarq_results")
    #         fname = fname + ".csv"
    #         download_col.download_button("Download data as csv file", convert_df(vennDf), file_name=fname)

    # st.subheader("Gene Selector")
    # with st.expander('Gene Selector'):
    #     genes = st.multiselect("Choose gene(s) of interest", fdf[identifier].unique())
    #     if not genes:
    #         st.stop()
    #     c3, c4 = st.columns(2)
    #
    #     for col, gene in zip(cycle([c3, c4]), genes):
    #         gene_df = fdf[fdf[identifier] == gene].copy()
    #         gene_df = gene_df[[identifier, 'library', 'contrast', 'LFC', 'fdr']].drop_duplicates()
    #         gene_df = gene_df.sort_values('contrast')
    #         fig = px.strip(gene_df, title=gene, x="contrast", y='LFC', color='library',
    #                        hover_data=[identifier, "library", "contrast", "fdr"], stripmode='overlay')
    #         fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
    #                           font=dict(size=16))
    #         fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='LightGrey')
    #         fig.add_hline(y=0, line_width=2, line_dash="dash", line_color="grey")
    #         col.plotly_chart(fig, use_container_width=True)


