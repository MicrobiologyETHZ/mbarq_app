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
from scripts.graphs import define_color_scheme
from scripts.graphs import plot_rank, plot_position, plot_heatmap


# Load your final results
def process_library_results(result_file, file_type='mbarq1', gene_id='Name'):
    df = pd.read_csv(result_file)
    if 'library' not in df.columns:
        library_name = st.text_input("Add library name?", value=result_file.name)
        df['library'] = library_name
    fixed_col_names = {'mbarq1': ['library', 'LFC',
                                  'neg_selection_fdr', 'pos_selection_fdr', 'contrast']}
    kegg_col_names = {'mbarq1': ['KEGG_Pathway']}
    position_col_names = {'mbarq1': ['Chromosome', 'Start', 'End']}
    missing_cols = [c for c in fixed_col_names[file_type] if c not in df.columns]
    if len(missing_cols) > 0:
        st.markdown(
            f"""The following columns are missing from the results files: {', '.join(missing_cols)}. 
            Please rename the columns/rerun mBARq and try again. Skipping {result_file.name}""")
        return pd.DataFrame()
    annot_columns = [c for c in df.columns if c in kegg_col_names[file_type]]
    if len(annot_columns) > 0:
        st.markdown(
            f"""The following KEGG annotation columns found: {', '.join(annot_columns)}"""
        )
        df[annot_columns] = df[annot_columns].replace(np.nan, '')
    else:
        st.markdown("""No KEGG annotation found. See ** tutorial on how to add KEGG annotation to your results """)
    gff_columns_found = all([c in df.columns for c in position_col_names[file_type]])
    st.markdown(f"""Gene Start/End Positions found: {gff_columns_found}""")
    if gene_id not in df.columns:
        st.markdown(f"""WARNING! No {gene_id} column found. Using {df.columns[0]} as gene names to display""")
        gene_id = df.columns[0]
    return df, gene_id


def subset_results(fdf, gene_id,  contrast_to_show, library_to_show, fdr_th, lfc_th):
    df = fdf[~fdf[gene_id].str.contains(':')].copy() # todo come up with cleverer way to filter these
    df = df[(df['contrast'].isin(contrast_to_show)) | (df['contrast'].isna())]
    if library_to_show != 'All':
        df = df[(df['library'] == library_to_show) | (df['library'].isna())]
    df['hit'] = ((abs(df['LFC']) > lfc_th) & (df['fdr'] < fdr_th))
    return df


def connect_to_string(hits_df, gene_name=''):
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = 'tsv-no-header'
    method = 'get_link'
    my_genes = list(hits_df[gene_name].unique()) if gene_name else list(hits_df.index)
    request_url = "/".join([string_api_url, output_format, method])
    species = st.number_input("NCBI species taxid", value=99287, help='Salmonella Typhimurium: 99287')
    params = {
        "identifiers": "\r".join(my_genes),  # your protein
        "species": species,  # species NCBI identifier
        "network_flavor": "confidence",  # show confidence links
        "caller_identity": "explodata"  # your app name
    }
    #
    if st.button('Get STRING network'):
        network = requests.post(request_url, data=params)
        network_url = network.text.strip()
        st.markdown(f"[Link to STRING network]({network_url})")
        sleep(1)


def display_kegg_selector(df, kegg_col='KEGG_Pathway', kegg_names_file="./examples/20-10-22-kegg-pathway-list-ko.csv"):
    if kegg_col not in df.columns:
        st.write("No KEGG Pathway information found.")
        return 'None'
    avail_pathways = [k.split(',') for k in df[kegg_col].unique()]
    avail_pathways = set([s for k in avail_pathways for s in k if s.startswith('ko')])
    kegg_names = pd.read_csv(kegg_names_file, index_col=0, header=None).to_dict()[1]
    avail_pathways = [f"{i}:{kegg_names[i]}" if i in kegg_names.keys() else i for i in avail_pathways]
    return st.selectbox('Choose KEGG Pathway to display', ["None"] + avail_pathways)


@st.cache
def convert_df(df):
    return df.to_csv(index=False).encode('utf-8')


#############
#    APP    #
#############

def app():
    colors, alphabetClrs, all_clrs = define_color_scheme()

    with st.container():
        st.subheader('Browse an example data set')
        data_type = st.radio('Choose dataset to show', ['Load my data', 'Look at an example'], index=0)
        if data_type == 'Load my data':
            result_files = st.file_uploader('Upload results file', accept_multiple_files=True)
            gene_id = st.text_input('Unique gene identifier used in the result files', value='Name')
            if len(result_files) > 0:
                result_dfs = []
                for uploaded_result in result_files:
                    st.write(f"Processing {uploaded_result.name}")
                    df, gene_id = process_library_results(uploaded_result, 'mbarq1', gene_id)
                    result_dfs.append(df)
                fdf = pd.concat(result_dfs)
            else:
                fdf = pd.DataFrame()
        else:
            result_files = [Path("examples/example_rra_results.csv")]
            gene_id = 'Name'
            st.subheader('Example results file')
            fdf, gene_id = process_library_results(result_files[0], 'mbarq1', gene_id)
            fdf = fdf[fdf.contrast == 'd1'].sort_values('number_of_barcodes', ascending=False)
            st.write(fdf.head(5))

    if not fdf.empty:

        fdf['fdr'] = np.where(fdf['LFC'] < 0, fdf['neg_selection_fdr'], fdf['pos_selection_fdr'])
        fdf['log10FDR'] = -10 * np.log10(fdf['fdr'])
        fdf = fdf.fillna({gene_id: 'N/A'})
        contrasts = fdf['contrast'].dropna().sort_values().unique()
        libraries = fdf['library'].dropna().sort_values().unique()
        st.subheader("Fitness Results")
        contrast_col, lfc_col, fdr_col, lfc_lib_col = st.columns(4)
        contrast_to_show = contrast_col.multiselect('Select a contrast', contrasts)
        library_to_show = lfc_lib_col.selectbox('Select library to show', ["All"] + list(libraries))
        fdr_th = fdr_col.number_input('FDR cutoff', value=0.05)
        lfc_th = lfc_col.number_input('Log FC cutoff (absolute)', min_value=0.0, step=0.5, value=1.0)

        # SUBSET DATAFRAME TO SPECIFIC CONTRAST AND LIBRARY
        subset_df = subset_results(fdf, gene_id, contrast_to_show, library_to_show, fdr_th, lfc_th)

        # SUBMIT SUBSET TO STRING
        st.markdown('#### Analyze with STRING-db')
        with st.expander('STRING-db'):
            st.markdown(f"Analyze hits for ``{', '.join(contrast_to_show)}`` contrast for ``{library_to_show}`` library(ies)")
            up = st.radio('Up or Down?', ('Upregulated Only', 'Downregulated Only', 'Both'), key='string')
            if up == 'Upregulated Only':
                string_df = subset_df[(subset_df['LFC'] > 0) & (subset_df['hit'] == True)]
            elif up == 'Downregulated Only':
                string_df = subset_df[(subset_df['LFC'] < 0) & (subset_df['hit'] == True)]
            st.markdown(f"There are {string_df[gene_id].nunique()} hits with FDR < {fdr_th} and absolute LFC > {lfc_th}")
            connect_to_string(string_df, gene_name=gene_id)

        # LOOK AT SUBSET BY RANK OR POSITION
        st.markdown('#### Rank Fitness Results')

        # SUBSET TO A SPECIFIC KEGG PATHWAY IF DESIRED
        kegg_col = 'KEGG_Pathway'
        kegg_to_show = display_kegg_selector(subset_df, kegg_col)
        if kegg_to_show != 'None':
            kegg_df = subset_df[subset_df[kegg_col].str.contains(kegg_to_show.split(":")[0])]
        else:
            kegg_df = subset_df.copy()
        # VISUALIZE SUBSET BY RANK OR POSITION
        with st.expander('Fitness Results'):
            if all([c in subset_df.columns for c in ['Chromosome', 'Start']]):
                graph_type = st.radio('By rank or by position?', ('Rank', 'Position'))
            else:
                graph_type = 'Rank'

            # DISPLAY RESULTS BY RANK
            st.subheader(f"{','.join(contrast_to_show)} - {kegg_to_show}")
            if graph_type == 'Rank':
                rank_df = kegg_df.sort_values('LFC').reset_index().reset_index().rename({'level_0': 'ranking'}, axis=1)
                hover_data = st.multiselect('Data to show on hover:', rank_df.columns, key='rank')
                hover_dict = {s: True for s in hover_data}
                fig = plot_rank(rank_df.dropna(subset=[gene_id, 'LFC']), colors, hover_dict)
                st.plotly_chart(fig, use_container_width=True)
            # DISPLAY RESULTS BY POSITION
            else:
                chr_to_show = st.selectbox('Select chromosome', kegg_df.Chromosome.unique())
                position_df = kegg_df[kegg_df.Chromosome == chr_to_show].copy()
                position_df = position_df[[gene_id, 'Start', 'LFC', 'contrast', 'library']]
                position_df = (position_df.groupby([gene_id, 'Start', 'contrast'])
                               .agg({'LFC': ['mean'], 'library': ['nunique']})
                               .reset_index())
                position_df.columns = [gene_id, 'Start', 'contrast', 'mean_LFC', 'number of libraries']
                position_df = (position_df.merge(kegg_df.drop(['Start'], axis=1), on=[gene_id, 'contrast'], how='left')
                               .dropna(subset=[gene_id, 'LFC']))
                hover_data = st.multiselect('Data to show on hover:', position_df.columns, key='position')
                hover_dict = {s: True for s in hover_data}
                fig = plot_position(position_df, colors, hover_dict)
                st.plotly_chart(fig, use_container_width=True)

        # DRAW HEATMAPS FOR THE SUBSET
        with st.expander('Show LFC Heatmaps'):
            c1, c2 = st.columns(2)
            heatDf = pd.DataFrame()
            heat_by = c1.radio('Show heatmap by:', ('Genes of Interest', 'KEGG Pathway'))
            if heat_by == 'KEGG Pathway':
                if not kegg_to_show or kegg_to_show == 'None':
                    st.markdown('No KEGG Annotation found :thinking_face:')
                else:
                    heatDf = kegg_df.copy()
            elif heat_by == 'Genes of Interest':
                show_genes = st.multiselect("Choose gene(s) of interest", subset_df[gene_id].unique())
                if not show_genes:
                    st.markdown("No genes selected :thinking_face:")
                else:
                    heatDf = subset_df[subset_df[gene_id].isin(show_genes)]
            if not heatDf.empty:
                #st.write(heatDf) # todo replace by aggrid
                absent = heatDf[heatDf.LFC.isna()]
                heatDf = (heatDf[[gene_id, 'LFC', 'contrast', 'library']]
                          .groupby(['contrast', 'Name'])
                          .agg({'LFC': ['mean']})
                          .reset_index())
                heatDf.columns = ['contrast', 'Name', 'LFC (mean)']
                heatDf = pd.concat([heatDf, absent[['Name']]]).drop_duplicates()
                heatDf = heatDf.pivot(index=gene_id, columns='contrast', values='LFC (mean)')

                fig = plot_heatmap(heatDf)
                st.plotly_chart(fig, use_container_width=False)
                st.markdown("### Download Files as csv")
                fname = st.text_input("File name", value=f"heatmap")
                fname = fname + ".csv"
                st.download_button("Download data as csv file", convert_df(heatDf.reset_index()), file_name=fname)
