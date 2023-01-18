import streamlit as st
from pathlib import Path
import pandas as pd
import numpy as np
import requests
from time import sleep
from scripts.graphs import define_color_scheme
from scripts.graphs import plot_rank, plot_position, plot_heatmap
import base64
import matplotlib
import seaborn as sns
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
import pandera as pa
from pandera.typing import Index, DataFrame, Series
from pandera.errors import SchemaError
from matplotlib import colors


class ResultSchema(pa.SchemaModel):
    LFC: Series[float] = pa.Field(coerce=True)
    neg_selection_fdr: Series[float] = pa.Field(coerce=True)
    pos_selection_fdr: Series[float] = pa.Field(coerce=True)
    contrast: Series[str] = pa.Field(coerce=True)
    library: Series[str] = pa.Field(coerce=True)


class ResultDataSet:
    def __init__(self, result_files=(), gene_id='locus_tag'):
        self.kind: str = 'result'
        self.result_files: str = result_files
        self.gene_id: str = gene_id
        self.results_df, self.gff_columns_found = self.combine_result_files() # todo figure out how to add typehints here
        self.subset_df = pd.DataFrame()
        self.string_df = pd.DataFrame()
        self.kegg_df = pd.DataFrame()

    # Load your final results
    def process_results_file(self, result_file):
        df = pd.read_csv(result_file)
        if 'library' not in df.columns:
            library_name = st.text_input("Add library name?", value=result_file.name)
            df['library'] = library_name
        gff_col_names = ['Chromosome', 'Start']
        gff_columns_found = all([c in df.columns for c in gff_col_names])
        if self.gene_id not in df.columns:
            st.markdown(f"""WARNING! No {self.gene_id} column found. Using {df.columns[0]} as gene names to display""")
            self.gene_id = df.columns[0]
        try:
            df = ResultSchema.validate(df)
        except SchemaError:
            st.write("RESULT FILE FORMAT ERROR 🙁")
            df = pd.DataFrame()
        return df, gff_columns_found

    def combine_result_files(self):
        if not self.result_files:
            return pd.DataFrame(), False
        else:
            results_df_list = []
            gff_col_list = []
            for uploaded_result in self.result_files:
                st.write(f"Processing {uploaded_result.name}")
                df, gff_columns_found = self.process_results_file(uploaded_result)
                results_df_list.append(df)
                gff_col_list.append(gff_columns_found)
            fdf = pd.concat(results_df_list)
            fdf['fdr'] = np.where(fdf['LFC'] < 0, fdf['neg_selection_fdr'], fdf['pos_selection_fdr'])
            fdf['log10FDR'] = -10 * np.log10(fdf['fdr'])
            fdf = fdf.fillna({self.gene_id: 'N/A'})  # todo do I need to do this?
            return fdf, all(gff_col_list)  # todo all or any?

    def summarize_libraries(self, library_to_show, lfc_low, lfc_hi, fdr_th):
        if not lfc_hi:
            self.results_df['hit'] = ((abs(self.results_df['LFC']) > lfc_low) & (self.results_df['fdr'] < fdr_th))
        else:
            self.results_df['hit'] = ((self.results_df['LFC'] > lfc_low)& (self.results_df['LFC'] < lfc_hi) & (self.results_df['fdr'] < fdr_th))
        if library_to_show != 'All':
            self.results_df = self.results_df[self.results_df['library'] == library_to_show]
            self.results_df['LFC_median'] = self.results_df['LFC']
            self.results_df['library_nunique'] = 1
            self.results_df['hit_sum'] = self.results_df['hit']
        else:
            df_grouped = self.results_df.groupby([self.gene_id, 'contrast']).agg({'LFC': ['median'],
                                                                                  'library': ['nunique'],
                                                                                  'hit': ['sum']}).reset_index()
            df_grouped.columns = [self.gene_id, 'contrast', 'LFC_new', 'library_nunique', 'hit_sum']
            self.results_df = self.results_df.merge(df_grouped, on=[self.gene_id, 'contrast'])
            self.results_df['LFC_median'] = self.results_df['LFC_new']

    def subset_results(self, contrast_to_show):
        #df = self.results_df[~fdf[gene_id].str.contains(':')].copy() # todo come up with cleverer way to filter these
        self.subset_df = self.results_df[(self.results_df['contrast'].isin(contrast_to_show))].copy()

    def connect_to_string(self, up, fdr_th, lfc_low, lfc_hi, species):
        if not lfc_hi:
            if up == 'Upregulated Only':
                self.string_df = self.subset_df[(self.subset_df['LFC'] > 0) & (self.subset_df['hit'] == True)]
            elif up == 'Downregulated Only':
                self.string_df = self.subset_df[(self.subset_df['LFC'] < 0) & (self.subset_df['hit'] == True)]
            else:
                self.string_df = self.subset_df[self.subset_df['hit'] == True]
            st.markdown(f"There are {self.string_df[self.gene_id].nunique()} hits with FDR < {round(fdr_th, 2)} and absolute LFC > {lfc_low}")
        else:
            self.string_df = self.subset_df[self.subset_df['hit'] == True]
            st.markdown(f"There are {self.string_df[self.gene_id].nunique()} hits with FDR < {round(fdr_th,2)} and within LFC range from {lfc_low} to {lfc_hi}")

        string_api_url = "https://version-11-5.string-db.org/api"
        output_format = 'tsv-no-header'
        method = 'get_link'
        my_genes = list(self.string_df[self.gene_id].unique())
        request_url = "/".join([string_api_url, output_format, method])
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

    def get_gene_to_kegg_map(self, kegg_id):
        if self.subset_df.contrast.nunique() > 1:
            st.write("⚠️WARNING: Multiple contrasts selected! The map will show median LFC across the selected contrasts")

        hitSummary = (self.subset_df.groupby([kegg_id]).hit.sum()
                      .reset_index()
                      .rename({'hit': 'hitSum'}, axis=1))

        self.subset_df = (self.subset_df.merge(hitSummary, on=kegg_id, how='outer'))
        self.subset_df['hitStar'] = self.subset_df['hitSum'].apply(lambda x: '*' if x > 0 else '')
        self.subset_df['NameForMap'] = self.subset_df[self.gene_id] + self.subset_df['hitStar'] + \
                                       "(" + self.subset_df['LFC_median'].round(2).astype(str) + ")"
        data_short = self.subset_df[set([self.gene_id, kegg_id, 'NameForMap', 'LFC_median'])].drop_duplicates().dropna()
        norm = matplotlib.colors.Normalize(vmin=-6, vmax=6, clip=True)
        mapper = matplotlib.cm.ScalarMappable(norm=norm, cmap=sns.diverging_palette(220, 20, as_cmap=True))
        data_short['hex'] = data_short.LFC_median.apply(mapper.to_rgba).apply(matplotlib.colors.to_hex)
        data_short['hex'] = data_short.hex.str.replace('#000000', "#7c7d83")
        ko_dict = data_short.set_index(kegg_id).to_dict()
        return ko_dict

    def get_kegg_path_df(self, pathGenes, kegg_id):
        if pathGenes:
            self.kegg_df = self.subset_df[self.subset_df[kegg_id].isin(pathGenes)].copy()
        else:
            self.kegg_df = self.subset_df.copy()


class DrawKeggMaps:
    def __init__(self, organism):
        self.organism = organism

    @st.cache
    def get_org_kegg_pathways(self):
        result = pd.read_table(io.StringIO(kegg_list("pathway", self.organism).read()), header=None)
        result.columns = [f'KEGG_Pathway', 'Pathway_Description']
        result[f'KEGG_Pathway'] = result[f'KEGG_Pathway'].str.split(":").str.get(1)
        result['Pathway_Description'] = result['Pathway_Description'].str.split(" - ").str.get(0)
        result['KEGG_Display'] = result[f'KEGG_Pathway'] + ":" + result['Pathway_Description']
        path_map = result.set_index('KEGG_Display').to_dict()
        return path_map['KEGG_Pathway']

    def displayPDF(self, file):
        # Opening file from file path
        with open(file, "rb") as f:
            base64_pdf = base64.b64encode(f.read()).decode('utf-8')
        # Embedding PDF in HTML
        pdf_display = F'<iframe src="data:application/pdf;base64,{base64_pdf}" width="700" height="1000" type="application/pdf"></iframe>'
        # Displaying File
        st.markdown(pdf_display, unsafe_allow_html=True)

    def display_kegg_map(self, pathwayName, ko_dict, title, numeric=False):
        """
        :
        """
        pathwayKGML = KGML_parser.read(kegg_get(pathwayName, "kgml"))


        pathGeneNames = [gene.name.split() for gene in pathwayKGML.genes]
        # todo add gene parsing function
        # pathGeneNames = set([extract_gene_name(gene, numeric=numeric) for sublist in pathGeneNames for gene in sublist])
        pathGeneNames = set([gene.split(":")[1] for sublist in pathGeneNames for gene in sublist]) # todo add optional number parsing
        canvas = KGMLCanvas(pathwayKGML, import_imagemap=True)
        not_found = []
        for element in pathwayKGML.genes:
            color = None
            name = None
            node_kos = [e.split(":")[1] for e in element.name.split()]
            for ko in node_kos:
                color = ko_dict['hex'].get(ko, color)
                name = ko_dict['NameForMap'].get(ko, name)
            for graphic in element.graphics:
                if color is not None:
                    graphic.bgcolor = color
                    graphic.name = name
                    not_found.append(0)
                else:
                    not_found.append(1)
        if sum(not_found)/len(not_found) > 0.85:
            st.warning(f'⚠️ {sum(not_found)} out of {len(not_found)} pathway genes not found in the dataset. Double check gene names match those used by KEGG')

        fname = f"{title}_map.pdf"
        canvas.draw(fname)
        k1, k2 = st.columns(2)
        st.info("Display works in Firefox only")
        if k1.button(f'Display {pathwayName} map'):
            self.displayPDF(fname)
        with open(fname, "rb") as f:
            k2.download_button(
                f"Download {pathwayName} map",
                data=f,
                file_name=fname,
            )
        return pathGeneNames

@st.cache
def convert_df(df):
    return df.to_csv(index=False).encode('utf-8')


#############
#    APP    #
#############

def app():
    st.markdown(""" # Differential Abundance """)

    # with st.expander('How this works: '):
    #     st.markdown("""
    #
    #     ### ADD TEXT HERE
    #
    #     - ADD TEXT HERE
    #     #
    #     """)

    colors, alphabetClrs, all_clrs = define_color_scheme()
    with st.container():
        st.subheader('Load the data file')
        data_type = st.radio('Choose dataset to show', ['Load my data', 'Look at an example'], index=0)
        if data_type == 'Load my data':
            result_files = st.file_uploader('Upload results file', accept_multiple_files=True)
            gene_id = st.text_input('Unique gene identifier used in the result files', value='Name')
        else:
            st.subheader('Example results file')
            result_files = [Path("examples/example_rra_results.csv")]
            gene_id = 'Name'
        rds = ResultDataSet(result_files, gene_id)
        show_sample = st.checkbox('Show sample of the dataset')
        if show_sample:
            try:
                st.write(rds.results_df.sample(5))
            except ValueError:
                st.write('Result table is empty')

    if not rds.results_df.empty:

        contrasts = rds.results_df['contrast'].sort_values().unique()
        libraries = rds.results_df['library'].sort_values().unique()
        libraries = ['All'] + list(libraries) if len(libraries) > 1 else libraries
        st.subheader("Fitness Results")
        # SUBSET DATAFRAME TO SPECIFIC CONTRAST AND LIBRARY
        contrast_col, lfc_col, fdr_col, lfc_lib_col = st.columns(4)
        contrast_to_show = contrast_col.multiselect('Select a contrast', contrasts, default=contrasts[0])
        library_to_show = lfc_lib_col.selectbox('Select library to show', libraries)
        fdr_th = fdr_col.number_input('FDR cutoff', value=0.05)
        type_lfc_th = lfc_col.radio('Absolute LFC cutoff or define range', ['Absolute', 'Range'])
        if type_lfc_th == 'Absolute':
            lfc_low = lfc_col.number_input('Log FC cutoff (absolute)', min_value=0.0, step=0.5, value=1.0)
            lfc_hi = None
        else:
            lfc_low = lfc_col.number_input('Min Log FC',  step=0.5, value=-5.0)
            lfc_hi = lfc_col.number_input('Max Log FC',  step=0.5, value=-1.0)

        rds.summarize_libraries(library_to_show,  lfc_low, lfc_hi, fdr_th)
        rds.subset_results(contrast_to_show)
        # SUBMIT SUBSET TO STRING
        st.markdown('#### Analyze with STRING-db')
        with st.expander('STRING-db'):
            st.info(f"❗Make sure STRING recognizes unique gene identifier (you've entered `{gene_id}`) for the taxon you specify")
            species = st.number_input("NCBI species taxid", value=99287, help='Salmonella Typhimurium: 99287')
            st.markdown(f"Analyze hits for ``{', '.join(contrast_to_show)}`` contrast for ``{library_to_show}`` library(ies)")
            if not lfc_hi:
                up = st.radio('Up or Down?', ('Upregulated Only', 'Downregulated Only', 'Both'), key='string')
            else:
                up = 'NA'
            rds.connect_to_string(up, fdr_th, lfc_low, lfc_hi, species)
        # SUBSET TO A SPECIFIC KEGG PATHWAY IF DESIRED/POSSIBLE
        st.markdown('#### KEGG Maps')
        st.write('')
        kegg_avail = st.checkbox('KEGG annotation available?')
        if kegg_avail:
            org_col, kegg_col = st.columns(2)
            organismId = org_col.text_input('Enter 3 letter organsim code', value='sey')
            kegg_options = [c for c in rds.results_df.columns if 'LFC' not in c and 'fdr' not in c]
            try:
                kix = kegg_options.index('locus_tag')
            except ValueError:
                kix = 0
            kegg_id = kegg_col.selectbox('Column corresponding to KEGG Entry names (usually locus_tag)', options=kegg_options, index=kix)
            ko_dict = rds.get_gene_to_kegg_map(kegg_id)
            km = DrawKeggMaps(organismId)
            with st.spinner(f"Loading the list of all KEGG pathways for {organismId}"):
                pathwayMap = km.get_org_kegg_pathways()

            pathwayDescription = st.selectbox('Select KEGG Pathway to explore', pathwayMap.keys())
            pathwayName = pathwayMap[pathwayDescription]

            # todo add checkbox (by default unchecked) asking user whether to display locus # only
            # if checkbox is checked, numeric = True, otherwise numeric = False
            #locus = st.checkbox("Locus Number")

            with st.spinner(f'Loading KEGG map for {pathwayName}...'):
                pathwayGenes = km.display_kegg_map(pathwayName, ko_dict,
                                                   f"{pathwayName}-{'-'.join(contrast_to_show)}") # add numeric
        else:
            st.write("KEGG maps are unavailable.")
            kegg_id = gene_id
            pathwayGenes = []
            pathwayName = 'All Genes'

        rds.get_kegg_path_df(pathwayGenes, kegg_id)

        # VISUALIZE SUBSET BY RANK OR POSITION
        st.markdown('#### Fitness results by rank or position')
        with st.expander('Fitness Results'):
            if rds.gff_columns_found:
                graph_type = st.radio('By rank or by position?', ('Rank', 'Position'))
            else:
                graph_type = 'Rank'
            # DISPLAY RESULTS BY RANK
            st.subheader(f"{','.join(contrast_to_show)} - {pathwayName}")
            if graph_type == 'Rank':
                rank_df = (rds.kegg_df[[gene_id, 'contrast', 'LFC_median', 'hit']].drop_duplicates()
                           .sort_values('LFC_median')
                           .reset_index()
                           .reset_index()
                           .rename({'level_0': 'ranking'}, axis=1)
                           .sort_values('ranking'))
                hover_data = st.multiselect('Data to show on hover:', rank_df.columns, key='rank')
                hover_dict = {s: True for s in hover_data}
                fig = plot_rank(rank_df.dropna(subset=[gene_id, 'LFC_median']), colors, hover_dict, gene_id)
                st.plotly_chart(fig, use_container_width=True)
            # DISPLAY RESULTS BY POSITION
            else:
                chr_to_show = st.selectbox('Select chromosome', rds.kegg_df.Chromosome.unique())
                position_df = rds.kegg_df[rds.kegg_df.Chromosome == chr_to_show].copy()
                hover_data = st.multiselect('Data to show on hover:', position_df.columns, key='position')
                hover_dict = {s: True for s in hover_data}
                fig = plot_position(position_df, hover_dict, gene_id)
                st.plotly_chart(fig, use_container_width=True)
        # DRAW HEATMAPS FOR THE SUBSET
        st.markdown('#### Heatmaps for genes/pathways of interest')
        with st.expander('Show LFC Heatmaps'):
            c1, c2 = st.columns(2)
            heat_by = c1.radio('Show heatmap by:', ( 'KEGG Pathway', 'Genes of Interest'))
            if heat_by == 'KEGG Pathway':
                if pathwayName == 'All Genes':
                    st.markdown('No KEGG pathway selected :thinking_face:')
                    heat_df = pd.DataFrame() # todo stopping here. Tired.
                else:
                    heatmap_title = pathwayName
                    heat_df = rds.results_df[rds.results_df[kegg_id].isin(pathwayGenes)]
                    absent = pd.DataFrame(
                        pd.Series(list(set(pathwayGenes) - set(heat_df[kegg_id].unique())), name=kegg_id))
                    heat_df = heat_df[[gene_id, kegg_id, 'LFC_median', 'contrast']].drop_duplicates()
                    heat_df = pd.concat([heat_df, absent], axis=0)
                    s = heat_df[gene_id].fillna('-').values
                    p = heat_df[kegg_id].values
                    heat_df['full_id'] = [
                        f"{kegg_id}: {gene_id}" if gene_id != '-' and gene_id != kegg_id else f"{kegg_id}" for
                        kegg_id, gene_id in zip(p, s)]
                    heat_df = heat_df.pivot(index='full_id', columns='contrast', values='LFC_median')
            elif heat_by == 'Genes of Interest':
                heatmap_title = ''
                show_genes = st.multiselect("Choose gene(s) of interest", rds.results_df[gene_id].unique())
                if not show_genes:
                    st.markdown("No genes selected :thinking_face:")
                    heat_df = pd.DataFrame()
                else:
                    heat_df = rds.results_df[rds.results_df[gene_id].isin(show_genes)][[gene_id, 'contrast', 'LFC_median']].drop_duplicates()
                    heat_df = heat_df.pivot(index=gene_id, columns='contrast', values='LFC_median')
            if not heat_df.empty:
                #st.write(heatDf) # todo replace by aggrid

                fig = plot_heatmap(heat_df)
                st.subheader(heatmap_title)
                st.plotly_chart(fig, use_container_width=False)
                st.markdown("### Download Files as csv")
                fname = st.text_input("File name", value=f"heatmap")
                fname = fname + ".csv"
                st.download_button("Download data as csv file", convert_df(heat_df.reset_index()), file_name=fname)
