from pathlib import Path
from typing import List, Union
import pandas as pd
import pandera as pa
import plotly.express as px
import streamlit as st
import yaml
from pandera.errors import SchemaError
from scripts.graphs import define_color_scheme
import numpy as np
from sklearn.decomposition import PCA
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas

import base64
import matplotlib
import seaborn as sns

import re

@st.cache
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv(index=False).encode('utf-8')

def define_color_scheme():
    alphabetClrs = px.colors.qualitative.Dark24
    clrs = ["#f7ba65", "#bf4713", "#9c002f", "#d73d00", "#008080", "#004c4c"]
    colors = {'grey': "#E2E2E2",
              'light_yellow': clrs[0],
              'darko': clrs[1],
              'maroon': clrs[2],
              'brighto': clrs[3],
              'teal': clrs[4],
              'darkteal': clrs[5]
              }
    sushi_colors = {'red': '#C0504D',
                    'orange': '#F79646',
                    'medSea': '#4BACC6',
                    'black': '#000000',
                    'dgreen': '#00B04E',
                    'lgreen': '#92D050',
                    'dblue': '#366092',
                    'lblue': '#95B3D7'}
   # all_clrs = [colors['brighto'], colors['teal'], colors['maroon']] + alphabetClrs[13:]
    all_clrs = ['#F79646', '#366092', '#00B04E', '#C0504D', colors['maroon'], colors['teal']] + alphabetClrs
    return colors, alphabetClrs, sushi_colors, all_clrs


class LibraryMap:

    def __init__(self, map_files: List = (),
                 map_df: pd.DataFrame = pd.DataFrame(),
                 attributes: List = (),
                 config_file: str = 'scripts/config.yaml'):
        self.map_files = map_files
        self.lib_map = map_df
        self.attributes = attributes
        self.color_by_cols = ('in CDS', 'library')
        self.stats = pd.DataFrame
        with open(config_file, 'r') as cf:
            config = yaml.load(cf, Loader=yaml.SafeLoader)['library_map']
        # Load column naming schema
        col_name_config = config['fixed_column_names']
        self.fixed_column_names = list(col_name_config.values())
        self.chr_col = col_name_config['chr_col']
        self.insertion_site_col = col_name_config['insertion_site_col']
        self.abundance_col = col_name_config['abundance_col']
        self.barcode_col = col_name_config['barcode_col']
        self.distance_col = col_name_config['distance_col']

    def load_map(self):
        map_dfs = []
        attr_names = []
        if self.map_files:
            for uploaded_map in self.map_files:
                st.write(f"Processing {uploaded_map.name}")
                df, df_name = pd.read_csv(uploaded_map), uploaded_map.name
                if 'library' not in df.columns:
                    library_name = st.text_input("Change library name?", value=df_name)
                    df['library'] = library_name
                missing_cols = [c for c in self.fixed_column_names if c not in df.columns]
                if len(missing_cols) > 0:
                    st.markdown(
                        f"""⚠️ The following columns are missing from the map files: {', '.join(missing_cols)}. 
                                    Please rename the columns/rerun mBARq and try again. Skipping {df_name} ⚠️""")
                    continue
                map_dfs.append(df)
            try:
                self.lib_map = pd.concat(map_dfs)
            except ValueError:
                st.error("No library map loaded")
        if not self.lib_map.empty:
            self.lib_map['library'] = 'example'
            self.lib_map['in CDS'] = self.lib_map[self.distance_col] == 0
            self.attributes = [c for c in self.lib_map.columns if c not in self.fixed_column_names]

    def validate_lib_map(self):
        lib_schema = pa.DataFrameSchema({
            self.chr_col: pa.Column(str, coerce=True),
            self.insertion_site_col: pa.Column(int, pa.Check(lambda x: x >= 0)),
            self.barcode_col: pa.Column(str, coerce=True),
            self.abundance_col: pa.Column(int, pa.Check(lambda x: x >= 0)),
            self.distance_col: pa.Column(float, nullable=True),
            'in CDS': pa.Column(bool),
            'library': pa.Column(str, coerce=True)
        }
        )
        try:
            self.lib_map = lib_schema.validate(self.lib_map)

        except SchemaError as err:
            st.error(f"""Schema Error: {err.args[0]}""")
            self.lib_map = pd.DataFrame()

    def graph_coverage_hist(self, chr_col_choice, num_bins, hist_col):
        df_to_show = self.lib_map[self.lib_map[self.chr_col] == chr_col_choice].sort_values(self.insertion_site_col)
        fig = px.histogram(df_to_show, x=self.insertion_site_col, nbins=int(num_bins),
                           labels={self.insertion_site_col: 'Position, bp'}, color_discrete_sequence=[hist_col])
        fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'},
                          autosize=True,
                          bargap=0.1,
                          font=dict(size=24))
        return fig

    def graph_insertions(self, chr_col_choice, color_by_choice, all_clrs):
        df_to_show = self.lib_map[self.lib_map[self.chr_col] == chr_col_choice].sort_values(self.insertion_site_col)
        fig = px.scatter(df_to_show, x=self.insertion_site_col, y=self.abundance_col, color=color_by_choice, log_y=True,
                         height=600, template='plotly_white',
                         color_discrete_sequence=all_clrs, hover_data=self.attributes,
                         labels={self.insertion_site_col: 'Position, bp',
                                 self.abundance_col: 'Read Counts'})
        fig.update_traces(marker=dict(size=8),
                          selector=dict(mode='markers'), opacity=0.7)
        fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'},
                          autosize=True,
                          font=dict(size=24))
        return fig

    def get_stats(self, name_col):
        table1 = (self.lib_map.groupby('library')
                  .agg({self.barcode_col: ['nunique'], self.distance_col: [lambda x: sum(x != 0)]})
                  .reset_index())
        table1.columns = ["Library", '# of insertions', '# of insertions outside of CDS']
        table1 = table1.set_index('Library')
        table1[f'# of {name_col}s with insertion'] = (self.lib_map[self.lib_map[self.distance_col] == 0]
                                                      .groupby('library')[name_col].nunique())
        table2 = (self.lib_map.groupby(['library', name_col])[self.barcode_col].count()
                  .reset_index().groupby('library')
                  .agg({self.barcode_col: ['median', 'max']})
                  .reset_index())
        table2.columns = ['Library', 'Median insertions per gene', 'Max insertions per gene']
        table2['Median insertions per gene'] = table2["Median insertions per gene"].astype(int)
        table2 = table2.set_index('Library')
        self.stats = table1.merge(table2, left_index=True, right_index=True)


class CountDataSet:
    def __init__(self, count_file, sample_data_file, config_file: str = 'scripts/config.yaml'):
        self.count_data = pd.read_csv(count_file)
        self.sample_data = pd.read_csv(sample_data_file).fillna('N/A')
        with open(config_file, 'r') as cf:
            config = yaml.load(cf, Loader=yaml.SafeLoader)['eda']
        # Load column naming schema
        self.col_name_config = config['fixed_column_names']
        self.fixed_column_names = list(self.col_name_config.values())
        self.barcode_col = self.col_name_config['barcode_col']
        self.gene_name_col = self.col_name_config['gene_name_col']
        self.sample_id_col = self.col_name_config['sample_id_col']
        self.valid = self._validate()

    def _validate(self):
        """
        First column of sample_data should be sampleIDs
        """
        st.write(f"Using {self.sample_data.columns[0]} to identify samples")
        self.sample_data = self.sample_data.rename({self.sample_data.columns[0]: self.sample_id_col}, axis=1)
        st.write(f"Using {self.count_data.columns[0]} to idnetify barcodes")
        st.write(f"Using {self.count_data.columns[1]} to idnetify genes")
        self.count_data = (self.count_data.rename({self.count_data.columns[0]: self.barcode_col,
                                                   self.count_data.columns[1]: self.gene_name_col}, axis=1)
                           .dropna(subset=[self.gene_name_col])
                           .drop_duplicates())
        samples_found = list(set(self.sample_data[self.sample_id_col].unique()).intersection(self.count_data.columns))
        if not samples_found or self.barcode_col in samples_found or self.gene_name_col in samples_found:
            st.error("No common samples found between sample data file and count table")
            return False
        self.sample_data = self.sample_data[self.sample_data[self.sample_id_col].isin(samples_found)]
        self.count_data = self.count_data[[self.barcode_col, self.gene_name_col] + samples_found]
        if self.count_data.empty:
            return False
        return True

    def normalize_counts(self):
        self.count_data = self.count_data.set_index([self.barcode_col, self.gene_name_col])
        self.count_data = self.count_data.loc[:, self.count_data.sum() > 0]
        self.count_data = np.log2((self.count_data / self.count_data.sum()) * 1000000 + 0.5).reset_index()

    def get_principal_components(self, numPCs, numGenes, chooseBy):
        """
        :param numPCs:
        :param numGenes:
        :param chooseBy:
        :return:
        """
        pcaDf = self.count_data.set_index(self.barcode_col).copy()
        pcaDf = pcaDf.drop(self.gene_name_col, axis=1)
        pcaSd = self.sample_data.set_index(self.sample_id_col).apply(lambda x: x.astype('category'))
        if numGenes:
            # calculate var for each, pick numGenes top var across samples -> df
            if chooseBy == 'variance':
                genes = pcaDf.var(axis=1).sort_values(ascending=False).head(int(numGenes)).index
                pcaDf = pcaDf.loc[genes].T
            else:
                pass
                # todo implement log2fc selection
        else:
            pcaDf = pcaDf.T
        pca = PCA(n_components=numPCs)
        principalComponents = pca.fit_transform(pcaDf)
        pcs = [f'PC{i}' for i in range(1, numPCs + 1)]
        pcDf = (pd.DataFrame(data=principalComponents, columns=pcs)
                .set_index(pcaDf.index))
        pcVar = {pcs[i]: round(pca.explained_variance_ratio_[i] * 100, 2) for i in range(0, numPCs)}
        pcDf = pcDf.merge(pcaSd, how="left", left_index=True, right_index=True)
        return pcDf, pcVar

    def pca_figure(self, pcDf, pcX, pcY, pcVarHi, pcVar, pcSym, expVars, colorSeq, w=None, h=None):
        h = h if h else 400
        w = w if w else 800
        fig = px.scatter(pcDf, x=pcX, y=pcY, color=pcVarHi, symbol=pcSym,
                         labels={pcX: f'{pcX}, {pcVar[pcX]} % Variance',
                                 pcY: f'{pcY}, {pcVar[pcY]} % Variance'},
                         color_discrete_sequence=colorSeq,
                         template='plotly_white',
                         height=h, width=w, hover_data=expVars, hover_name=pcDf.index)
        fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'},
                          autosize=True,
                          font=dict(size=16))
        fig.update_traces(marker=dict(size=20,
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

    def barcode_abundance_plot(self, geneDf, groupBy, colorBy, colorSeq, box=True):
        if box:
            fig = px.box(geneDf, x=groupBy, y='log2CPM', color=colorBy,
                         hover_data=geneDf.columns, points='all',
                         color_discrete_sequence=colorSeq, )
        else:
            fig = px.violin(geneDf, x=groupBy, y='log2CPM', color=colorBy,
                            hover_data=geneDf.columns, points='all',
                            color_discrete_sequence=colorSeq, )
        fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                          font=dict(size=16))
        fig.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='LightGrey')
        return fig


class ResultDataSet:
    def __init__(self, result_files=(), config_file="scripts/config.yaml",
                 gene_id='locus_tag'):
        self.result_files: str = result_files
        self.gene_id: str = gene_id
        self.results_df = pd.DataFrame()
        self.subset_df = pd.DataFrame()

        with open(config_file, 'r') as cf:
            config = yaml.load(cf, Loader=yaml.SafeLoader)['results']
        col_name_config = config['fixed_column_names']
        self.lfc_col = col_name_config['lfc_col']
        self.fdr_col = col_name_config['fdr_col']
        self.fdr_col2 = col_name_config['fdr_col2']
        self.contrast_col = col_name_config['contrast_col']
        self.library_col = col_name_config['library_col']
        self.string_df = pd.DataFrame()
        self.kegg_df = pd.DataFrame()
        self.colors, self.alphabetClrs, self.sushi_colors, self.all_clrs = define_color_scheme()

    def load_results(self):
        results_df_list = []
        for uploaded_result in self.result_files:
            st.write(f"Processing {uploaded_result.name}")
            df = pd.read_csv(uploaded_result)
            if 'library' not in df.columns:
                library_name = st.text_input("Add experiment name", value=uploaded_result.name.split("_rra")[0])
                df['library'] = library_name
            if self.gene_id not in df.columns:
                st.warning(f""" No {self.gene_id} column found. Using {df.columns[0]} as gene names to display""")
                self.gene_id = df.columns[0]
            results_df_list.append(df)
        try:
            fdf = pd.concat(results_df_list)
            fdf['fdr'] = np.where(fdf[self.lfc_col] < 0, fdf[self.fdr_col], fdf[self.fdr_col2])
            fdf['-log10FDR'] = -1 * np.log10(fdf['fdr'])
            fdf = fdf.fillna({self.gene_id: 'N/A'})
            self.results_df = fdf
        except ValueError:
            st.error('No result files loaded')

    def validate_results_df(self):
        results_schema = pa.DataFrameSchema({
            self.lfc_col: pa.Column(float, coerce=True),
            self.fdr_col: pa.Column(float, coerce=True),
            self.fdr_col2: pa.Column(float, coerce=True),
            self.contrast_col: pa.Column(str, coerce=True),
            self.library_col: pa.Column(str, coerce=True),
            'fdr': pa.Column(float),
            '-log10FDR': pa.Column(float)})

        try:
            self.results_df = results_schema.validate(self.results_df)
        except SchemaError as err:
            st.error(f"""Schema Error: {err.args[0]}""")
            self.results_df = pd.DataFrame()

    def identify_hits(self, library_to_show, lfc_low, lfc_hi, fdr_th):
        if not lfc_hi:
            self.results_df['hit'] = ((abs(self.results_df[self.lfc_col]) > lfc_low)
                                      & (self.results_df['fdr'] < fdr_th))
        else:
            self.results_df['hit'] = ((self.results_df[self.lfc_col] > lfc_low) &
                                      (self.results_df[self.lfc_col] < lfc_hi) &
                                      (self.results_df['fdr'] < fdr_th))

        if library_to_show != 'All':
            self.results_df = self.results_df[self.results_df[self.library_col] == library_to_show]
            self.results_df['LFC_median'] = self.results_df['LFC']
            self.results_df['library_nunique'] = 1
            self.results_df['hit_sum'] = self.results_df['hit']
        else:
            df_grouped = (self.results_df.groupby([self.gene_id, self.contrast_col])
                          .agg({self.lfc_col: ['median'], self.library_col: ['nunique'],
                                'hit': ['sum']})
                          .reset_index())
            df_grouped.columns = [self.gene_id, self.contrast_col, 'LFC_median', 'library_nunique', 'hit_sum']

            self.results_df = self.results_df.merge(df_grouped, on=[self.gene_id, self.contrast_col], how='left')

    def graph_by_rank(self, contrast=(), kegg=False):
        rank_df = self.kegg_df if kegg else self.results_df
        if contrast:
            rank_df = rank_df[rank_df[self.contrast_col].isin(contrast)]
        rank_df = (rank_df[[self.gene_id, self.contrast_col, 'LFC_median', 'hit', 'fdr']].drop_duplicates()
                   .sort_values('LFC_median')
                   .reset_index()
                   .reset_index()
                   .rename({'level_0': 'ranking'}, axis=1)
                   .sort_values('ranking'))

        hover_dict = {self.gene_id:True, self.contrast_col: True, 'LFC_median': True, 'fdr': True, 'ranking': False}

        fig = px.scatter(rank_df, x='ranking', y='LFC_median', color='hit', symbol=self.contrast_col,
                         height=500,
                         color_discrete_map={
                             True: self.colors['darko'],
                             False: self.colors['grey']},
                         hover_name=self.gene_id,
                         hover_data=hover_dict,
                         labels={"ranking": '', 'LFC_median': 'LFC'}
                         )
        fig.add_hline(y=0, line_width=2, line_dash="dash", line_color="grey")
        fig.update_xaxes(showticklabels=False)
        fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                          font=dict(size=18))
        fig.update_traces(marker=dict(size=14, opacity=0.8),
                          selector=dict(mode='markers'))
        return fig

    def graph_heatmap(self, genes):
        heat_df = (self.results_df[self.results_df[self.gene_id].isin(genes)][[self.gene_id, self.contrast_col, 'LFC_median']]
                   .drop_duplicates()
                   .pivot(index=self.gene_id, columns=self.contrast_col, values='LFC_median'))
        heat_df.index.name = 'Gene'
        fig = px.imshow(heat_df, color_continuous_scale=px.colors.diverging.Geyser,
                        color_continuous_midpoint=0,
                        width=1000, height=900)
        fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                          font=dict(size=10))
        return fig


    def display_pathway_heatmap(self, pathway_gene_names, kegg_id):
        if kegg_id not in self.results_df.columns:
            st.error(f"{kegg_id} not found in the results table")
        else:
            heat_df = self.results_df[self.results_df[kegg_id].isin(pathway_gene_names)]
            absent = pd.DataFrame(
                pd.Series(list(set(pathway_gene_names) - set(heat_df[kegg_id].unique())), name=kegg_id))
            heat_df = heat_df[[self.gene_id, kegg_id, 'LFC_median', self.contrast_col]].drop_duplicates()
            heat_df = pd.concat([heat_df, absent], axis=0)
            s = heat_df[self.gene_id].fillna('-').values
            p = heat_df[kegg_id].values
            heat_df['Gene'] = [
                f"{kegg_id}: {self.gene_id}" if self.gene_id != '-' and self.gene_id != kegg_id else f"{kegg_id}" for
                kegg_id, self.gene_id in zip(p, s)]
            heat_df = heat_df.pivot(index='Gene', columns=self.contrast_col, values='LFC_median')
            fig = px.imshow(heat_df, color_continuous_scale=px.colors.diverging.Geyser,
                            color_continuous_midpoint=0,
                            width=1000, height=900)
            fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'}, autosize=True,
                              font=dict(size=10))
            return fig
    # def subset_results(self, contrast_to_show):
    #     #df = self.results_df[~fdf[gene_id].str.contains(':')].copy() # todo come up with cleverer way to filter these
    #     self.subset_df = self.results_df[(self.results_df['contrast'].isin(contrast_to_show))].copy()


class KeggMapsDataset:

    def __init__(self, kegg_id, organism, results_df, gene_id=''):
        self.organism = organism
        self.kegg_id = kegg_id
        self.gene_id = gene_id if gene_id else kegg_id
        self.results_df = results_df
        self.gene_to_pathway = {}


    def validate_df(self):
        # kegg_id in results_df columns
        # other columns?
        # library
        # only 1 contrast
        pass

    @st.cache
    def get_org_kegg_pathways(self):
        result = pd.read_table(io.StringIO(kegg_list("pathway", self.organism).read()), header=None)
        result.columns = [f'KEGG_Pathway', 'Pathway_Description']
        result[f'KEGG_Pathway'] = result[f'KEGG_Pathway'].str.split(":").str.get(1)
        result['Pathway_Description'] = result['Pathway_Description'].str.split(" - ").str.get(0)
        result['KEGG_Display'] = result[f'KEGG_Pathway'] + ":" + result['Pathway_Description']
        path_map = result.set_index('KEGG_Display').to_dict()
        return path_map['KEGG_Pathway']


    def get_gene_to_pathway_dict(self):
        """
        Take the results df and convert to dictionary with color and names for each gene to display
        """
        if self.gene_id != self.kegg_id:
            self.results_df[self.gene_id] = self.results_df[self.gene_id].fillna(self.kegg_id)

        if self.results_df.library.nunique() > 1:
            self.results_df['hit'] = (self.results_df.groupby([self.kegg_id]).hit.sum()
                                      .reset_index())['hit']
        self.results_df['hitStar'] = self.results_df['hit'].apply(lambda x: '*' if x > 0 else '')
        self.results_df['NameForMap'] = self.results_df[self.gene_id] + self.results_df['hitStar'] + \
                                       "(" + self.results_df['LFC_median'].round(2).astype(str) + ")"
        self.results_df['NameForMapNum'] = self.results_df[self.gene_id].apply(self.parse_number_out) + self.results_df['hitStar'] + \
                                        "(" + self.results_df['LFC_median'].round(2).astype(str) + ")"
        data_short = (self.results_df[{self.gene_id, self.kegg_id, 'NameForMap', 'NameForMapNum', 'LFC_median'}]
                      .drop_duplicates()
                      .dropna())
        norm = matplotlib.colors.Normalize(vmin=-6, vmax=6, clip=True) # todo remove hard min and max values
        mapper = matplotlib.cm.ScalarMappable(norm=norm, cmap=sns.diverging_palette(220, 20, as_cmap=True)) # Try matplotlib RdYlBu
        data_short['hex'] = data_short.LFC_median.apply(mapper.to_rgba).apply(matplotlib.colors.to_hex)
        data_short['hex'] = data_short.hex.str.replace('#000000', "#7c7d83")
        self.gene_to_pathway = data_short.set_index(self.kegg_id).to_dict()


    def parse_number_out(self, gene_name: Union[str, None]) -> Union[str, None]:
        """
        :param gene_name: KEGG gene name in the following format: organism:gene_name,
                examples: eco:b3451 sey:SL1344_1569 ece:Z876
        :param numbers_only: Whether to return full gene name or numbers only

        :return: Parsed gene name
        :rtype: str

        """

            # match = re.findall(r"(?<=[a-zA-Z_])\d+", string)
            # positive lookbehind (?<=)
            # searches from behind --> if match takes the part
            # after the match
            # \D means "not a digit"
            # \d matches a digit (equivalent to [0-9])
            # and matches the previous token between one and unlimited times,
            # as many times as possible, giving back as needed (greedy)
        if not gene_name:
            return
        match = re.findall(r"(?<=[\D])\d+", gene_name)
        if match:
            # takes the last match, because a list is provided
            number = match[-1]
            return number if number.isnumeric() else gene_name
        return gene_name

    def displayPDF(self, file):
        # Opening file from file path
        with open(file, "rb") as f:
            base64_pdf = base64.b64encode(f.read()).decode('utf-8')
        # Embedding PDF in HTML
        pdf_display = F'<iframe src="data:application/pdf;base64,{base64_pdf}" width="700" height="1000" type="application/pdf"></iframe>'
        # Displaying File
        st.markdown(pdf_display, unsafe_allow_html=True)

    def display_kegg_map(self, pathway_name, title, numeric=False):
        """
        ko_dict: {hex: {ko: color}, NameForMap: {ko: name}}

        """
        pathway_kgml = KGML_parser.read(kegg_get(pathway_name, "kgml"))
        pathway_gene_names = [gene.name.split() for gene in pathway_kgml.genes]
        pathway_gene_names = set([gene.split(":")[1] for sublist in pathway_gene_names for gene in sublist])
        canvas = KGMLCanvas(pathway_kgml, import_imagemap=True)
        not_found = []
        for element in pathway_kgml.genes:
            color = None
            name = None
            node_kos = [e.split(":")[1] for e in element.name.split()]
            for ko in node_kos:
                color = self.gene_to_pathway['hex'].get(ko, color)
                name = self.gene_to_pathway['NameForMapNum'].get(ko, name) if numeric else self.gene_to_pathway['NameForMap'].get(ko, name)
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
        if k1.button(f'Display {pathway_name} map'):
            self.displayPDF(fname)
        with open(fname, "rb") as f:
            k2.download_button(
                f"Download {pathway_name} map",
                data=f,
                file_name=fname,
            )
        return pathway_gene_names

