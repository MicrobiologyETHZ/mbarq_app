from pathlib import Path
from typing import List, Union
import pandas as pd
import pandera as pa
import plotly.express as px
import streamlit as st
import yaml
from pandera.errors import SchemaError
import numpy as np
from sklearn.decomposition import PCA
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas

import base64
import matplotlib
import seaborn as sns

import re

@st.cache_data
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv(index=False).encode('utf-8')


def define_color_scheme():
    alphabet_clrs = px.colors.qualitative.Dark24
    app_colors = {'grey': "#E2E2E2",
                  'red': '#C0504D',
                  'light_yellow': "#f7ba65",
                  'darko': "#bf4713",
                  'maroon': "#9c002f",
                  'brighto': "#d73d00",
                  'teal': "#008080",
                  'darkteal': "#004c4c",
                    'orange': '#F79646',
                    'medSea': '#4BACC6',
                    'black': '#000000',
                    'dgreen': '#00B04E',
                    'lgreen': '#92D050',
                    'dblue': '#366092',
                    'lblue': '#95B3D7'}
    all_clrs = ['#F79646', '#366092', '#8fce00', '#C0504D', app_colors['teal'], app_colors['maroon']] + alphabet_clrs
    return alphabet_clrs, app_colors, all_clrs


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
        opt_name_config = config['optional_column_names']
        self.fixed_column_names = list(col_name_config.values())
        self.optional_column_names = list(opt_name_config.values())
        self.chr_col = col_name_config['chr_col']
        self.insertion_site_col = col_name_config['insertion_site_col']
        self.abundance_col = col_name_config['abundance_col']
        self.barcode_col = col_name_config['barcode_col']
        self.distance_col = col_name_config['distance_col']

    def load_map(self, silent=False):
        map_dfs = []
        if self.map_files:
            for uploaded_map in self.map_files:
                if not silent:
                    st.write(f"_Processing {uploaded_map.name}_")
                df, df_name = pd.read_csv(uploaded_map), uploaded_map.name
                if 'library' not in df.columns:
                    library_name = df_name
                    if not silent:
                        library_name = st.text_input("Change library name?", value=library_name)
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
            self.lib_map['in CDS'] = self.lib_map[self.distance_col] == 0
        self.attributes = [c for c in self.lib_map.columns if c not in self.fixed_column_names
                           and c not in self.optional_column_names + ['library', 'in CDS']]

    def validate_lib_map(self):
        lib_schema = pa.DataFrameSchema({
            self.chr_col: pa.Column(str, coerce=True),
            self.insertion_site_col: pa.Column(int, pa.Check(lambda x: x >= 0)),
            self.barcode_col: pa.Column(str, coerce=True),
            self.abundance_col: pa.Column(int, pa.Check(lambda x: x >= 0)),
            self.distance_col: pa.Column(int, coerce=True, nullable=True),
            'in CDS': pa.Column(bool),
            'library': pa.Column(str, coerce=True)
        }
        )
        try:
            self.lib_map = lib_schema.validate(self.lib_map)

        except SchemaError as err:
            # st.error(f"""Schema Error: {err.args[0]}""")
            self.lib_map = pd.DataFrame()

    def graph_coverage_hist(self, chr_col_choice, num_bins, hist_col):
        df_to_show = self.lib_map[self.lib_map[self.chr_col] == chr_col_choice].sort_values(self.insertion_site_col)
        fig = px.histogram(df_to_show, x=self.insertion_site_col, nbins=int(num_bins),
                           labels={self.insertion_site_col: 'Position, bp'}, 
                           color_discrete_sequence=[hist_col])
        fig.update_layout(bargap=0.1)
        fig.update_xaxes(showline=True, linewidth=1, linecolor='black',
                         tickfont=dict(size=24, color='black'), title_font=dict(size=30, color='black'))
        fig.update_yaxes(showline=True, linewidth=1, linecolor='black',
                         tickfont=dict(size=24, color='black'), title_font=dict(size=30, color='black'))
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
        fig.update_xaxes(showline=True, linewidth=1, linecolor='black',
                         tickfont=dict(size=18, color='black'), title_font=dict(size=24, color='black'))
        fig.update_yaxes(showline=True, linewidth=1, linecolor='black',
                         tickfont=dict(size=18, color='black'), title_font=dict(size=24, color='black'))
        return fig

    def get_stats(self):
        table1 = (self.lib_map.groupby('library')
                  .agg({self.barcode_col: ['nunique'], self.distance_col: [lambda x: sum(x != 0)]})
                  .reset_index())
        table1.columns = ["Library", '# of insertions', '# of insertions outside of CDS']
        table1 = table1.set_index('Library')
        for att in self.attributes:
            table1[f'# of {att}s with insertion'] = (self.lib_map[self.lib_map[self.distance_col] == 0]
                                                     .groupby('library')[att].nunique())
        table2 = (self.lib_map.groupby(['library', self.attributes[0]])[self.barcode_col].count()
                  .reset_index().groupby('library')
                  .agg({self.barcode_col: ['median', 'max']})
                  .reset_index())
        table2.columns = ['Library', f'Median insertions per {self.attributes[0]}', f'Max insertions per {self.attributes[0]}']
        table2[f'Median insertions per {self.attributes[0]}'] = table2[f"Median insertions per {self.attributes[0]}"].astype(int)
        table2 = table2.set_index('Library')
        self.stats = table1.merge(table2, left_index=True, right_index=True)


class CountDataSet:
    def __init__(self, count_file, sample_data_file, config_file: str = 'scripts/config.yaml', silent=False):
        self.count_file = count_file
        self.sample_data_file = sample_data_file
        self.count_data = pd.read_csv(self.count_file)
        self.sample_data = pd.read_csv(self.sample_data_file).fillna('N/A')
        with open(config_file, 'r') as cf:
            config = yaml.load(cf, Loader=yaml.SafeLoader)['eda']
        # Load column naming schema
        self.col_name_config = config['fixed_column_names']
        self.fixed_column_names = list(self.col_name_config.values())
        self.barcode_col = self.col_name_config['barcode_col']
        self.gene_name_col = self.col_name_config['gene_name_col']
        self.sample_id_col = self.col_name_config['sample_id_col']
        self.valid = self._validate(silent=silent)
        self.norm_counts = pd.DataFrame()

    def _validate(self, silent=False):
        """
        First column of sample_data should be sampleIDs
        """
        if not silent:
            st.write(f"_Using {self.sample_data.columns[0]} to identify samples_")
        self.sample_data = self.sample_data.rename({self.sample_data.columns[0]: self.sample_id_col}, axis=1)
        if not silent:
            st.write(f"_Using {self.count_data.columns[0]} to identify barcodes_")
            st.write(f"_Using {self.count_data.columns[1]} to identify genes_")
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
        self.norm_counts = self.count_data.set_index([self.barcode_col, self.gene_name_col]).copy()
        self.norm_counts = self.norm_counts.loc[:, self.norm_counts.sum() > 0]
        self.norm_counts = np.log2((self.norm_counts / self.norm_counts.sum()) * 1000000 + 0.5).reset_index()

    def get_principal_components(self, numPCs, numGenes, chooseBy):
        """
        :param numPCs:
        :param numGenes:
        :param chooseBy:
        :return:
        """
        pcaDf = self.norm_counts.set_index(self.barcode_col).copy()
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
        pcDf = pcDf[~pcDf.isna().any(axis=1)]
        return pcDf, pcVar

    def pca_figure(self, pc_df, pc_x_axis, pc_y_axis, highlight_var, percent_variance,
                   symbol_var, experiment_vars, color_sequence, w=None, h=None,
                   font_size=24):
        h = h if h else 400
        w = w if w else 800
        font_size = max(font_size, 8)
        fig = px.scatter(pc_df, x=pc_x_axis, y=pc_y_axis, color=highlight_var, symbol=symbol_var,
                         labels={pc_x_axis: f'{pc_x_axis}, {percent_variance[pc_x_axis]} % Variance',
                                 pc_y_axis: f'{pc_y_axis}, {percent_variance[pc_y_axis]} % Variance'},
                         color_discrete_sequence=color_sequence,
                         template='plotly_white',
                         height=h, width=w, hover_data=experiment_vars, hover_name=pc_df.index)
        fig.update_traces(marker=dict(size=20,
                                      line=dict(width=2,
                                                color='DarkSlateGrey'), opacity=0.9),
                          selector=dict(mode='markers'))
        fig.update_xaxes(showline=True, linewidth=2, linecolor='black',
                         tickfont=dict(size=font_size-6, color='black'), title_font=dict(size=font_size, color='black'))
        fig.update_yaxes(showline=True, linewidth=2, linecolor='black',
                         tickfont=dict(size=font_size-6, color='black'), title_font=dict(size=font_size, color='black'))
        fig.update_layout(legend=dict(font=dict(size=font_size)), legend_title=dict(font=dict(size=font_size)))
        percent_variance_df = pd.DataFrame.from_dict(percent_variance, orient='index').reset_index()
        percent_variance_df.columns = ['PC', '% Variance']
        fig2 = px.line(percent_variance_df, x='PC', y='% Variance', markers=True,
                       labels={'PC': ''})
        fig2.update_traces(marker=dict(size=12,
                                       line=dict(width=2,
                                                 color='DarkSlateGrey')))

        pc_summary = pc_df.groupby(highlight_var)[[c for c in pc_df.columns if c.startswith('PC')]].median()
        fig3 = px.imshow(pc_summary)
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
                 gene_id='Name'):
        self.result_files: str = result_files
        self.gene_id: str = gene_id
        self.results_df = pd.DataFrame()
        self.subset_df = pd.DataFrame()
        self.hit_df = pd.DataFrame()
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
        self.alphabet_clrs, self.app_colors, self.all_clrs = define_color_scheme()

    def load_results(self, silent=False):
        results_df_list = []
        for uploaded_result in self.result_files:
            if not silent:
                st.write(f"_Processing {uploaded_result.name}_")
            df = pd.read_csv(uploaded_result)
            if 'library' not in df.columns:
                library_name = uploaded_result.name.split("_rra")[0]
                if not silent:
                    library_name = st.text_input("Add experiment name", value=library_name)
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
        except KeyError:
            # todo rethink validation and column generation
            st.error(f'Could not find one of the following columns: '
                     f'{", ".join([self.lfc_col, self.fdr_col, self.fdr_col2])}. Wrong file format?')
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


    def identify_hits(self, library_to_show, lfc_low, lfc_hi, fdr_th):
        self.hit_df = self.results_df.copy()
        if not lfc_hi:
            self.hit_df['hit'] = ((abs(self.hit_df[self.lfc_col]) > lfc_low)
                                      & (self.hit_df['fdr'] < fdr_th))
        else:
            self.hit_df['hit'] = ((self.hit_df[self.lfc_col] > lfc_low) &
                                      (self.hit_df[self.lfc_col] < lfc_hi) &
                                      (self.hit_df['fdr'] < fdr_th))
        if 'LFC_median' in self.hit_df.columns:
            self.hit_df = self.hit_df.drop('LFC_median', axis=1)
        if library_to_show != 'All':
            self.hit_df = self.hit_df[self.hit_df[self.library_col] == library_to_show]
            self.hit_df['LFC_median'] = self.hit_df['LFC']
            self.hit_df['library_nunique'] = 1
            self.hit_df['hit_sum'] = self.hit_df['hit']
        else:
            # How to define hits from multiple libraries
            # Median LFC > cutoff + hit in at least 1 library
            df_grouped = (self.hit_df.groupby([self.gene_id, self.contrast_col])
                          .agg({self.lfc_col: ['median'], self.library_col: ['nunique'],
                                'hit': ['sum']})
                          .reset_index())
            df_grouped.columns = [self.gene_id, self.contrast_col, 'LFC_median', 'library_nunique', 'hit_sum']

            if not lfc_hi:
                df_grouped['hit_all'] = ((abs(df_grouped['LFC_median']) > lfc_low)
                                        & (df_grouped['hit_sum'] > 0))
            else:
                df_grouped['hit_all'] = ((df_grouped['LFC_median'] > lfc_low) &
                                        (df_grouped['LFC_median'] < lfc_hi) &
                                      (df_grouped['hit_sum'] > 0))
            self.hit_df = self.hit_df.merge(df_grouped, on=[self.gene_id, self.contrast_col], how='left')

    def graph_by_rank(self,  contrast=(), kegg=False, multiple_libs=False):
        rank_df = self.kegg_df if kegg else self.hit_df
        if contrast:
            rank_df = rank_df[rank_df[self.contrast_col].isin(contrast)]
        hit_col = 'hit_all' if multiple_libs else 'hit'
        rank_df = (rank_df[[self.gene_id, self.contrast_col, 'LFC_median', hit_col, 'fdr']].drop_duplicates()
                   .sort_values('LFC_median')
                   .reset_index()
                   .reset_index()
                   .rename({'level_0': 'ranking'}, axis=1)
                   .sort_values('ranking'))

        hover_dict = {self.gene_id: True, self.contrast_col: True, 'LFC_median': True, 'fdr': True, 'ranking': False}

        fig = px.scatter(rank_df, x='ranking', y='LFC_median', color=hit_col, symbol=self.contrast_col,
                         height=500,
                         color_discrete_map={
                             True: self.app_colors['darko'],
                             False: self.app_colors['grey']},
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

    def graph_heatmap(self, genes, font_size=24):
        heat_df = (self.hit_df[self.hit_df[self.gene_id].isin(genes)][[self.gene_id, self.contrast_col, 'LFC_median']]
                   .drop_duplicates()
                   .pivot(index=self.gene_id, columns=self.contrast_col, values='LFC_median'))
        heat_df.index.name = 'Gene'
        fig = px.imshow(heat_df, color_continuous_scale=px.colors.diverging.Geyser,
                        color_continuous_midpoint=0,
                        width=1000, height=900)

        fig.update_xaxes(showline=True, linewidth=2, linecolor='black',
                         tickfont=dict(size=font_size - 6, color='black'),
                         title_font=dict(size=font_size, color='black'), tickangle=90)
        fig.update_yaxes(showline=True, linewidth=2, linecolor='black',
                         tickfont=dict(size=font_size - 6, color='black'),
                         title_font=dict(size=font_size, color='black'))
        fig.update_layout(legend=dict(font=dict(size=font_size - 2, color='black')),
                          coloraxis_colorbar=dict(title=dict(text='LFC', font=dict(size=font_size, color='black')),
                                                  tickfont=dict(size=font_size - 2, color='black')))

        return fig

    def display_pathway_heatmap(self, pathway_gene_names, kegg_id):

        if kegg_id not in self.results_df.columns:
            st.error(f"{kegg_id} not found in the results table")
        else:
            heat_df = self.hit_df[self.hit_df[kegg_id].isin(pathway_gene_names)]
            absent = pd.DataFrame(
                pd.Series(list(set(pathway_gene_names) - set(heat_df[kegg_id].unique())), name=kegg_id))

            heat_df = heat_df[list({self.gene_id, kegg_id, 'LFC_median', self.contrast_col})].drop_duplicates()
            heat_df = pd.concat([heat_df, absent], axis=0)
            gene_names = heat_df[self.gene_id].fillna('-').values
            kegg_tags = heat_df[kegg_id].values

            heat_df['Gene'] = [
                f"{tag}: {name}" if name != '-' and name != tag else f"{tag}" for
                tag, name in zip(kegg_tags, gene_names)]


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

    @st.cache_data
    def get_org_kegg_pathways(_self, organism):
        result = pd.read_table(io.StringIO(kegg_list("pathway", organism).read()), header=None)
        result.columns = [f'KEGG_Pathway', 'Pathway_Description']
        #result[f'KEGG_Pathway'] = result[f'KEGG_Pathway'].str.split(":").str.get(1)
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
        data_short = (self.results_df[list({self.gene_id, self.kegg_id, 'NameForMap', 'NameForMapNum', 'LFC_median'})]
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
        with open(fname, "rb") as f:
            k1.download_button(
                f"Download {pathway_name} map",
                data=f,
                file_name=fname,
            )
        return pathway_gene_names

