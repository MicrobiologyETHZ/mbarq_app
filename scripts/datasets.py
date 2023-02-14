from pathlib import Path
from typing import List
import pandas as pd
import pandera as pa
import plotly.express as px
import streamlit as st
import yaml
from pandera.errors import SchemaError
from scripts.graphs import define_color_scheme
import numpy as np

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
        fig.update_traces(marker=dict(size=8, line=dict(width=2, color='DarkSlateGrey')),
                          selector=dict(mode='markers'), opacity=0.9)
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



class ResultDataSet:
    def __init__(self, result_files=(), config_file = "scripts/config.yaml",
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
       # self.string_df = pd.DataFrame()
       # self.kegg_df = pd.DataFrame()

    def load_results(self):
        results_df_list = []
        for uploaded_result in self.result_files:
            st.write(f"Processing {uploaded_result.name}")
            df = pd.read_csv(uploaded_result)
            if 'library' not in df.columns:
                library_name = st.text_input("Add library name?", value=uploaded_result.name)
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


    # def subset_results(self, contrast_to_show):
    #     #df = self.results_df[~fdf[gene_id].str.contains(':')].copy() # todo come up with cleverer way to filter these
    #     self.subset_df = self.results_df[(self.results_df['contrast'].isin(contrast_to_show))].copy()