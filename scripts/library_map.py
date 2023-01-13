from pathlib import Path
from typing import List

import pandas as pd
import pandera as pa
import plotly.express as px
import streamlit as st
import yaml
from pandera.errors import SchemaError

from scripts.graphs import define_color_scheme

with open('scripts/config.yaml', 'r') as cf:
    config = yaml.load(cf, Loader=yaml.SafeLoader)['library_map']

# Load column naming schema
col_name_config = config['fixed_column_names']
FIXED_COLUMN_NAMES = list(col_name_config.values())
CHR_COL = col_name_config['chr_col']
INSERTION_SITE_COL = col_name_config['insertion_site_col']
ABUNDANCE_COL = col_name_config['abundance_col']
BARCODE_COL = col_name_config['barcode_col']
DISTANCE_COL = col_name_config['distance_col']


class LibraryMap:

    def __init__(self, map_files: List):
        self.map_files = map_files
        self.lib_map = pd.DataFrame()
        self.attributes = []
        self.color_by_cols = ('in CDS', 'library')
        self.stats = pd.DataFrame

    def _read_in_library_map(self, uploaded_map):
        return pd.read_csv(uploaded_map), uploaded_map.name

    def _load_library_map(self, df, df_name):
        df['in CDS'] = df.distance_to_feature == 0
        if 'library' not in df.columns:
            library_name = st.text_input("Change library name?", value=df_name)
            df['library'] = library_name
        if 'number_of_reads' in df.columns:  # todo: this is for backwards compatibility remove this in the future
            df = df.rename({'number_of_reads': 'abundance_in_mapping_library'}, axis=1)
        missing_cols = [c for c in FIXED_COLUMN_NAMES if c not in df.columns]
        if len(missing_cols) > 0:
            st.markdown(
                f"""⚠️ The following columns are missing from the map files: {', '.join(missing_cols)}. 
                Please rename the columns/rerun mBARq and try again. Skipping {df_name.name} ⚠️""")
            return pd.DataFrame(), []
        attrs = [c for c in df.columns if c not in FIXED_COLUMN_NAMES]
        return df, attrs

    def process_library_maps(self):
        map_dfs = []
        attr_names = []
        for uploaded_map in self.map_files:
            st.write(f"Processing {uploaded_map.name}")
            df, attrs = self._load_library_map(*self._read_in_library_map(uploaded_map))
            map_dfs.append(df)
            attr_names.extend(attrs)
        self.lib_map = pd.concat(map_dfs)
        self.attributes = set(attr_names)

    def validate_lib_map(self):
        lib_schema = pa.DataFrameSchema({
            CHR_COL: pa.Column(str, coerce=True),
            INSERTION_SITE_COL: pa.Column(int, pa.Check(lambda x: x >= 0)),
            BARCODE_COL: pa.Column(str, coerce=True),
            ABUNDANCE_COL: pa.Column(int, pa.Check(lambda x: x >= 0)),
            DISTANCE_COL: pa.Column(float, nullable=True)
        }
        )

        try:
            self.lib_map = lib_schema.validate(self.lib_map)

        except SchemaError as err:
            st.error(f"""Schema Error: {err.args[0]}""")
            self.lib_map = pd.DataFrame()

    def graph_coverage_hist(self, chr_col_choice, num_bins, hist_col):
        df_to_show = self.lib_map[self.lib_map[CHR_COL] == chr_col_choice].sort_values(INSERTION_SITE_COL)
        fig = px.histogram(df_to_show, x=INSERTION_SITE_COL, nbins=int(num_bins),
                           labels={INSERTION_SITE_COL: 'Position, bp'}, color_discrete_sequence=[hist_col])
        fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'},
                          autosize=True,
                          bargap=0.1,
                          font=dict(size=24))
        return fig

    def graph_insertions(self, chr_col_choice, color_by_choice, all_clrs):
        df_to_show = self.lib_map[self.lib_map[CHR_COL] == chr_col_choice].sort_values(INSERTION_SITE_COL)
        fig = px.scatter(df_to_show, x=INSERTION_SITE_COL, y=ABUNDANCE_COL, color=color_by_choice, log_y=True,
                         height=600, template='plotly_white',
                         color_discrete_sequence=all_clrs, hover_data=self.attributes,
                         labels={INSERTION_SITE_COL: 'Position, bp',
                                 ABUNDANCE_COL: 'Read Counts'})
        fig.update_traces(marker=dict(size=8, line=dict(width=2, color='DarkSlateGrey')),
                          selector=dict(mode='markers'), opacity=0.9)
        fig.update_layout({'paper_bgcolor': 'rgba(0,0,0,0)', 'plot_bgcolor': 'rgba(0,0,0,0)'},
                          autosize=True,
                          font=dict(size=24))
        return fig

    def get_stats(self, name_col):
        table1 = (self.lib_map.groupby('library')
                  .agg({BARCODE_COL: ['nunique'], DISTANCE_COL: [lambda x: sum(x != 0)]})
                  .reset_index())
        table1.columns = ["Library", '# of insertions', '# of insertions outside of CDS']
        table1 = table1.set_index('Library')
        table1[f'# of {name_col}s with insertion'] = (self.lib_map[self.lib_map[DISTANCE_COL] == 0]
                                                      .groupby('library')[name_col].nunique())
        table2 = (self.lib_map.groupby(['library', name_col])
                  .barcode.count().reset_index().groupby('library')
                  .agg({BARCODE_COL: ['median', 'max']})
                  .reset_index())
        table2.columns = ['Library', 'Median insertions per gene', 'Max insertions per gene']
        table2['Median insertions per gene'] = table2["Median insertions per gene"].astype(int)
        table2 = table2.set_index('Library')
        self.stats = table1.merge(table2, left_index=True, right_index=True)


@st.cache
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv(index=False).encode('utf-8')


def app():
    st.markdown(""" # Library Map """)

    with st.expander('How this works: '):
        st.markdown("""
    
        ### Visualize insertion position along the genome.
    
        - Takes in library map **CSV** file (`*.annotated.csv`) produced by `mbarq map`. Expects to find the following columns in the file: `barcode`, `abundance_in_mapping_library`, `insertion_site`, `chr`, `distance_to_feature`, `strand`. 
        - Can load more than one library file at the same time to compare! 
        - You can select which sequence (e.g. chromosome or plasmids) to display, and color the insertions by distance to feature, multimapping or library.
        - You can click on the figure legend to only show a specific subset of data (i.e. if looking at multiple libraries, double clicking on the specific library name will show data for that library only)
    
        """)

    with st.container():
        # Get the data
        st.subheader('Load your own data or browse the example data set')
        data_type = st.radio('Choose dataset to show', ['Look at an example', 'Load my data'], index=1, key='lib')
        if data_type == 'Load my data':
            map_files = st.file_uploader('Upload library map file', accept_multiple_files=True)
        else:
            map_files = [Path('examples/example_library_mapping_file.annotated.csv')]
            st.subheader('Example mapping file')  # todo update example mapping file not to include multimap column
            example_df = pd.read_csv(map_files[0])
            st.write(example_df.sample(5))
            st.download_button(
                label="Download example data as CSV",
                data=convert_df(example_df),
                file_name='example_library_mapping_file.csv',
                mime='text/csv',
            )

    if map_files:
        # Process and validate library files
        lm = LibraryMap(map_files)
        lm.process_library_maps()
        lm.validate_lib_map()
        if lm.lib_map.empty:
            st.error(f"""⚠️ Something went wrong when processing library map files. 
                            Please check the file formats and try again ⚠️""")
            st.stop()
        if st.checkbox("Show sample of the Library Map?"):
            st.write(lm.lib_map.sample(10))
        # Graph coverage map or individual insertion abundance
        with st.container():
            # Define colors
            colors, alphabetClrs, all_clrs = define_color_scheme()
            graph_type = st.radio("Choose graph", ['Coverage Histogram', 'Individual Insertions'])
            c1, c2, c3 = st.columns(3)
            chr_col_choice = c1.selectbox('Choose sequence to display', lm.lib_map[CHR_COL].unique())
            if graph_type == 'Individual Insertions':
                color_by_choice = c2.selectbox('Color by', lm.color_by_cols)
                fig = lm.graph_insertions(chr_col_choice, color_by_choice, all_clrs)
            else:
                num_bins = c2.number_input('Number of bins', value=100)
                hist_col = c3.text_input('Color hex', value=colors['teal'])
                fig = lm.graph_coverage_hist(chr_col_choice, num_bins, hist_col)
            st.plotly_chart(fig, use_container_width=True)

        # Generate summary stats for the libraries
        with st.container():
            name_col = st.selectbox("Choose attribute column", lm.attributes)
            lm.get_stats(name_col)
            st.markdown("#### Insertion Summary")
            st.table(lm.stats)
