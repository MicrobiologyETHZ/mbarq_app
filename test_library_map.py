import pandera as pa
import yaml
from scripts.library_map import LibraryMap

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


def test_get_stats():
    input_schema = pa.DataFrameSchema({
            CHR_COL: pa.Column(str, coerce=True),
            INSERTION_SITE_COL: pa.Column(int, pa.Check(lambda x: x >= 0)),
            BARCODE_COL: pa.Column(str, coerce=True),
            ABUNDANCE_COL: pa.Column(int, pa.Check(lambda x: x >= 0)),
            DISTANCE_COL: pa.Column(float, nullable=True),
            'in CDS': pa.Column(bool),
            'library': pa.Column(str, coerce=True)
        }
        )

    output_schema = pa.DataFrameSchema(columns={"Library": pa.Column(str),
                                       '# of insertions': pa.Column(int),
                                       '# of insertions outside of CDS': pa.Column(int),
                                       'Median insertions per gene': pa.Column(int),
                                       'Max insertions per gene': pa.Column(int)},
                                       checks=pa.Check(lambda df: df['Median insertions per gene'] <= df['Max insertions per gene']))
    pass
    # todo add other checks

