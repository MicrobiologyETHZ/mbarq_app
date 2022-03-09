import pandas as pd
from pathlib import Path

def proces_gmt(gmt_file, path_name=1):
    genes = []
    paths = []
    with open(gmt_file, 'r', errors='ignore') as fh:
        for line in fh.readlines():
            for gene in line.split('\t')[path_name:]:
                genes.append(gene.strip("\n"))
                paths.append(":".join(line.split('\t')[0:path_name]))
    return pd.DataFrame([genes, paths], index=['Name', 'KEGG_Pathway']).T


def merge_with_results(results_files, gmt_df, left_on):
    for results_file in results_files:
        resDf = pd.read_csv(results_file).rename({left_on: 'Name'}, axis=1)
        resDf = resDf.merge(gmt_df, how='left', on='Name')
        resDf.to_csv(Path(results_file).with_suffix('.kegg.csv'), index=False)



if __name__ == "__main__":
    gmt_file = "/Users/ansintsova/git_repos/mbarq_app/data/SL1344_test/04-03-2022-SL1344-KEGG-API.gmt"
    res_files = ["/Users/ansintsova/git_repos/mbarq_app/data/SL1344_test/library_10_1-unfiltered-results.csv"]
    gmt_df = proces_gmt(gmt_file, path_name=2)

    merge_with_results(res_files, gmt_df, left_on='id')
