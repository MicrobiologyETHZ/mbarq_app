import streamlit as st
from scripts.datasets import ResultDataSet
from pathlib import Path
import requests
from time import sleep

def app():
    st.markdown(""" # Analyze fitness results with STRING-db """)

    with st.expander('How this works: '):
        an_url = "https://mbarq.readthedocs.io/en/latest/analysis.html"
        st.markdown(f"""#### Fitness data: 
        - a `csv` file produced by `mbarq analyze` command. To learn more about how to use `mbarq analyze`, please read [here]({an_url}).
        - First column must be a gene identifier (for example, locus tag). 
        - Must also include `LFC` and `contrast` columns, where `LFC` is log2 fold change in gene abundance for a specific treatment compared to control, and `contrast` specifies the treatment.  
        - You can define hits by setting LFC and FDR cutoffs, and submit gene identifiers of the hits to STRING search.
        - Make sure that the gene identifier you used for analysis is recognized by STRING. 
        - If you load the library map on the **Data Upload** page, you would be able to choose which identifier to use for STRING (for example, by default library map will have Name, locus tag and ID). 
        """)

    with st.container():
        # Get the data
        if 'results_ds' in st.session_state.keys():
            rds = st.session_state['results_ds']
        else:
            st.info('Browse example results file or upload your data in **⬆️ Data Upload**')
            result_files = [Path("examples/string_example_rra_results.csv")]
            gene_id = 'Name'
            rds = ResultDataSet(result_files=result_files, gene_id=gene_id)
            rds.load_results()
            rds.validate_results_df()

        if st.checkbox('Show sample of the dataset'):
            try:
                st.write(rds.results_df.sample(5))
            except ValueError:
                st.write('Result table is empty')

    if not rds.results_df.empty:
        if 'attributes' in st.session_state.keys():
            gene_identifier = st.selectbox('Choose gene identifier',
                                           st.session_state['attributes'])
        else:
            gene_identifier = rds.gene_id
        st.info(
            f"❗Make sure STRING recognizes unique gene identifier (you've entered `{gene_identifier}`) for the taxon you specify")
        species = st.number_input("NCBI species taxid", value=99287, help='Salmonella Typhimurium: 99287')
        contrasts = rds.results_df[rds.contrast_col].sort_values().unique()
        libraries = rds.results_df[rds.library_col].sort_values().unique()
        if len(libraries) > 1:
            libraries = ['All'] + list(libraries)
            library_to_show = st.selectbox('Select experiment to show', libraries)
        else:
            library_to_show = libraries[0]
        # SUBSET DATAFRAME TO SPECIFIC CONTRAST AND LIBRARY
        contrast_col, lfc_col1, lfc_col2, fdr_col = st.columns(4)
        contrast_to_show = contrast_col.selectbox('Select a contrast', contrasts)
        fdr_th = fdr_col.number_input('adjusted p value cutoff', value=0.05)
        type_lfc_th = lfc_col1.radio('Absolute LFC cutoff or define range', ['Absolute', 'Range'])
        if type_lfc_th == 'Absolute':
            lfc_low = lfc_col2.number_input('LFC cutoff (absolute)', min_value=0.0, step=0.5, value=1.0)
            lfc_hi = None
        else:
            lfc_low = lfc_col2.number_input('Min LFC',  step=0.5, value=-5.0)
            lfc_hi = lfc_col2.number_input('Max LFC',  step=0.5, value=-1.0)

        rds.identify_hits(library_to_show,  lfc_low, lfc_hi, fdr_th)
        string_df = rds.hit_df[rds.hit_df[rds.contrast_col] == contrast_to_show]
        st.markdown(
            f"Analyze hits for ``{contrast_to_show}`` contrast for ``{library_to_show}`` experiment(s)")
        up = st.radio('Up or Down?', ('Upregulated Only', 'Downregulated Only', 'Both'), key='string')
        if up == 'Upregulated Only':
            string_df = string_df[(string_df['LFC'] > 0) & (string_df['hit'] == True)]
        elif up == 'Downregulated Only':
            string_df = string_df[(string_df['LFC'] < 0) & (string_df['hit'] == True)]
        else:
            string_df = string_df[string_df['hit'] == True]

        if lfc_hi:
            st.markdown(
                f"There are {string_df[gene_identifier].nunique()} hits with FDR < {round(fdr_th, 2)} and within LFC range from {lfc_low} to {lfc_hi}")
        else:
            st.markdown(
                f"There are {string_df[gene_identifier].nunique()} hits with FDR < {round(fdr_th, 2)} and absolute LFC > {lfc_low}")


        string_api_url = "https://version-11-5.string-db.org/api"
        output_format = 'tsv-no-header'
        method = 'get_link'
        my_genes = list(string_df[gene_identifier].unique())
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

            network_url = network
            if network_url.status_code == 200:
                st.markdown(f"[Link to STRING network]({network_url.text.strip()})")
                sleep(1)
            elif network_url.status_code == 400:
                st.markdown(f"STRING does not recognize unique gene identifier provided")
            else:
                st.markdown(f"HTTP request error")


app()

