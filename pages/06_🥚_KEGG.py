import streamlit as st
from scripts.datasets import ResultDataSet, KeggMapsDataset
from pathlib import Path
st.set_page_config(layout='wide')
import requests
from time import sleep
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
import pandas as pd
import plotly.express as px


def app():
    st.markdown(""" ## Visualize fitness results with KEGG pathways """)
    with st.expander('How this works: '):
        an_url = "https://mbarq.readthedocs.io/en/latest/analysis.html"
        st.markdown(f"""
        #### Fitness data: 
        - a `csv` file produced by `mbarq analyze` command. To learn more about how to use `mbarq analyze`, please read [here]({an_url}).
        - First column must be a gene identifier (for example, locus tag). 
        - Must also include `LFC` and `contrast` columns, where `LFC` is log2 fold change in gene abundance for a specific treatment compared to control, and `contrast` specifies the treatment.  
        - If you organism has KEGG annotation, you can provide a three letter organism identifier and load any of the KEGG metabolic maps available.
        - You can choose a metabolic pathway of interest, and look at LFC of genes in that pathway. Genes identified as hits will have a * next to their name.
        - Make sure that the gene identifier you used for analysis is recognized by KEGG. 
        - If you load the library map on the **Data Upload** page, you would be able to choose which identifier to use for KEGG (for example, by default library map will have Name, locus tag and ID). 
        
        """)

    with st.container():
        # Get the data
        if 'results_ds' in st.session_state.keys():
            rds = st.session_state['results_ds']
        else:
            st.info('Browse example results file or upload your data in **⬆️ Data Upload**')
            result_files = [Path("examples/example_rra_results.csv")]
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
        organism_id = st.text_input('Enter 3 letter organism code', value='sey')
        kegg_options = [c for c in rds.results_df.columns if 'LFC' not in c and 'fdr' not in c]
        try:
            kix = kegg_options.index('locus_tag')
        except ValueError:
            kix = 0
        kegg_id = st.selectbox('Column corresponding to KEGG Entry names (usually locus_tag)',
                                     options=kegg_options, index=kix)
        st.info(
            f"❗Make sure KEGG recognizes unique gene identifier (you've entered `{kegg_id}` for taxon `{organism_id}`) ")

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
        fdr_th = fdr_col.number_input('FDR cutoff', value=0.05)
        type_lfc_th = lfc_col1.radio('Absolute LFC cutoff or define range', ['Absolute', 'Range'])
        if type_lfc_th == 'Absolute':
            lfc_low = lfc_col2.number_input('Log FC cutoff (absolute)', min_value=0.0, step=0.5, value=1.0)
            lfc_hi = None
        else:
            lfc_low = lfc_col2.number_input('Min Log FC', step=0.5, value=-5.0)
            lfc_hi = lfc_col2.number_input('Max Log FC', step=0.5, value=-1.0)
        rds.identify_hits(library_to_show, lfc_low, lfc_hi, fdr_th)
        kegg_df = rds.hit_df[rds.hit_df[rds.contrast_col] == contrast_to_show].copy()
        kmd = KeggMapsDataset(kegg_id, organism_id, kegg_df, rds.gene_id)
        kmd.get_gene_to_pathway_dict()

        with st.spinner(f"Loading the list of all KEGG pathways for {organism_id}"):
            pathway_map = kmd.get_org_kegg_pathways()
        pathway_description = st.selectbox('Select KEGG Pathway to explore', pathway_map.keys())
        pathway_name = pathway_map[pathway_description]
        numeric = True if st.checkbox("Display locus numbers only") else False
        if st.button("Draw map"):
            with st.spinner(f'Drawing {pathway_name} for {contrast_to_show}'):
                pathway_gene_names = kmd.display_kegg_map(pathway_name, f"{pathway_name}-{contrast_to_show}", numeric)
            st.subheader(pathway_description.split(":")[1])
            fig = rds.display_pathway_heatmap(pathway_gene_names, kegg_id)
            st.plotly_chart(fig, use_container_width=True)



app()