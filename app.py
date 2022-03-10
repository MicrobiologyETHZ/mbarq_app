import streamlit as st
import tarfile
from pathlib import Path
import yaml


# Custom imports
from multipage import MultiPage
from pages import PCA, DiffAb, Pathway,  Expression, Library, Home
st.set_page_config(page_title="mBARq App", layout='wide')

st.title("mBARq Data Visualization and Mining")

# Create an instance of the app
app = MultiPage()

pages = {'Home' :('Home', Home.app),
        'Library': ('Library Map', Library.app),
         'PCA': ('PCA', PCA.app),
         'Expression': ('Barcode Abundance', Expression.app),
         'DiffAb': ('Differential Abundance', DiffAb.app),
         'Pathway':('Metabolic Pathways', Pathway.app)}

for page_name, page in pages.items():
    app.add_page(page[0], page[1])





#


#app.add_page("Analyse featureCounts Data", featureCounts_analysis.app)
#app.add_page("Explore DE analysis results", explore_results.app)

# The main app
app.run()