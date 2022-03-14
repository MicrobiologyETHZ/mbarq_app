import streamlit as st
from PIL import Image

from pathlib import Path

# Custom imports
from multipage import MultiPage
from pages import Home, Library, PCA, Expression, DiffAb, Pathway

st.set_page_config(page_title="mBARq App", layout='wide',
                   #page_icon=Image.open("/Users/ansintsova/git_repos/mbarq_app/data/images/image.png")
                   )
st.title("mBARq Data Visualization and Mining")

# Create an instance of the app
app = MultiPage()

pages = {'Home': ('Home', Home.app),
         'Library': ('Library Map', Library.app),
         'PCA': ('PCA', PCA.app),
         'Expression': ('Barcode Abundance', Expression.app),
         'DiffAb': ('Differential Abundance', DiffAb.app),
         'Pathway': ('Metabolic Pathways', Pathway.app)}

for page_name, page in pages.items():
    app.add_page(page[0], page[1])

app.run()
