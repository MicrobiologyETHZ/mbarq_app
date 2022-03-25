import streamlit as st
from PIL import Image

from pathlib import Path

# Custom imports
from multipage import MultiPage
from pages import Home, Library, Exploratory, DiffAb, Pathway

st.set_page_config(page_title="mBARq App", layout='wide',
                   page_icon=Image.open("images/image.png")
                   )
#st.title("mBARq Data Visualization and Mining")
st.image("images/mbarq-logo.png")
# Create an instance of the app
app = MultiPage()

pages = {'Home': ('Home', Home.app),
         'Library': ('Library Map', Library.app),
         'Exploratory': ('Exploratory Analysis', Exploratory.app),
         'DiffAb': ('Differential Abundance', DiffAb.app),
         'Pathway': ('Metabolic Pathways', Pathway.app)}

for page_name, page in pages.items():
    app.add_page(page[0], page[1])

app.run()
