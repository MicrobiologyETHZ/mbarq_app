import streamlit as st
import base64
import pandas as pd
import matplotlib
import seaborn as sns
from Bio import SeqIO
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
import numpy as np
#from fpdf import FPDF


def displayPDF(file):
    # Opening file from file path
    with open(file, "rb") as f:
        base64_pdf = base64.b64encode(f.read()).decode('utf-8')
    # Embedding PDF in HTML
    pdf_display = F'<iframe src="data:application/pdf;base64,{base64_pdf}" width="700" height="1000" type="application/pdf"></iframe>'
    # Displaying File
    st.markdown(pdf_display, unsafe_allow_html=True)


def display_kegg_selector(df, kegg_col='KEGG_Pathway', kegg_names_file="./examples/20-10-22-kegg-pathway-list-ko.csv"):
    if kegg_col not in df.columns:
        st.write("No KEGG Pathway information found.")
        return 'None'
    avail_pathways = [k.split(',') for k in df[kegg_col].unique()]
    avail_pathways = set([s for k in avail_pathways for s in k if s.startswith('ko')])
    kegg_names = pd.read_csv(kegg_names_file, index_col=0, header=None).to_dict()[1]
    avail_pathways = [f"{i}:{kegg_names[i]}" if i in kegg_names.keys() else i for i in avail_pathways]
    return st.selectbox('Choose KEGG Pathway to display', ["None"] + avail_pathways)





def app():
    result_file = st.file_uploader('Upload results file', accept_multiple_files=False)
    if result_file:
        data = pd.read_csv(result_file)
        data['hit'] = ((abs(data.LFC) > 1) & (data.fdr < 0.05))
        hitSummary = (data.groupby(['Name', 'contrast']).hit.sum()
                      .reset_index()
                      .rename({'hit': 'hitSum'}, axis=1))
        data = (data.merge(hitSummary, on=['Name', 'contrast'], how='outer'))
        data['hitStar'] = data['hitSum'].apply(lambda x: '*' if x > 0 else '')
        data['NameForMap'] = data['Name'] + data['hitStar'] + " (" + data['LFC_median'].round(2).astype(str) + ")"
        data['KEGG_Pathway'] = data["KEGG_Pathway"].fillna('-')
        pathwayName = display_kegg_selector(data)
        if pathwayName != 'None':
            pathwayName = pathwayName.split(":")[0]
            day = st.selectbox('Select contrast', list(data['contrast'].unique()))
            pDf = data[(data.KEGG_Pathway.str.contains(pathwayName)) & ((data.contrast == day) | (data.contrast.isnull()))].copy()
            # minima = -max(abs(pDf.LFC_median.min()), abs(pDf.LFC_median.max()))
            # maxima = max(abs(pDf.LFC_median.min()), abs(pDf.LFC_median.max()))
            norm = matplotlib.colors.Normalize(vmin=-6, vmax=6, clip=True)
            mapper = matplotlib.cm.ScalarMappable(norm=norm, cmap=sns.diverging_palette(220, 20, as_cmap=True))
            pDf['hex'] = pDf.LFC_median.apply(mapper.to_rgba).apply(matplotlib.colors.to_hex)
            ko_map = (pDf[['Name', 'LFC_median', 'KEGG_ko', 'hex', 'NameForMap']]
                      .replace('-', np.nan)
                      .dropna(subset=['KEGG_ko'])
                      .drop_duplicates())
            new_cols = ko_map.KEGG_ko.str.split(",", expand=True)
            ko_map = pd.concat([ko_map, new_cols], axis=1)
            ko_map = (ko_map.melt(id_vars=['LFC_median', 'NameForMap', 'hex'], value_vars=new_cols.columns, value_name='KO')[
                          ['NameForMap', 'LFC_median', 'KO', 'hex']]
                      .dropna(subset=['KO']))
            ko_map['hex'] = ko_map.hex.str.replace('#000000', "#7c7d83")
            ko_dict = ko_map.set_index("KO").to_dict()

            pathway = KGML_parser.read(kegg_get(pathwayName, "kgml"))
            canvas = KGMLCanvas(pathway, import_imagemap=True)

            for element in pathway.orthologs:
                color = None
                name = None
                node_kos = element.name.split()
                for ko in node_kos:
                    color = ko_dict['hex'].get(ko, color)
                    name = ko_dict['NameForMap'].get(ko, name)
                for graphic in element.graphics:
                    if color is not None:
                        graphic.bgcolor = color
                        graphic.name = name
                    else:
                        graphic.bgcolor = '#FFFFFF'
            c1, c2 = st.columns(2)
            st.markdown("⚠️ Displaying the map is only possible in Firefox")
            fname = f"{pathwayName}_{day}_map.pdf"
            canvas.draw(fname)
            if c1.button(f'Display {pathwayName} map'):
                displayPDF(fname)

            with open(fname, "rb") as f:
                c2.download_button(
                    f"Download {pathwayName} map",
                    data=f,
                    file_name=fname,
                )




        # kegg_col = 'KEGG_Pathway'
        # if kegg_col not in subset_df.columns:
        #     st.write("No KEGG Pathway information found.")
        #     kegg_to_show = 'None'
        #     kegg_df = subset_df.copy()
        # else:
        #     st.markdown('#### Analyze KEGG Maps')
        #     kegg_to_show = display_kegg_selector(subset_df, kegg_col)
        #     if kegg_to_show != 'None':
        #         kegg_df = subset_df[subset_df[kegg_col].str.contains(kegg_to_show.split(":")[0])]
        #         pathwayName = kegg_to_show.split(":")[0]
        #         contrast = contrasts[0]
        #         c1, c2 = st.columns(2)
        #         if c1.button(f'Display {pathwayName} map'):
        #             fname = display_kegg_map(kegg_df, pathwayName, contrast )
        #             displayPDF(fname)
        #             c1.markdown("⚠️ Displaying the map is only possible in Firefox")
        #             with open(fname, "rb") as f:
        #                 c2.download_button(
        #                     f"Download {pathwayName} map",
        #                     data=f,
        #                     file_name=fname,
        #                 )
        #     else:
        #         kegg_df = subset_df.copy()
        #         st.write("No KEGG pathway selected")