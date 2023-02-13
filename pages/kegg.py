import streamlit as st

st.markdown("# KEGG Integration")
if 'count_df' in st.session_state.keys():
    df = st.session_state['count_df']
    st.write(df.head())