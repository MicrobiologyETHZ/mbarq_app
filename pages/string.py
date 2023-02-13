import streamlit as st

if 'results_df' in st.session_state.keys():
    st.write(st.session_state['results_df'].head())
else:
    st.stop()