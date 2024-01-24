# mBARq app

Streamlit application to run a versatile and user-friendly framework for the analysis and interpretation of RB-TnSeq and other barcoded sequencing data.

# Changes for deploy:
- Versions relaxed in the requirements file:
    - altair == 4.2
    - streamlit == 1.23.1
    - biopython == 1.79

- The old requirements file is preserved in `requirements_bu20240115.txt`

# Instructions to deploy in the server

- Access the micro-flask server `micro-flask-dev.ethz.ch`
- Activate mbarq conda environment, installed using the `requirements.txt` file
```
conda activate mbarq
```
- Run streamlit specifying the port
```
streamlit run streamlit_app.py --server.port 55556
```
- Run it in a screen session to keep it open