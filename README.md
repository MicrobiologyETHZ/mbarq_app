# mBARq app

Streamlit application to run a versatile and user-friendly framework for the analysis and interpretation of RB-TnSeq and other barcoded sequencing data.

# Changes for deploy:
- Versions relaxed in the requirements file:
    - altair == 4.2
    - streamlit == 1.23.1
    - biopython == 1.79

- The old requirements file is preserved in `requirements_bu20240115.txt`

# How to run in docker
```bash
docker build -t mbarq_app:0.1 .
docker run -p 8501:8501 mbarq_app:0.1
# check if the app is healthy
curl -i http://localhost:8501
HTTP/1.1 200 OK
Server: TornadoServer/6.4
Content-Type: text/html
Date: Fri, 22 Mar 2024 11:03:55 GMT
Accept-Ranges: bytes
Etag: "85c873604085732a292417b9cda592560938848813798b7e98d66d9233d2cb072ae85d199d602a71e2e1f49164da4f0d8d2fafb374dc2a2314fa439bdc4b76bd"
Last-Modified: Fri, 22 Mar 2024 10:41:59 GMT
Cache-Control: no-cache
Content-Length: 891
Vary: Accept-Encoding
```

# Instructions to deploy in the server

## Prerequisite
* [Docker](https://docs.docker.com/engine/install/) !!! Ensure you can do docker without sudo in your local env

## By Script

1. run the script
```bash
chmod +x buildDockerToServer.sh
./buildDockerToServer.sh
```

2. login the server and run the bash script to redeploy docker
```bash
cd /docker-run
./clean_docker_and_restart-mbarq.sh
```
3. check if https://mbarq.microbiomics.io/ is running
