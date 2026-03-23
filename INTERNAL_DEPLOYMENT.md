# Internal Deployment Guide (mbarq.microbiomics.io)

This document describes the current internal deployment flow for `mbarq_app`, based on:

- local script: `buildDockerToServer.sh`
- server directory structure under `/science/instances/mbarq.microbiomics.io`

## Overview

`buildDockerToServer.sh` does **three** things:

1. Builds Docker image locally as `mbarq_app_amd:0.1`
2. Saves image into a jar file named `mbarq_app_$USER.jar`
3. Uploads both files to the server path:
   - jar file
   - redeploy script (`redeploy-mbqrq.sh` from local repo)

After upload, redeploy still needs to be executed on the server.

## Prerequisites

- Local machine:
  - Docker installed and working
  - SSH access to `monod-biol`
  - `scp` available
- Server:
  - target directory exists: `/science/instances/mbarq.microbiomics.io`
  - user has permissions for files in that directory
  - user can run `sudo docker ...`

## Deploy Steps

### 1) Run deployment build/upload script locally

From repo root:

```bash
chmod +x buildDockerToServer.sh
./buildDockerToServer.sh
```

This uploads:

- `mbarq_app_$USER.jar`
- `redeploy-mbqrq.sh` (local file name)

to:

- `/science/instances/mbarq.microbiomics.io`

### 2) Login server and go to deployment directory

```bash
ssh "$USER"@monod-biol
cd /science/instances/mbarq.microbiomics.io
ls -lah
```

Expected files include:

- `mbarq_app_$USER.jar`
- redeploy shell script (see naming note below)

### 3) Execute redeploy script on server

If the script on server is named `redelopy-mbarq.sh`:

```bash
chmod +x redelopy-mbarq.sh
./redelopy-mbarq.sh
```

Current script content used on server:

```bash
sudo docker stop $(sudo docker ps -a -q --filter name=mbarq_app)
sudo docker rm $(sudo docker ps -a -q)
sudo docker rmi $(sudo docker images -a -q)
sudo docker load < mbarq_app_$USER.jar
sudo docker run -d -p 55556:8501 --name mbarq_app mbarq_app_amd:0.1
```

### 4) Verify service

Basic container check:

```bash
sudo docker ps | rg mbarq_app
```

App endpoint check (from any host with access):

```bash
curl -I https://mbarq.microbiomics.io/
```

## Naming Notes (Important)

There are currently two similar script names:

- local upload script references: `redeploy-mbqrq.sh`
- server example file name: `redelopy-mbarq.sh`

Use one consistent name to reduce operational risk. Recommended:

- standardize to `redeploy-mbarq.sh` in both local repo and server path.

## Rollback (Simple)

If a new deployment fails:

1. Keep previous jar files with timestamped names (recommended for future updates)
2. Load previous jar:

```bash
sudo docker load < <previous_jar_file>.jar
```

3. Re-run redeploy script.

