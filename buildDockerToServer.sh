#!/bin/bash
set -euo pipefail

IMAGE_TAG="mbarq_app_amd:0.1"
JAR_FILE="mbarq_app_${USER}.jar"
REMOTE_HOST="${USER}@monod-biol"
REMOTE_DIR="/science/instances/mbarq.microbiomics.io"
REDEPLOY_SCRIPT="redeploy-mbqrq.sh"

echo "------Build Docker------"
docker build -t "${IMAGE_TAG}" . --platform=linux/amd64

echo "------Save Docker Image To Jar------"
docker save "${IMAGE_TAG}" > "${JAR_FILE}"

echo "------Copy Docker Image + Redeploy Script To Server (takes time)------"
scp "./${JAR_FILE}" "${REMOTE_HOST}:${REMOTE_DIR}/"
scp "./${REDEPLOY_SCRIPT}" "${REMOTE_HOST}:${REMOTE_DIR}/"

echo "------DONE------"