#!/bin/bash
echo "------Build Docker------"
docker build -t mbarq_app_amd:0.1 . --platform=linux/amd64
echo "------Save Docker Image To Jar------"
docker save mbarq_app_amd:0.1 > mbarq_app_"$USER".jar
echo "------Copy Docker Image To Server. It takes time------"
scp ./mbarq_app_"$USER".jar "$USER"@micro-flask-dev.ethz.ch:/docker-run/
echo "------DONE------"