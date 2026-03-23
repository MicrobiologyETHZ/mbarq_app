sudo docker stop $(sudo docker ps -a -q --filter name=mbarq_app)
sudo docker rm $(sudo docker ps -a -q)
sudo docker rmi $(sudo docker images -a -q)
sudo docker load < mbarq_app_$USER.jar
sudo docker run -d -p 55556:8501 --name mbarq_app mbarq_app_amd:0.1