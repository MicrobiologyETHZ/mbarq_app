FROM python:3.9-slim

WORKDIR /app

RUN apt-get update

COPY requirements.txt requirements.txt

RUN pip3 install -r requirements.txt

copy . .

EXPOSE 8501

ENTRYPOINT ["streamlit", "run", "Home.py", "--server.port=8501", "--server.address=0.0.0.0"]