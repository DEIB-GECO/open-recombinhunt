FROM python:3.11.13-slim-bookworm

WORKDIR /app

RUN apt-get update && apt-get install -y build-essential
RUN pip install --upgrade pip

COPY ./requirements.txt /app/requirements.txt

RUN pip install -r requirements.txt

COPY ./ /app/

EXPOSE 8501

HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

ENTRYPOINT ["python", "-u", "-m", "streamlit", "run", "app/streamlit_app.py", "--server.enableCORS=false","--server.enableXsrfProtection=false","--server.enableWebsocketCompression=true"]
