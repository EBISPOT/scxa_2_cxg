
FROM python:3.11-buster

RUN pip install poetry

RUN apt-get update && apt-get install -y graphviz graphviz-dev

ADD . /opt/
RUN cd /opt && poetry install

WORKDIR /opt

CMD ["bash", "-c", "poetry run python src/bulk_experiments.py --study_filter E --chunk_size 10 --download"]


