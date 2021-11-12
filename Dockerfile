# The build-stage image:
FROM continuumio/miniconda3

RUN conda install python=3.7 -y
RUN conda install rdkit -c rdkit -y
RUN conda install xtb-python -c conda-forge -y

COPY requirements.txt .

COPY xtbservice ./xtbservice
RUN pip install --no-cache-dir -r requirements.txt

COPY README.md .

CMD gunicorn -w $WORKERS xtbservice.xtbservice:app -b 0.0.0.0:$PORT -k uvicorn.workers.UvicornWorker -t $TIMEOUT --keep-alive $TIMEOUT
