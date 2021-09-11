# The build-stage image:
FROM continuumio/miniconda3 AS build

# Install the package as normal:
COPY environment.yml .
RUN conda env create -f environment.yml

# Install conda-pack:
RUN conda install -c conda-forge conda-pack

COPY install_packages.sh .
RUN ./install_packages.sh

COPY requirements.txt .

COPY xtbservice ./xtbservice

COPY README.md .

RUN pip install --no-cache-dir -r requirements.txt

# Use conda-pack to create a standalone enviornment
# in /venv:
RUN conda-pack -n xtbservice -o /tmp/env.tar && \
  mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
  rm /tmp/env.tar

# We've put venv in same path it'll be in final image,
# so now fix up paths:
RUN /venv/bin/conda-unpack


# The runtime-stage image; we can use Debian as the
# base image since the Conda env also includes Python
# for us.
FROM debian:buster AS runtime

# Copy /venv from the previous stage:
COPY --from=build /venv /venv

CMD gunicorn -w 4 xtbservice.xtbservice:app -b 0.0.0.0:$PORT -k uvicorn.workers.UvicornWorker