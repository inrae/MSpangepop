FROM python:3.11-slim

WORKDIR /app

RUN apt-get update && apt-get install -y \
    build-essential \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

COPY . /app

RUN pip install --upgrade pip
RUN pip install msprime tskit pandas h5py numpy matplotlib ipython pyyaml biopython pytest bitarray

EXPOSE 5000

ENV PYTHONUNBUFFERED=1

# Add maintainer information using LABEL
LABEL maintainer="Lucien Piat" \
      description="This image is used for the MSpangepop workflow at INRAe."
