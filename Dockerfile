FROM python:3.9-slim

LABEL description="methylink"
LABEL maintainer="Mei Wu"

ENV WORK_DIR=/work

RUN apt-get update && apt-get install -y \
    zlib1g-dev libcurl4-gnutls-dev liblzma-dev libncurses5-dev libncursesw5-dev libssl-dev libbz2-dev

COPY . ${WORK_DIR}

WORKDIR ${WORK_DIR}

RUN python -m pip install .