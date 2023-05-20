FROM python:3.8-alpine

LABEL description="methylink"
LABEL maintainer="Mei Wu"

ARG WORK_DIR=project

WORKDIR ${WORK_DIR}
COPY . /${WORK_DIR}

RUN pip install --upgrade --no-cache-dir .