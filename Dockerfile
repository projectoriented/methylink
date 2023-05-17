FROM python:3.8-alpine

LABEL description="methylink"
LABEL maintainer="Mei Wu"

ARG WORK_DIR=project

COPY . /${WORK_DIR}

RUN cd /${WORK_DIR} && pip install --upgrade --no-cache-dir .