# Stage: build
FROM python:3.9-slim-bullseye AS build
COPY . /opt/labelmerge/
RUN cd /opt/labelmerge \
    && pip install --prefer-binary --no-cache-dir \
        poetry==1.2.2 \
    && poetry build -f wheel
    
# Stage: runtime
# NOTE: g++ required to install wheel (snakebids)
FROM python:3.9-slim-bullseye AS runtime
COPY --from=build /opt/labelmerge/dist/*.whl /opt/labelmerge/
RUN apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
      g++=4:10.2.1-1 \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && WHEEL=`ls /opt/labelmerge | grep whl` \
    && pip install /opt/labelmerge/$WHEEL \
    && rm -r /opt/labelmerge \
    && apt-get purge -y -q g++ \
    && apt-get --purge -y -qq autoremove
ENTRYPOINT ["labelmerge"]
