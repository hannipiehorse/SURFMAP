FROM lopesi2bc/surfmap:2.2.0

WORKDIR /surfmap
COPY . /surfmap/

USER root
RUN /surfmap/.env/bin/pip uninstall -y surfmap || true
RUN /surfmap/.env/bin/pip install -e . --no-deps

# Combine Python installs to reduce layers
RUN /surfmap/.env/bin/pip install matplotlib "pdb2pqr==3.6.2"
ENV MPLBACKEND=Agg

ENV PATH="/surfmap/.env/bin:${PATH}"
USER surfmap
