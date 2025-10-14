FROM lopesi2bc/surfmap:2.2.0

WORKDIR /surfmap
COPY . /surfmap/

# Do the editable install as root so we can overwrite the preinstalled package
USER root
RUN /surfmap/.env/bin/pip uninstall -y surfmap || true
RUN /surfmap/.env/bin/pip install -e . --no-deps

# Add matplotlib for 3D plots; use headless backend in containers
RUN /surfmap/.env/bin/pip install matplotlib
ENV MPLBACKEND=Agg

# QoL: make the venv's python first on PATH
ENV PATH="/surfmap/.env/bin:${PATH}"

# OPTIONAL: switch back to the original non-root user if the image has one. Many images use a user named 'surfmap'. If this line fails during build, just remove it and keep running as root (you can still map UID/GID at runtime).
USER surfmap