FROM jupyter/minimal-notebook:latest

# Copy required files
COPY start.sh .
COPY environment.yml .

# Switch to root to install packages & configure permissions
USER root

# Set required permissions
RUN chmod +x start.sh

# Install system libraries required for python and R installations
RUN apt-get update && apt-get install -y --no-install-recommends build-essential apt-utils ca-certificates zlib1g-dev gfortran locales libxml2-dev libcurl4-openssl-dev libssl-dev libzmq3-dev libreadline6-dev xorg-dev libcairo-dev libpango1.0-dev libbz2-dev liblzma-dev libffi-dev libsqlite3-dev libopenmpi-dev libhdf5-dev libjpeg-dev libblas-dev liblapack-dev libpcre2-dev libgit2-dev libgmp-dev

# Install common linux tools
RUN apt-get update && apt-get install -y --no-install-recommends wget curl htop less nano vim emacs git

# Configure default locale
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8

# Switch back to non-privileged user
USER jovyan

# Create the environment
RUN conda env create -f environment.yml && conda clean -a

# Activate the environments
RUN echo "source activate {{  cookiecutter.project_name }}" >> ~/.bashrc
ENV PATH /home/user/miniconda/envs/{{ cookiecutter.project_name }}/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name {{ cookiecutter.project_name }} > {{ cookiecutter.project_name }}_environment.yml
