# {{ cookiecutter.project_name }}

A single cell analysis project.

# Building

`docker build -t schillerlab/{{ cookiecutter.project_name }}:latest`

# Usage

1. Start the container `docker run -it -p 8888:8888 schillerlab/{{ cookiecutter.project_name }}:latest /bin/bash`
2. Run bash start.sh

To run the container with root rights add `-e GRANT_SUDO=yes --user root` to the run commands.
