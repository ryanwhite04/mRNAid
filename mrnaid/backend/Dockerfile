FROM python:3.9

RUN apt-get update -y \
        && apt-get install -y --no-install-recommends \
            libgtextutils-dev libtbb-dev \
            autoconf automake libcurl4-gnutls-dev libncurses5-dev \
            build-essential pkg-config \
        && apt-get -y clean && apt-get -y autoclean && apt-get -y autoremove


# Set environment variables to prevent prompts
ENV DEBIAN_FRONTEND=noninteractive

# Create a non-root user
RUN useradd -m -s /bin/bash docker
# -m means create a home directory

ENV PATH="/home/docker/.local/bin:${PATH}"

# Create a new directory for the app and set ownership
WORKDIR /usr/src/app
RUN chown docker:docker /usr/src/app

# Copy the files into the container
COPY --chown=docker:docker . .
USER docker
RUN pip install -r ./requirements.txt
# Switch to the non-root user

# open port 5000
EXPOSE 5000
# Your CMD or ENTRYPOINT instruction would go here
# CMD ["python", "your_script.py"]
