FROM nvidia/cuda:12.2.0-runtime-ubuntu22.04
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*
#streamlit environment variables
ENV STREAMLIT_SERVER_MAX_UPLOAD_SIZE=100000
# install miniconda
ENV CONDA_DIR /opt/conda
ENV CONDA_ENV=bioconda_env
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py311_24.1.2-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
RUN conda config --set channel_priority strict
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
WORKDIR /app
COPY conda.yaml .
RUN conda env create -n $CONDA_ENV -f conda.yaml
RUN conda install -y -c conda-forge fa2==0.3.5 lightning==2.2.1
COPY requirements.txt .
RUN pip3 install --use-pep517 -r requirements.txt
RUN pip3 install pyg_lib torch-scatter torch-sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-2.2.1+cu121.html
COPY . .
RUN rm -rf tests
CMD ["streamlit", "run", "Dashboard.py", "--server.port=8501", "--server.address=0.0.0.0"]