FROM nvidia/cuda:12.2.0-runtime-ubuntu22.04

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    python3.10 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . .
ENV STREAMLIT_SERVER_MAX_UPLOAD_SIZE=100000
RUN pip3 install -r requirements.txt
RUN pip3 install pyg_lib torch-scatter torch-sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-2.1.0+cu121.html
CMD ["streamlit", "run", "Dashboard.py", "--server.port=8501", "--server.address=0.0.0.0"]