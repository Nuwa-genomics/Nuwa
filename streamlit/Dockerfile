FROM nvidia/cuda:11.2.2-cudnn8-devel-ubuntu20.04

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    python3.8 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . .
ENV STREAMLIT_SERVER_MAX_UPLOAD_SIZE=100000
RUN pip3 install -r requirements.txt
CMD ["streamlit", "run", "Dashboard.py", "--server.port=8501", "--server.address=0.0.0.0"]