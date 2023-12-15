---
sort: 1
---

# Installation

## 1. Install Docker Compose for command line 🐋

See the [tutorial](https://docs.docker.com/compose/install/) to install for your platform.

```tip
Deleting docker containers and images won't affect any saved data in the 'streamlit-volume' directory.

However, it's recommended to back up saved data within this directory.
```

## 2. Bring up containers

First, pull the repo from github:
```bash
git clone https://github.com/ch1ru/nuwa
```

If you have a Nvidia GPU:
```bash
docker-compose -f cuda.docker-compose.yml up -d --build
```

Otherwise:
```bash
docker-compose -f cpu.docker-compose.yml up -d --build
```
The web app can now be accessed at http://localhost

```tip
By default the web server serves on port 80. This can be changed in the .env file if this conflicts with any existing service.
```

## 3. Stopping containers

To stop containers simply run:
```bash
docker-compose -f cuda.docker-compose.yml down #use the yaml file you used to build containers
```

```tip
To see running containers run `docker ps`
```

