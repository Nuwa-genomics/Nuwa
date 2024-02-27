# Installation

## 1. Install Docker Compose 🐋

See the [tutorial](https://docs.docker.com/compose/install/) to install for your platform.

```tip
Deleting docker containers and images won't affect any saved data in the 'streamlit-volume' directory. However, it's recommended to back up saved data within this directory.
```

## 2. Clone repo to get the latest version

```bash
git clone https://github.com/nuwa-genomics/nuwa
```

Alternatively download and extract the zip file of a specific version from the [releases page](https://github.com/nuwa-genomics/Nuwa/releases).

## 3. Bring up containers

```note
## If using a GPU
If you are planning to use a GPU for faster training:
- Make sure cuda drivers are installed on the host machine.
- Install and configure [Nvidia container toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html) for docker
```

**Then bring up containers:**

for GPU:
```bash
docker-compose -f cuda.docker-compose.yml up -d --build
```

Or for CPU only:
```bash
docker-compose -f cpu.docker-compose.yml up -d --build
```

**The web app can now be accessed at http://localhost**

```tip
By default the web server runs on port 80. This can be changed in the .env file if this conflicts with any existing service.
```

## 3. Stopping containers

To stop containers simply run:
```bash
docker-compose -f cuda.docker-compose.yml down #use the yaml file you used to build containers
```

```tip
To see running containers run `docker ps`
```

