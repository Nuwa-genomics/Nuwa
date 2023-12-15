# Nuwa

[![Stars](https://img.shields.io/github/stars/ch1ru/nuwa?logo=GitHub&color=yellow)](https://github.com/ch1ru/nuwa/stargazers)
![CI](https://github.com/ch1ru/nuwa/actions/workflows/run_tests.yml/badge.svg?branch=main)

A bioinformatics web tool built with scanpy for genomics data processing and analysis.

## What it does?

Nuwa aims to integrate several deep learning models in a visual, easy to use interface with other filtering and data analysis familiar to most scanpy users.

## Quick start

```bash
#clone repo
git clone https://github.com/ch1ru/Nuwa.git && cd Nuwa
#If you have a Nvidia GPU
docker-compose -f cuda.docker-compose.yml up -d --build
#Or if you have a CPU
docker-compose -f cpu.docker-compose.yml up -d --build
```

## Usage

Navigate to http://localhost to use the app. Files are mounted to the streamlit-volume directory in the installation path.

## Features

- Shortcodes (Toasts card, mermaid)
- Pages Plugins (emoji, gist, avatar, mentions)
- Auto generate sidebar
- [Attribute List Definitions](https://kramdown.gettalong.org/syntax.html#attribute-list-definitions) (Primer/css utilities, Font Awesome 4)
- Service worker (caches)
- SEO (404, robots.txt, sitemap.xml)
- Canonical Link (Open Graph, Twitter Card, Schema data)


## License

Nuwa is under MIT license, a short and simple permissive license with conditions only requiring preservation of copyright and license notices. Licensed works, modifications, and larger works may be distributed under different terms and without source code. [See licence](https://github.com/ch1ru/Nuwa/blob/main/LICENSE)
