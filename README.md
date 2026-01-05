# bioscanflow

Repository for developing documentation for metabarcoding-based biomonitoring workflows within the framework of Biodiversity Genomics Europe project.
Readthedocs page: https://bioscanflow.readthedocs.io/en/latest/index.html# 
__________________________

'Dockerfile' and 'bioscanflow_env.yml' are for building [Docker image](https://registry.hub.docker.com/r/pipecraft/bioscanflow) for bioinformatics processing of the sequenceing data. 
__________________________

# for developers


Install sphinx and rtd-theme:

```bash
pip install -U sphinx
pip install sphinx-rtd-theme
pip install sphinxcontrib.youtube # for youtube videos
```

Build for testing (need to be in docs folder)

For Windows:

```bash
./make.bat html
```

For Linux and macOS:

```bash
make html
```






