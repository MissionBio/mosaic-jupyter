# Vignettes

### About

This repository maintains a set of vignettes showcasing
the capabilities of the Mission Bio analysis package, [Mosaic](https://github.com/MissionBio/mosaic)

### List of vignettes
- [Basic usage of Mosaic](https://missionbio.github.io/mosaic-vignettes/basics/basics.html)

### Interactive app
A user interface for Mosaic has been built using [Streamlit](https://www.streamlit.io/)

After installing [Mosaic](https://github.com/MissionBio/mosaic), install streamlit
in the appropriate environment using:

```
pip install streamlit
```

the app can be launched using:

```
streamlit run ./streamlit/mosaic.py
```

For s3 access, credentials have to be added to `/.aws/credentials`

The downloaded h5 files must be stored under `./streamlit/h5/downloads/`
to be accessible throught the UI.

