# Spatiotemporal Agricultural Drought Monitoring in Lampung Province (2015–2025)

[![Google Earth Engine](https://img.shields.io/badge/Google_Earth_Engine-GEE-34A853?style=flat&logo=googleearth)](https://earthengine.google.com/)
[![Python](https://img.shields.io/badge/Python-3.10+-3776AB?style=flat&logo=python&logoColor=white)](https://www.python.org/)
[![MODIS](https://img.shields.io/badge/Dataset-MODIS_VCI-FF6F00)](https://modis.gsfc.nasa.gov/)
[![CHIRPS](https://img.shields.io/badge/Dataset-CHIRPS_Rainfall-008080)](https://www.chc.ucsb.edu/data/chirps)

An end-to-end remote sensing and geospatial data analytics workflow for monitoring agricultural drought dynamics across Lampung Province, Indonesia. This project integrates atmospheric precipitation deficit metrics (**CHIRPS SPI-3**) with high-resolution vegetation condition metrics (**MODIS VCI**) using cropland masking and dual temporal-spatial gap-filling techniques.

---

## Executive Summary

Lampung Province is a vital agricultural hub in Indonesia. However, intense climate phenomena such as El Niño severely disrupt local rice production and food security. 

This project establishes a robust geospatial framework combining **Google Earth Engine (GEE)** for cloud-native spatial data extraction and **Google Colab (Python)** for temporal correlation analysis, sign-inversion corrections, and publication-ready cartographic outputs.

---

## System Architecture & Workflow

```text
[CHIRPS Daily (1981-2025)] ──► GEE: Thom's MLE Gamma Fitting ──► SPI-3 Climatology & Z-Scores
                                                                         │
[MODIS MOD13Q1 (250m)]    ──► GEE: Dual Gap-Filling & Cropland ──► Monthly VCI Stack
                                    Masking (MCD12Q1)                    │
                                                                         ▼
                                                          Google Colab (Python)
                                                           ├── Data Integration & Inversion Fix
                                                           ├── Lag Correlation (0–3 Months)
                                                           ├── Dual-Axis Time-Series Plotting
                                                           └── Spatial Raster Hazard Mapping
