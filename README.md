## perNMF
> **Work in progress**

This repository contains the R implementation and reproducibility code for **Personalized Nonnegative Matrix Factorization (perNMF)**.

**perNMF** is designed for nonnegative data collected from multiple individuals or sources, where both shared and source-specific latent structures may be present.

## Data

The raw dataset is not distributed with this repository.

The real-data analysis expects preprocessed user-level matrices withthe following structure:
- one matrix per user,
- rows = image labels,
- columns = posts,
- entries = label counts across images within each post.

