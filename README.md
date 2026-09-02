## Nonnegative Matrix Factorization with Shared and Individual-Specific Components
> **Work in progress**

This repository contains the R implementation and reproducibility code for **Nonnegative Matrix Factorization with Shared and Individual-Specific Components**.

This method is designed for nonnegative data collected from multiple individuals or sources, where both shared and source-specific latent structures may be present.

## Data

The raw dataset is not distributed with this repository.

The real-data analysis expects preprocessed user-level matrices withthe following structure:
- one matrix per user,
- rows = image labels,
- columns = posts,
- entries = label counts across images within each post.

## Reproducibility
### Simulation
The simulation workflow consists of three steps:

```text
simul_train.r
    ↓
simul_test.r
```

Run all commands from the root directory of the repository.

### Real-Data Analysis

The real-data workflow consists of two main steps:

```text
real_train.r
     ↓
real_test.r
```

Run the scripts from the root directory of the repository.

## Directory and codes

```
.
+-- real
|        +-- plot_recon_real.png
+-- simul
|        +-- plot_recon_simul.png
+-- README.md
+-- generate.r
+-- real_test.r
+-- real_train.r
+-- simul_test.r
+-- simul_train.r
+-- util.r
```

## Acknowledgement
This repository was developed with support from the 서울시립대학교 데이터 사이언스 플러스 차세대 융합인재 양성사업단 - http://dsplus.uos.ac.kr/
