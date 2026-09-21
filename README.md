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
(국문) 이 논문은 정부(과학기술정보통신부)의 재원으로 한국연구재단의 지원을 받아 수행된 연구임 (NO. RS-2022-NR068754)

(영문) This research was supported by the National Research Foundation of Korea(NRF) grant funded by the Korea government(MSIT) (NO. RS-2022-NR068754)
