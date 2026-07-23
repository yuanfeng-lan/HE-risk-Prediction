
# HE-risk-Prediction

## XGBoost-SHAP Interpretable Model for Hepatic Encephalopathy Risk Prediction in Cirrhosis

---

### 📖 Project Overview

This repository contains the complete analysis code and models for the manuscript:

**"XGBoost-SHAP Interpretable Modeling Identifies and Validates an Eight-Gene Biomarker for Hepatic Encephalopathy Risk Prediction in Cirrhosis"**

The project integrates four cirrhotic transcriptomic cohorts (GSE41919, GSE57193, GSE139602, GSE15654) to develop an interpretable machine learning model for early risk stratification of hepatic encephalopathy (HE) in cirrhotic patients. The final model is based on an eight-gene signature (*PRB2*, *TUBA1C*, *NPC2*, *LRRC32*, *TLN1*, *SOX9*, *SERPINA3*, and *RNASE4*) and incorporates both progression and prognostic risk prediction.

---

### 📁 Repository Structure

```
HE-risk-Prediction/
├── 00_rawdata/               # Raw data files (GEO cohorts)
├── 06_RF/                    # Random forest feature selection
├── 07_xgboost/               # XGBoost model training & SHAP analysis
├── 09_SEM/                   # Structural equation modeling & mediation analysis
├── 10_GSEA/                  # Gene set enrichment analysis
├── models/                   # Pre-trained XGBoost models
│   ├── final_xgboost_model_GSE139602.Rdata
│   └── final_xgboost_model_GSE15654_cox.Rdata
├── renv.lock                 # R package environment lock file
├── session_info.log          # R session information
├── README.md                 # This file
└── LICENSE                   # MIT License
```

---

### 🔧 Requirements & Environment Setup

#### Prerequisites

- **R** (version ≥ 4.3.0)
- **renv** package for environment management

#### Quick Setup

Clone the repository and restore the R environment:

```bash
git clone https://github.com/yuanfeng-lan/HE-risk-Prediction.git
cd HE-risk-Prediction
```

Then, in R:

```r
# Install renv if not already installed
install.packages("renv")

# Restore the exact package versions used in this study
renv::restore()
```

This will install all required R packages with their exact versions as recorded in `renv.lock`.

Alternatively, you can view the session information:

```r
# Load the saved session info
readLines("session_info.log")
```

---

### 📊 Workflow Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                         Data Processing                        │
├─────────────────────────────────────────────────────────────────┤
│  GSE41919 (HE cohort)    │    GSE57193 (HE cohort)             │
│  GSE139602 (Progression) │    GSE15654 (Prognosis)             │
└──────────────────────────┴─────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Feature Identification                      │
├─────────────────────────────────────────────────────────────────┤
│  • Differential expression analysis (limma)                    │
│  • Mfuzz trajectory analysis                                   │
│  • Prognosis-associated DEGs                                   │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│              Machine Learning Feature Selection                │
├─────────────────────────────────────────────────────────────────┤
│  • LASSO Regression     • RFE (Recursive Feature Elimination)  │
│  • Random Forest        → 8 final marker genes                 │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    XGBoost Model Development                   │
├─────────────────────────────────────────────────────────────────┤
│  Model 1: Progression Model (GSE139602, 5 genes)               │
│  Model 2: Prognostic Model (GSE15654, 3 genes)                 │
│  • SHAP for interpretability                                   │
│  • AUC-weighted fusion → integrated risk score                 │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                  Functional & Mechanistic Analysis             │
├─────────────────────────────────────────────────────────────────┤
│  • SEM / Mediation analysis                                    │
│  • GSEA with gradient permutation testing (1k/10k/100k)        │
│  • Functional characterization (metabolic/synaptic/immune)     │
└─────────────────────────────────────────────────────────────────┘
```

---

### 🚀 Running the Analysis

All analysis scripts are organized sequentially:

Run each script sequentially in R:

```r
source("00_data_preprocessing.R")
source("06_RF_feature_selection.R")
# ... and so on
```

---

### 🧠 Models

Two pre-trained XGBoost models are available in the `models/` directory:

| Model | Cohort | Features | Objective |
|-------|--------|----------|-----------|
| `final_xgboost_model_GSE139602.Rdata` | GSE139602 | 5 genes (PRB2, TUBA1C, LRRC32, NPC2, TLN1) | Progression prediction |
| `final_xgboost_model_GSE15654_cox.Rdata` | GSE15654 | 3 genes (SOX9, SERPINA3, RNASE4) | Survival prognosis |

To load a model:

```r
model <- readRDS("models/final_xgboost_model_GSE139602.Rdata")
```

---

### 📝 Key Methods Summary

| Component | Method | Package |
|-----------|--------|---------|
| Feature selection | LASSO, RFE, Random Forest | glmnet, caret, randomForest |
| Classification | XGBoost | xgboost |
| Interpretability | SHAP | shapviz |
| Differential expression | limma | limma |
| Trajectory analysis | Mfuzz | Mfuzz |
| SEM | Structural Equation Modeling | lavaan |
| Mediation | Bootstrap (2000 resamples) | lavaan / mediation |
| GSEA | fgsea with gradient permutation | fgsea |
| Survival analysis | Kaplan-Meier, log-rank | survival, survminer |

---

### 📂 Data Sources

All datasets used are publicly available from the Gene Expression Omnibus (GEO):

| Dataset | Description | Reference |
|---------|-------------|-----------|
| GSE41919 | HE vs. non-HE cirrhosis (n=11) | PMID: 31571389 |
| GSE57193 | HE vs. non-HE cirrhosis (n=8) | PMID: 31571389, 25064044 |
| GSE139602 | Progression stages (n=39) | PMID: 35540106, 40437093 |
| GSE15654 | Prognosis cohort (n=216) | PMID: 23333348, 28256512, 31344396 |

---


---

### 🙏 Acknowledgments
