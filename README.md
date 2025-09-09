<img width="1536" height="1024" alt="ChatGPT Image Sep 9, 2025, 03_12_04 PM" src="https://github.com/user-attachments/assets/6818e589-fb93-45ce-be07-f7f2a2b10ae0" />

# 🌿🐝 _Osmia lignaria_ & Climate Change in California

Modeling the future of a wild solitary bee and its native plant resources

## 🌍 Context

Pollinators like _Osmia lignaria_ are crucial for California’s ecosystems and agriculture. Climate change threatens both the bee and its plant resources. Understanding these dynamics helps inform conservation strategies, agricultural resilience, and climate policy in California.

## 📖 Project Description

This repository hosts the code of analyses for studying the effects of climate change on the wild solitary bee Osmia lignaria and its plant resources in California, as part of California’s Fifth Climate Change Assessment. The input data is not included as the volume is high. But all data is publicly available and reference in the main document and appendices with the respective dois (see citation). 

We use species distribution modeling (SDM), climate projections (LOCA2), and ecological trait data to:

- Estimate current and future distributions of O. lignaria and its resource plants.
- Quantify changes in climatic suitability under different Shared Socioeconomic Pathways (SSP245, SSP370, SSP585).
- Assess ecological risks to pollination services and biodiversity conservation.
🔗 Repository URI: github.com/baltazargch/Osmia_lignaria_CalAdapt

## 🗂 Repository Structure

<details>
  <summary><strong>Main structure</strong></summary>

<pre>
Osmia_lignaria_CalAdapt/
├── Python/
├── R/
├── Apendix_II.Rmd
├── LICENSE
├── Osmia_lignaria_CalAdapt.Rproj
├── README.md
├── download_paralell.sh
└── requerements.txt
</pre>
</details>   

## ⚙️ Methods

- Occurrence Data: GBIF + curated observations
- Climate Data: LOCA2 (Cal-Adapt, 3 km resolution)
- Modeling Frameworks:
  - biomod2 for ensemble SDMs
  - ENMeval for MaxEnt tuning
- Bioclimatic Variables: custom terra-compatible functions for BIO1–BIO19
- Thresholding: Minimum Training Presence (MTP) & 10th percentile presence
- Evaluation: AUC, TSS, Kappa, cross-validation

## 📊 Outputs

- 📑 Tables: Model performance metrics, percent change in suitable area
- 🗺 Maps: Current & projected distributions of O. lignaria and plants
- 📈 Figures: Climate suitability trends, resource plant vulnerabilities
- 📘 Appendices: Detailed methodological documentation

## 🚀 Getting Started

Prerequisites
- R ≥ 4.3
- Key packages: `install.packages(c("terra", "sf", "tidyverse", "biomod2", "ENMeval", "rgbif"))`

### Run the models

- Check the scripts in R/ for:
- Data preparation
- Model training and evaluation
- Ensemble forecasting
- Visualization

## ✒️ Citation

If you use this code or results, please cite:

Gaiarsa, M.P., et al. (2025). POLLINATION IN A CHANGING WORLD: THE EFFECTS OF HABITAT AND TEMPERATURE ON THE HEALTH OF CALIFORNIA'S WILD POLLINATORS. California’s Fifth Climate Change Assessment.

## 🤝 Contributing

Contributions, suggestions, and collaborations are welcome! Please open an issue or submit a pull request.

## 📬 Contact

Author: Baltazar González Chávez
Affiliation: University of California, Merced
Email: bgonzalezchavez@ucmerced.edu

## ✨ Let’s build knowledge for pollinator conservation in a changing climate.
