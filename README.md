<div align="center">

<!-- Replace with your own banner/logo path -->
<img src="banner.svg" alt="Neuronal Network Modeling Logo">

*Modeling and Predicting Neuronal Network Behavior Using Gene Expression Data*

<!-- Badges (edit links as needed) -->
[![Final Report](https://img.shields.io/badge/Final%20Report-PDF-red?logo=adobeacrobatreader&logoColor=white)](FinalReport.pdf)
![R](https://img.shields.io/badge/R-4.3.1-276DC3)
![RStudio](https://img.shields.io/badge/RStudio-2023.06.0%2B421-75AADB)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

</div>

---

### 🧠 Overview

This project explores how **single-neuron gene expression profiles**, paired with anatomical and functional properties, can be used to **model and predict neuronal network behavior** in a **neurosymbolic framework**.

**Research Question:**  
How can the transcriptional profiles of individual neurons, characterized by their anatomical and functional properties, be used to model and predict neuronal network behavior in a neurosymbolic framework?

**Dataset:**  
Refine.bio experiment: *Correlating anatomy and function with gene expression in individual neurons...*  
https://www.refine.bio/experiments/SRP094496/correlating-anatomy-and-function-with-gene-expression-in-individual-neurons-by-combining-in-vivo-labeling-patch-clamp-and-single-cell-rna-seq

---

## Quick Start

You’ll need R/RStudio, required packages, and the dataset.

### (1) 🧰 Install Dependencies

Make sure you have:

- **R 4.3.1**
- **RStudio 2023.06.0+421**

### (2) 📦 Install Packages

At minimum, install `BiocManager`:

```r
install.packages("BiocManager")
```

(If your scripts install additional Bioconductor packages, run them after installing `BiocManager`.)

### (3) 📥 Get the Dataset

Download the dataset from Refine.bio and place it in the directory expected by the project (e.g., a `data/` folder in the repo):

```bash
# example (adjust to your repo's expected structure)
mkdir -p data
# download manually from the link in the README and move files into ./data
```

### (4) ▶️ Run the Project

Open the project in **RStudio**, then run the main analysis script(s) (example names below—adjust to your repo):

```r
# Example (adjust to your actual filenames)
source("src/main.R")
# or open an RMarkdown notebook and Knit:
# rmarkdown::render("notebooks/analysis.Rmd")
```

## 🧩 Tech Stack

| Component        | Technology |
| --------------- | ---------- |
| 🧪 Analysis      | R          |
| 🧰 IDE           | RStudio    |
| 🧬 Bio Tools     | BiocManager |



## 👥 Authors

- Ethan Fan — [GitHub](https://github.com/Ethan12103)  
- Matheus Kunzler Maldaner — [GitHub](https://github.com/matheusmaldaner)  
- Nathan Gilman — [GitHub](https://github.com/npgilman)  
- Rama Janco — [GitHub](https://github.com/ramawama)  
