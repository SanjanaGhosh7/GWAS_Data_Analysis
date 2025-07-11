# 🧬 GWAS_Data_Analysis


**Author**: Sanjana Ghosh  
**Role**: Bioinformatics Postgraduate | Aspiring Genomics Data Scientist  
**Tools**: Python, R, Pandas, NumPy, ggplot2, Jupyter, RStudio

---

## 📌 Project Overview

This project explores genome-wide association data from the **GWAS Catalog** to identify genes significantly associated with human traits and diseases. Conducted using Python and R, this analysis focuses on frequency calculations, gene–trait associations, SNP mapping, population bias identification, odds ratios, and data visualization.

---

## 📂 Data Source

- Original full dataset:  
  🔗 https://www.ebi.ac.uk/gwas/api/search/downloads/full

---

📌 NOTE: `GWAS_Analysis.ipynb` is the original step-by-step analysis of the GWAS Catalog dataset with all the comments and outputs, ideal for learning, debugging, and understanding the full flow interactively. And `GWAS_Analysis.py` is  the converted and cleaned Python script version of the notebook. It contains the complete GWAS data analysis pipeline for reproducible command-line execution.
Best suited for automation, integration into larger workflows, or running as a script.

## 🔍 Analysis Performed

1. **Gene Frequency Analysis** – Identified top 10 genes with the highest association counts.
2. **Trait Diversity** – Counted unique traits and mapped gene–trait relationships.
3. **SNP & Chromosome Mapping** – Extracted top genes’ SNP-Risk Alleles and chromosome IDs.
4. **Population Insight** – Identified the most common population associated with each top gene.
5. **Odds Ratio Filtering** – Identified top genes with the strongest effect sizes.
6. **Data Visualization** – Plotted lollipop charts, stacked bars, and OR-based bar plots using `ggplot2`.
   
---

## 📊 Visualizations

- Lollipop plot: Top genes by frequency across traits  
- Stacked bar chart: SNP–chromosome mapping  
- Bar plot: Genes with highest odds ratios  

Plots are saved in the `Outputs/plots/` folder.

---

## 👩‍💻 Skills Demonstrated

- Genomic data wrangling with Python  
- Statistical filtering and SNP analysis  
- Data visualization using `ggplot2` in R  
- Handling multi-gene entries and real-world large bioinformatics datasets  
- Scientific documentation and reproducible research 

---

## 🚀 Future Plans

- Extend this project to include enrichment analysis using `Reactome` or `KEGG`
- Integrate other omics layers from TCGA for multi-omics insights
- Apply ML/AI to predict trait-gene linkages from variant patterns

---

## 📬 Feel free to Connect or Collaborate

**Sanjana Ghosh**  
[LinkedIn](www.linkedin.com/in/sanjana-ghosh-2a5b7c11d) | sanjanaghosh150@gmail.com

