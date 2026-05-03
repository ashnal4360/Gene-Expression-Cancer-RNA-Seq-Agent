# Gene-Expression-Cancer-RNA-Seq-Agent

## Project Overview
A Streamlit-based gene expression analysis app that compares two cancer classes from RNA-seq data, ranks the most differentially expressed genes, and visualizes the results with a volcano-style plot.

## Why This Project
This project demonstrates bioinformatics analysis, Python data wrangling, statistical testing, and interactive data visualization in one reproducible workflow.

## Dataset
Gene Expression Cancer RNA-Seq dataset from Kaggle, containing gene-expression features and cancer class labels. 
The two datasets can be downloaded from here: https://www.kaggle.com/datasets/waalbannyantudre/gene-expression-cancer-rna-seq-donated-on-682016

## What It Does
- Loads large RNA-seq CSV files.
- Lets the user choose two cancer classes.
- Runs two-sample t-tests across gene features.
- Ranks the strongest signals by p-value.
- Generates a downloadable results table and plot.

## Key Results
The app identifies genes that differ strongly between cancer groups, making it a useful exploratory screening tool.

## Tech Stack
Python, Streamlit, Pandas, NumPy, SciPy, Plotly.

## How to Run
1. Install dependencies.
2. Put `data.csv` and the labels file in the project folder.
3. Run `streamlit run app.py`.

## Next Improvements
- Add multiple-testing correction.
- Add normalization options.
- Add gene annotation and pathway enrichment.
