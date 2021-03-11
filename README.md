# CMU-02718-HW5
Fall 2020 - Computational Medicine course project - SARS-CoV-2 Subunit Vaccine Design

### Summary
The goal of this project is to design a subunit vaccine targeting the spike glycoprotein of SARS-CoV-2. Specifically, this project focuses on selecting a set of peptides which could theoretically be delivered using an RNA vaccine.

Work involved cleaning, standardizing, and encoding MHC binding data, training regression models (Random Forest) to predict pIC50 values for a given HLA allele and peptide k-mer, and designing a novel approach for taking binding predictions from SARS-CoV-2 spike glycoprotein (`.fasta` format) and selecting a subset of peptides to maximize coverage of HLA allele reference sets.

Analysis was performed using Jupyter Notebook and Python.



### Project Structure
Vaccine design functions and helper functions are implemented as Python package `vaccine` under /packages.
