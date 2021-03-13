# CMU-02718-HW5
Fall 2020 - Computational Medicine course project - SARS-CoV-2 Subunit Vaccine Design

### Summary
The goal of this project is to design a subunit vaccine targeting the spike glycoprotein of SARS-CoV-2. Specifically, this project focuses on selecting a set of peptides which could theoretically be delivered using an RNA vaccine.

Work involved cleaning, standardizing, and encoding MHC binding data, training regression models (Random Forest) to predict pIC50 values for a given HLA allele and peptide k-mer, and designing a novel approach for taking binding predictions from SARS-CoV-2 spike glycoprotein (`.fasta` format) and selecting a subset of peptides to maximize coverage of HLA allele reference sets.

Analysis was performed using Jupyter Notebook and Python.



### Project Structure
Vaccine design functions and helper functions are implemented as Python package `vaccine` under /packages.


### Algorithm
#### Binding predictions

- For each MHC class, train separate regression models (Random Forest) for each HLA allele in the reference set. MHC I data is limited to peptides of length 9; MHC II data is limited to peptides of length 15.

- For each MHC class, subdivide the spike sequence into overlapping k-mers of matching length (9 or 15).
Use the regression models to predict pIC50 binding values for each (k-mer, allele) combination.

- Create data structure P to track the set of alleles covered by each peptide k-mer. Sort P by number of alleles covered, descending.


#### Minimum Coverage
- Generate a peptide set which covers as many alleles in a single MHC class at least once. Order alleles so that the least-covered allele is covered first. Order peptides so that the peptide with the possible allele coverage is selected first.

```
for design limit D, set of alleles A, set of peptides P, and pIC50 threshold p_min:

  limit P to peptides with pIC50 > p_min

  U: set of uncovered alleles from A covered by at least one peptide in P
  
  S: set of selected peptides


  while |S| ≤ D and |U| > 0 and |P| > 0:

    select allele a which is covered by fewest number of peptides in P
    select first p in P which covers a
    
    add p to S
    remove peptides from P which overlap with p
    remove from U each allele covered by p
    remove from U any alleles with nonzero coverage in revised P
    
  return S
```
