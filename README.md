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

#### Balanced Coverage
- Minimum coverage won't result in equal coverage of every allele in a given MHC class. With remaining design space, choose peptides to equalize coverage as much as possible. Order alleles so least-covered is improved first (non-zero only). 


```
for design limit D, set of selected peptides S, and pIC50 value p_min:

	limit P to peptides with pIC50 > p_min and which are not in S

	C: set of alleles covered by S

	L: set of alleles for which no additional coverage can be found in P, initially empty

	
	while |S| ≤ D and |P| > 0 and |L| ≠ |C|:
		
		select allele a in C with the least coverage which is not in L
		search for peptide p in P which covers a
		if p exists:
			add p to S
			remove peptides from P which overlap with p
		else:
			add a to L

	return S
```

Note: Implementation attempts to increase coverage for each allele in C once, even if that allele remains the least-covered in C and peptides remain in P which could cover that allele. Implementation should be revised to match pseudocode.

#### Final Design
- A final set of peptides must be selected from the design sets which have been generated. 
- To balance coverage of MHC classes, half of the peptides in the final design will be chosen from the design sets of each MHC class. For each MHC class, selection considers all design sets generated from each combination of (D,p_min). The goal is to select peptides with the highest pIC50 values to maximize coverage.


```
Fc: final design set for a single MHC class

Df: user-defined size of final design set

p_current = max{p_min for design sets}


while |Fc| < |Df|:

	Smin : set of min coverage peptides from design set using p_current as p_min

	Sbal : set of peptides added to balanced design set

	
	if Fc does not contain all peptides in Smin:

		add peptide in Smin not in Fc to Fc


	else if Fc does not contain all peptides in Sbal:

		add peptide in Sbal not in Fc to Fc
	
	else:

		set p_current to next highest p_min in design sets
	
return Fc
```

Note: When choosing from Smin or Sbal, implementation favors peptides which covers the most alleles.

