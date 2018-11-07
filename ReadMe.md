# About this project

This repository **meetU-MasterStudents/2018---2019-Equipe-6** takes part in the [Meet-U](http://meet-u.org/edition_2018.html) project. The aim of the 2018's edition is to predict the 3D structure of a protein taken as an entry the query. Our team is focused on the first part of the work, therefore the goal of this repository is to find the best profile-profile comparison between template sequences and a target sequence by transforming them into profiles and then optimizing a profile-profile alignment.

We choosed to use C++ and Python.

# Team

Our team is composed of 5 students from UPMC university. Three of us compose the biological subteam and 2 of us compose the development subteam.

- [Yasser Mohseni Behbahani](https://github.com/yassermb)
- [Bénédicte Colnet](https://github.com/BenedicteColnet)
- [Gabriela Lobinska](https://github.com/gabriela3001)
- [Irene Mauricette Mendy](https://github.com/reinamauricette) 
- [Amandine Sandri](https://github.com/amandinesandri) 


# Algorithms and Packages used

## 1) Basic Local ALignment Search Tool (BLAST) algorithm is used to search sequence databases for optimal alignement to the [query](https://github.com/meetU-MasterStudents/2018---2019-Equipe-6/blob/master/query.fasta)

**Installation and main instructions** can be used at the following [link](https://www.ncbi.nlm.nih.gov/books/NBK52640/)
**Useful Instruction** can be used at the following [link](https://www.ncbi.nlm.nih.gov/books/NBK52640/)

## 2) MUltiple Sequene Comparison by Log-Expectation MUSCLE that creates multiple Alignments of our sequence database previously found. 
PARAMETRES DE MUSCLE ?? 

**Instructions** can be used at the following [link](https://petrov.stanford.edu/software/src/muscle3.6/muscle3.6.html)
**Installation** is done following below command line:

```bash
sudo apt-get install muscle
```

## 3) Command lines to build the package
### Build the packages
**Packages** are quite simple to use. 
Start by writting the below command line to build the package:
```bash
sudo g++ -o profile main.cpp Profile.cpp -std=c++11
```

### Execute the packages
Next, use the below command line to execute the packages: 
```bash
./profile 
```

When alignment has finished, ``profile`` will display the matrix of frequency of each amino acid in the query and template sequences.

# Useful links 

[Pairwise Sequence Alignment PSA](https://www.ebi.ac.uk/Tools/psa/)

Multiple Sequence Alignment MSA](https://www.ebi.ac.uk/Tools/msa/clustalo/)

[To download PDB FTP](http://www.wwpdb.org/download/downloads)
[To download formatted data records for a list of input UIDs](https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch)

# References

 1. GUOLI WANG, ROLAND L. and DUNBRACK JR. "Scoring profile-to-profile sequence alignments" Prot Sci. 2004 Jun; 13(6) 1612-1626.
 2. 

# Authors
```
This file is part of Equipe 6.

```
