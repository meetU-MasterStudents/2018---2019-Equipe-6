# [EQUIPE AMONT] Some words about strategy

Our team proposes a pipeline that 
- First uses PSI-BLAST on the query to catch the best hits over the data based asked bu the user.
- Second, we uses MUSCLE (MUltiple Sequence Comparison by Log- Expectation) for creating multiple alignments of proteins sequences using BLAST result.
- Third, we build a profile from all the HOMSTRAD data base and a profile from the query you want to test. This profil is a PSSM profile. Then, we compare each profile of the data with the dot product. This gives a score for each comparison.

A unique function can launch directly the process over the benchmark list and gives back accuracy plot (variation on parameter on database, blast, and so on is possible). the accuracy plot takes into account the Fold and Super Famille prediction. The launch is multi-thread to easily scan all data.

The aval team will have output files and PSSM profiles available.



# About this project

This repository **meetU-MasterStudents/2018---2019-Equipe-6** takes part in the [Meet-U](http://www.meet-u.org/) project. The aim of the 2018's edition is to predict the 3D structure of a protein taken as an entry the query. Our team is focused on the first part of the work, therefore the goal of this repository is to find the best profile-profile comparison between template sequences and a target sequence by transforming them into profiles and then optimizing a profile-profile alignment.

We chose to use C++ and Python.

# Team

Our team is composed of 5 students from UPMC university. 

- [Yasser Mohseni Behbahani](https://github.com/yassermb)
- [Bénédicte Colnet](https://github.com/BenedicteColnet)
- [Gabriela Lobinska](https://github.com/gabriela3001)
- [Irene Mauricette Mendy](https://github.com/reinamauricette) 
- [Amandine Sandri](https://github.com/amandinesandri) 

# Flowchart

![Flowchart](https://docs.google.com/drawings/d/1QtJI3bWGgn3PLK5nTrLXih2F6LsaX6BcbD77oNbhQek/export/png)

# Algorithms and Packages used

### 1) BLAST+ (version 2.7.1+)
Basic Local ALignment Search Tool ([BLAST](https://en.wikipedia.org/wiki/BLAST)) algorithm is used to search sequence databases for optimal alignement to the [query](https://github.com/meetU-MasterStudents/2018---2019-Equipe-6/blob/master/query.fasta)

**Installation and main instructions** can be used at the following [link](https://www.ncbi.nlm.nih.gov/books/NBK52640/)

If the local database is used for BLAST algorithms please first build the local database using following command (This is an example to build the local swissprot database):
```bash
makeblastdb -in uniprot_sprot.fasta -parse_seqids -dbtype prot
```

### 2) MUSCLE
Multiple Sequene Comparison by Log-Expectation (MUSCLE) that creates multiple Alignments using the output of BLAST.

**Instructions** can be used at the following [link](https://petrov.stanford.edu/software/src/muscle3.6/muscle3.6.html)
**Installation** is done following below command line:

```bash
sudo apt-get install muscle
```

### 3) Run

Here are the command lines to build the package. 

**Compilation of cpp files** Only needed once. 
Start by writting the below command line to build the package:
```bash
sudo g++ -o Profile main.cpp Profile.cpp -std=c++11
```

**Execute the packages**
Next, use the below command line to execute the packages: 
```bash
Benchmark.py
```
These command launches the analysis on the benchmark data. Therefore you should upload the test-dataset from the common folder of Meet-U and gives the absolute path to the Benchmark.py.
Then differents parameters are possible.

# Examples

Here below are examples of output files you should get after running the previous commands:
- [Example of Query Sequence](https://github.com/meetU-MasterStudents/2018---2019-Equipe-6/edit/master/Example/Query.?)

- [Exmple of MSA from Blast](https://github.com/meetU-MasterStudents/2018---2019-Equipe-6/edit/master/Example/MSA.blast_out.?)

- [Example of Query Profiles](https://github.com/meetU-MasterStudents/2018---2019-Equipe-6/edit/master/Example/query.aamtx.?)

- [Example of Homstrad Profiles](https://github.com/meetU-MasterStudents/2018---2019-Equipe-6/edit/master/Example/template.aamtx.?)

- [Example of scores from profiles-profiles alignments](https://github.com/meetU-MasterStudents/2018---2019-Equipe-6/edit/master/Example/score.foldrec.?)


# Useful links 

[Pairwise Sequence Alignment PSA](https://www.ebi.ac.uk/Tools/psa/)

[Multiple Sequence Alignment MSA](https://www.ebi.ac.uk/Tools/msa/clustalo/)

[To download PDB FTP](http://www.wwpdb.org/download/downloads)

[To download Uniprot FTP](http://BLABLABLA on va trouver)

[To download formatted data records for a list of input UIDs](https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch)

# References

 1. GUOLI WANG, ROLAND L. and DUNBRACK JR. "Scoring profile-to-profile sequence alignments" Prot Sci. 2004 Jun; 13(6) 1612-1626.
  

# Authors
```
This file is part of Equipe 6.

```
