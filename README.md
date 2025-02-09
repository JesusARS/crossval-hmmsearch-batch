# Scripts for Cross-Validation and Batch HMMsearch

This repository contains scripts designed for automating the processing of biological sequences using the HMMER software.  
The `k_fold_cross_validation.py` script allows the creation of Hidden Markov Models (HMM) and performs k-fold cross-validation on a set of target sequences.  
The `batch_hmmsearch.py` script runs an HMM scan in batch mode on a set of `.FASTA` files, isolates the detected sequences, and quantifies them per file.

## System Requirements

- Operating System: Linux  
- External software dependencies:  
  - [hmmer](http://hmmer.org/) (tested on version 3.3.2)  
  - [clustalo](http://www.clustal.org/omega/) (tested on version 1.2.4)  
  - [seqkit](https://bioinf.shenwei.me/seqkit/) (tested on version 2.1.0)  
  - [cd-hit](https://github.com/weizhongli/cdhit/tree/master) (tested on version 4.8.1)  
- Python packages:  
  - [Biopython](https://biopython.org/) (tested on version 1.84)  

**Note**:  
1. Install `hmmer`, `clustalo`, and `seqkit` separately and ensure they are accessible in the system PATH. Also, make sure `Biopython` is installed in your Python environment.  
2. The `HmmDataProcessor` class, contained in this repository, is required for script execution as it provides functions for processing `hmmsearch` results used by both scripts.  

## Scripts

### 1. k-Fold Cross-Validation  

This script performs k-fold cross-validation on biological sequence datasets with different levels of redundancy, which are evaluated as a hyperparameter. A negative dataset for testing is also required.  
For each redundancy level, it generates a folder containing redundancy-reduced data, sequence partitions for k iterations, multiple sequence alignments, HMMs, and `hmmsearch` results for each iteration.  
Additionally, it creates a `results.csv` file that quantifies true positives (TP), false negatives (FN), true negatives (TN), and false positives (FP) for all generated models.  

#### Usage  

```bash
python3 k_fold_cross_validation.py <job_name> -ts <training_data.fasta> -ns <negative_data.fasta> -k <k_folds> -r <redundancy>
```

#### Example  

```bash
python3 k_fold_cross_validation.py GAF_domain -ts GAF_training_data.fasta -ns negative_data.fasta -k 5 -r 100,90,80,70,60
```

#### Parameters  
- job_name: Job name.  
- -ts, --training_data_file: Training data file in FASTA format.  
- -ns, --negative_set_file: Negative set file in FASTA format.  
- -k, --k_folds: Number of partitions for cross-validation (default: 5).  
- -r, --redundancy: Redundancy levels to evaluate as a hyperparameter, separated by commas (default: 100).  

### 2. Batch HMM Search  

This script runs `hmmsearch` with a given HMM on a batch of FASTA files specified in a CSV file.  
Each analyzed file generates a folder containing the `hmmsearch` result and a FASTA file with sequences found below the specified E-value threshold.  
Additionally, it creates a `quantification.csv` file with the number of detected proteins per FASTA file.  

#### Usage  

```bash
python3 batch_hmmsearch.py <job_name> <HMM.hmm> <FASTA_files.csv> <E_value>
```

#### Example  

```bash
python3 batch_hmmsearch.py GAF HMM.hmm FASTA_files.csv 1e-5
```

#### Parameters  
- job_name: Job name.  
- hmm_path: Path to the HMM file.  
- FASTA_files: CSV file containing paths to FASTA files.
- E_value: E-value threshold for searches.

