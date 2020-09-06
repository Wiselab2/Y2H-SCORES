# Instructions to run the toy dataset for Y2H-SCORES

In this repository you will find the instructions to run the toy dataset of the software Y2H-SCORES published in  Velásquez-Zapata *et.al.* (2020) Y2H-SCORES: A statistical framework to infer protein-protein interactions from next-generation yeast-two-hybrid sequencing data.   

**File description**

 This repository should contain the following files:

* bait1\_input\_arguments.csv, bait2\_input\_arguments.csv: minimum configuration files to run Y2H-SCORES. You can use the same configuration files as the ones used to run the NGPINT pipeline, keeping the same file structure. The full paths to the input files and the sample names that are written in this file are used to run Y2H-SCORES.
* Two folders called bait1 and bait2 with the corresponding bait#\_final\_report.csv and bait#\_salmon\_counts.mat obtained as outputs from the NGPINT pipeline. These files are inputs to calculate the Y2H-SCORES.
* Text file with the full paths to the configuration files that contain the input arguments that were used to run NGPINT. This file should be indicated with the argument -fofn. For this example the corresponding file was called fofn\_for\_compute\_scores.txt. This file provides the path information to *compute\_scores_newInt.py* so it can read the configuration files.
**Running Y2H-SCORES with the toy example dataset**

To run Y2H-SCORES with the toy dataset in Linux/Mac use the terminal and in Windows the Anaconda Prompt :

```
cd Y2H-SCORES
python compute_scores_newInt.py -fofn toy_dataset/fofn_for_compute_scores.txt -out_dir output_toy_dataset/  -threshold_p_val 0.01  -threshold_fold_change 2  -threshold_enrichment_score 0.5
```
In this example we are running the Y2H-SCORES with the arguments:

```
-fofn toy_dataset/fofn_for_compute_scores.txt contains the text file with the full paths to the configuration files bait1\_input\_arguments.csv and bait2\_input\_arguments.csv.

-out_dir output_toy_dataset/ corresponds to the name of the folder that will be generated to save the outputs from Y2H-SCORES. The software will create this folder, it is not necessary to make it beforehand.

-threshold_p_val 0.01 is the desired threshold for the p-values that are used to calculate the specificity score. 

-threshold_fold_change 2 the desired threshold for the fold-changes that are used to calculate the specificity score. Larger values indicate more stringent scores.

-threshold_enrichment_score 0.5 the desired threshold for the p-values that are used to calculate the enrichment score. Smaller values indicate more stringent scores.
```
After running the Y2H-SCORES software there will be a new folder in *toy\_example*  called *output\_toy_dataset*. This folder contains the *Total\_scores.csv* file with the compilation of the three scores and the columns: prey, bait, Enrichment score, Specificity score, *In-frame* score, In-frame Transcripts, sum of scores and borda scores as an ensemble. To prioritize interactions use the ensemble by borda counts.