# Instructions to run Y2H-SCORES

This repository contains the software Y2H-SCORES published in  Velásquez-Zapata *et.al.* (2021) Next-generation yeast-two-hybrid analysis with Y2H-SCORES identifies novel interactors of the MLA immune receptor.   

**Overview**

This manual describes how to operate Y2H-SCORES to prioritize candidate protein-protein interactors reported by NGPINT (Banerjee et. al., 2020). Y2H-SCORES takes the output files from NGPINT (salmon counts and fusion counts) and computes three scores to rank protein-protein interactor candidates in a Y2H Next Generation Interaction Screening (NGIS) experiment.  

Y2H-SCORES is an aggregate of three experimental outcomes: 1) the enrichment score measures the non-selected population as a baseline to detect which preys are significantly enriched under selection, 2) the specificity score uses selected samples as a control baseline to measure the specificity of a prey, and 3) the *in-frame* score makes use of the fusion read information to identify *in-frame* enrichment of the prey fragments under selection. 

Y2H-SCORES has been published in PLoS Computational Biology, it can be accessed from https://doi.org/10.1371/journal.pcbi.1008890. Please refer to the paper for more details about each module. The entire pipeline has been coded in R.

**Citations**

If you use Y2H-SCORES, please cite  

* Valeria Velásquez-Zapata, J Mitch Elmore, Sagnik Banerjee, Karin S Dorman, Roger P Wise (2021) Next-generation yeast-two-hybrid analysis with Y2H-SCORES identifies novel interactors of the MLA immune receptor. PLoS Comput Biol 17(4): e1008890. https://doi.org/10.1371/journal.pcbi.1008890

* Sagnik Banerjee, Valeria Velásquez-Zapata, Gregory Fuerst, J Mitch Elmore, Roger P Wise (2020) NGPINT: a next-generation protein–protein interaction software. Briefings in Bioinformatics: bbaa351, https://doi.org/10.1093/bib/bbaa351

**Software requirements**

* R
* R packages: reshape2, tidyverse, psych, DESeq2, mass, optparse

**Downloading Y2H-SCORES**

For Linux/Mac operating systems, in the terminal: 

```
git clone https://github.com/vvelasqz/Y2H-SCORES.git
```

In Windows download the folder from https://github.com/Wiselab2/Y2H-SCORES/tree/master/Software.

**Running Y2H-SCORES**

Y2H-SCORES is coupled with the software NGPINT so it is designed to work with the outputs from this program. To run Y2H-SCORES, use the run\_scores.R, a script that process the input count tables and generate the three scores. It requires the following files:

1. Text file with the full paths to the configuration files that contain the input arguments used to run NGPINT. This file should be indicated with the argument --fofn.

2. Outputs from the pipeline NGPINT. Y2H-SCORES uses the BaitName\_final\_report.csv and BaitName\_salmon\_count.matrix files. These files contain the total counts and the fusion counts that will be used to calculate the Y2H-SCORES. Each configuration file in --fofn should contain an output\_directory cell with the location of the outputs of NGPINT. 

3. Configuration files from the NGPINT software.  They should contain all the paths to the input files and the sample names.

4. The specificity score is calculated by separating the baits in groups of 10. The user can provide one text file with the baits in each group separated by comma. and indicated with the argument --spec_groups. If no file is provided, the baits will be grouped randomly.

  To run Y2H-SCORES in Linux/Mac use the terminal and in Windows the Anaconda Prompt: 

```
cd Y2H-SCORES/Software
Rscript run_scores.R --fofn --out_dir --spec_groups --spec_p_val --spec_fold_change --enrich_p_val --enrich_fold_change --normalized 
```

With the arguments:

```
--fofn: text file with the full paths to the configuration files

--out_dir: full path to the output directory to save the calculations from the Y2H-SCORES.  
          This should be different from the output\_directory in the configuration files.
          
--spec_groups: text file with baits to be analyzed together to calculate the specificity score. Default NULL 
			   Each line should contain each group, separating baits by comma. 
			   If no file is provided, the baits will be grouped randomly.

--spec_p_val: threshold for the p-values that are used to calculate the specificity score.
              Values between [0,1]. Default value of 1. Smaller values indicate more stringent scores.

--spec_fold_change: threshold for the fold-changes that are used to calculate the specificity score.  
                    This should be larger or equal to zero. Default of zero. 
                    Larger values indicate more stringent scores.

-enrich_p_val: the desired threshold for the p-values that are used to calculate the enrichment score.  
               Values between [0,1]. Default value of 1. Smaller values indicate more stringent scores.
               
--enrich_fold_change: threshold for the fold-changes that are used to calculate the enrichment score.  
                      This should be larger or equal to zero. Default of zero. 
                      Larger values indicate more stringent scores.
                    
--normalized: Boolean value (T or F) indicating if the counts are normalized. 
			  Default false (F), then the program will implement library size normalizaiton
```

**Output**

After running Y2H-SCORES there should be an output directory named as indicated with the --out\_dir argument. This folder should have a file called Total\_scores.csv with the compilation of the three scores for every prey/bait combination.  The table should have the following columns: prey, bait, Enrichment score, Specificity score, *In-frame* score, *In-frame* prey transcripts (list of prey transcripts that have the maximum *in-frame* score for each prey), Sum scores and Borda score.

To prioritize interactors use the ensemble scores given by the Borda score column. Alternatively, each of the three scores and their sum can be used as a guide. Higher values indicate more likely interactors. 

## Toy example

 The *toy\_example* directory contain a minimum dataset that is required to run Y2H-SCORES. 

**File description**

 This repository should contain the following files:

* bait1\_input\_arguments.csv, bait2\_input\_arguments.csv: minimum configuration files to run Y2H-SCORES. You can use the same configuration files as the ones used to run the NGPINT pipeline, keeping the same file structure. The full paths to the input files and the sample names that are written in this file are used to run Y2H-SCORES.

* Two folders called bait1 and bait2 with the corresponding bait#\_final\_report.csv and bait#\_salmon\_counts.mat obtained as outputs from the NGPINT pipeline. These files are inputs to calculate the Y2H-SCORES.

* Text file with the full paths to the configuration files that contain the input arguments that were used to run NGPINT. This file should be indicated with the argument -fofn. For this example the corresponding file was called fofn\_for\_compute\_scores.txt. This file provides the path information to *compute\_scores_newInt.py* so it can read the configuration files.

  **Running Y2H-SCORES with the toy example dataset**

To run Y2H-SCORES with the toy dataset in Linux/Mac use the terminal and in Windows the Anaconda Prompt :

```
cd Y2H-SCORES/Software
Rscript run_scores.R --fofn "toy_dataset/fofn_for_compute_scores.txt" --out_dir "output_toy_dataset/" --spec_p_val 0.5 --spec_fold_change 2 --enrich_p_val 0.5 --enrich_fold_change 0 --normalized F
```

In this example we are running the Y2H-SCORES with the arguments:

```
--fofn toy_dataset/fofn_for_compute_scores.txt contains the text file with the full paths to the configuration files 
  bait1\_input\_arguments.csv and bait2\_input\_arguments.csv.

--out_dir output_toy_dataset/ corresponds to the name of the folder that will be generated to save the outputs from Y2H-SCORES. The software will create this folder, it is not necessary to make it beforehand.

--spec_p_val 0.5 is the desired threshold for the p-values that are used to calculate the specificity score. 

--spec_fold_change 2 the desired threshold for the fold-changes that are used to calculate the specificity score. Larger values indicate more stringent scores.

--enrich_p_val 0.5 the desired threshold for the p-values that are used to calculate the enrichment score.  
  Smaller values indicate more stringent scores.
  
--enrich_fold_change 0 the desired threshold for the fold-changes that are used to calculate the enrichment score. Larger values indicate more stringent scores.

--normalized F This ask the program to normalize the raw counts using library size method.
  
```

After running the Y2H-SCORES software there will be a new folder in *toy\_example*  called *output\_toy_dataset*. This folder contains the *Total\_scores.csv* file with the compilation of the three scores and the columns: prey, bait, Enrichment score, Specificity score, *In-frame* score, In-frame prey transcripts, Sum of scores and Borda scores.

---

## Running Y2H-SCORES with Containers (Docker & Apptainer / SIF)

Container images with R and all required CRAN/Bioconductor dependencies pre-installed are automatically built and published to the GitHub Container Registry (GHCR). The container entrypoint runs `run_scores` directly (`Rscript /opt/Y2H-SCORES/run_scores.R`), so you only need to pass the `run_scores.R` command-line flags.

> **Important path notes when using containers:**
> * **Trailing slash on `--out_dir`:** Always include a trailing `/` on `--out_dir` (e.g., `--out_dir "output_toy_dataset/"`), because `run_scores.R` concatenates output filenames directly onto the `--out_dir` string.
> * **Paths inside `--fofn` and bait config CSVs:** All paths listed in your `--fofn` file and in the `output_directory` row of each `*_input_arguments.csv` file must be accessible inside the container (relative to your working directory or inside a mounted path).

---

### 1. Using Docker (Pre-built Image from GHCR)

#### Pull the image
```bash
docker pull ghcr.io/vvelasqz/y2h-scores:latest
```

#### Check CLI options
```bash
docker run --rm ghcr.io/vvelasqz/y2h-scores:latest --help
```

#### Run on the toy dataset
From the `Software/` directory, mount the current directory into `/data` and set `/data` as the container working directory (`-w /data`) so the relative paths inside `toy_dataset/fofn_for_compute_scores.txt` resolve properly:

```bash
cd Y2H-SCORES/Software

docker run --rm \
  -v "$(pwd):/data" \
  -w /data \
  ghcr.io/vvelasqz/y2h-scores:latest \
  --fofn "toy_dataset/fofn_for_compute_scores.txt" \
  --out_dir "output_toy_dataset/" \
  --spec_p_val 0.5 \
  --spec_fold_change 2 \
  --enrich_p_val 0.5 \
  --enrich_fold_change 0 \
  --normalized F
```

#### Run on your own data
If your configuration and NGPINT output files use absolute host paths (e.g., under `/path/to/project`), mount that directory at the same path inside the container (or use paths relative to `/data` with `-w /data`):

```bash
docker run --rm \
  -v "/path/to/project:/path/to/project" \
  -w "/path/to/project" \
  ghcr.io/vvelasqz/y2h-scores:latest \
  --fofn "/path/to/project/fofn_for_compute_scores.txt" \
  --out_dir "/path/to/project/scores_output/" \
  --spec_p_val 0.05 \
  --spec_fold_change 1 \
  --enrich_p_val 0.05 \
  --enrich_fold_change 0 \
  --normalized F
```

---

### 2. Building from the `Dockerfile` Locally

If you want to build the Docker image locally from `Software/Dockerfile`:

```bash
# From the repository root:
docker build -t y2h-scores:local ./Software

# Or from inside Software/:
cd Y2H-SCORES/Software
docker build -t y2h-scores:local .
```

Run the locally built image the same way:

```bash
cd Y2H-SCORES/Software

docker run --rm \
  -v "$(pwd):/data" \
  -w /data \
  y2h-scores:local \
  --fofn "toy_dataset/fofn_for_compute_scores.txt" \
  --out_dir "output_toy_dataset/" \
  --spec_p_val 0.5 \
  --spec_fold_change 2 \
  --enrich_p_val 0.5 \
  --enrich_fold_change 0 \
  --normalized F
```

---

### 3. Obtaining the Apptainer / Singularity `.sif` Container

You can either pull the pre-built native `.sif` image directly from GHCR via ORAS (fastest) or build a `.sif` locally from the Docker image.

#### Option A: Pull the pre-built `.sif` via ORAS (Recommended)
```bash
apptainer pull y2h_scores.sif oras://ghcr.io/vvelasqz/y2h-scores-sif:latest
```

#### Option B: Build a `.sif` from the GHCR Docker image
```bash
apptainer build y2h_scores.sif docker://ghcr.io/vvelasqz/y2h-scores:latest
```

#### Option C: Build a `.sif` from a local Docker daemon image
If you built `y2h-scores:local` locally with Docker:
```bash
apptainer build y2h_scores.sif docker-daemon://y2h-scores:local
```

---

### 4. Executing Y2H-SCORES with Apptainer and the `.sif` File

By default, Apptainer automatically bind-mounts and preserves your current working directory on the host, while `run_scores.R` resolves its internal helper scripts from `/opt/Y2H-SCORES` inside the SIF.

#### Check CLI options
```bash
apptainer run y2h_scores.sif --help
```

#### Run the `.sif` on the toy dataset
```bash
cd Y2H-SCORES/Software

apptainer run y2h_scores.sif \
  --fofn "toy_dataset/fofn_for_compute_scores.txt" \
  --out_dir "output_toy_dataset/" \
  --spec_p_val 0.5 \
  --spec_fold_change 2 \
  --enrich_p_val 0.5 \
  --enrich_fold_change 0 \
  --normalized F
```

*(You can also execute the `.sif` file directly: `./y2h_scores.sif --fofn ...`)*

#### Run the `.sif` on HPC / external directories
If your NGPINT outputs and configuration files live on a shared cluster filesystem outside your home/current directory (such as `/scratch` or `/work`), pass `--bind` (`-B`) so Apptainer mounts that path inside the container:

```bash
apptainer run \
  --bind /scratch/my_y2h_project:/scratch/my_y2h_project \
  y2h_scores.sif \
  --fofn "/scratch/my_y2h_project/fofn_for_compute_scores.txt" \
  --out_dir "/scratch/my_y2h_project/scores_output/" \
  --spec_p_val 0.05 \
  --spec_fold_change 1 \
  --enrich_p_val 0.05 \
  --enrich_fold_change 0 \
  --normalized F
```

#### Run directly without pre-downloading a `.sif` file
Apptainer can also execute the remote ORAS or Docker image directly (caching it in `~/.apptainer/cache`):

```bash
apptainer run oras://ghcr.io/vvelasqz/y2h-scores-sif:latest --help
```



