# Sequence-data-analysis-for-noise-control
The source code for the paper "Quantitative Control of Noise in Mammalian Gene Expression by Dynamic Histone Regulations".

**If you want to run the code, please read the following descrptions and pay attention to the file path. The file path in the code should be modified according to your file path.**

## build_reference
It includes the code to build the HeLa-AB1 reference genome which for ATAC-seq, ChIP-seq and scATAC-seq data processing. The metadata directory includes some input and intermediate files which are used in the code.

### Procedure
step1: run the code in `get_reference.sh`

## ChIP-seq
It includes the code to process the ChIP-seq data. 

### Procedure
step1: run the code in `chip_analysis.sh`

step2: run the code in `calculate_depth.py`

## ATAC-seq
It includes the code to process the ATAC-seq data. 

### Procedure
step1: run the code in `atac_analysis.sh`

step2: run the code in `calculate_depth.py`

## scATAC-seq
It includes the code to process the scATAC-seq data. The metadata directory includes some input and intermediate files which are used in the code.

### Procedure
step1: run the code in `scatac_analysis.sh`

step2: run the code in `calculate_depth.py`

## scRNA-seq
It includes the code to process the scRNA-seq data. The metadata directory includes some input and intermediate files which are used in the code.

### Procedure
step1: run the code in `new_ref.sh`

step2: run the code in `10x_scRNA.sh`

step3: run the code in `qc.R`

step4: run the code in `combine_repeat.py`

step5: run the code in `F9_scRNA_analysis.m`, `Hela_scRNA_analysis.m` and `run_saver.R`
