## Install
Clone this directory.
`git clone git@github.com:bihealth/atlatl.git`

Create a fresh conda environment.
`conda create -n atlatl --file envs/atlatl.yml -y`

Then activate the environment.
`conda activate atlatl`

cd into it.
`cd atlatl`

Then install atlatl via pip. The installation will be kept in this directory.
`pip install -e .`

Done, now you can test and use it.

## Run Tests of atlatl
#### Input Preprocessing
`atlatl workflow_prep --workdir tests/data/ --alleles tests/data/input/construct_a.bed tests/data/input/construct_b.bed --sequences tests/data/input/sequences.fasta --reference tests/data/input/ref.fasta --add_annotations tests/data/input/annotations.bed --cores 2`

#### Read Mapping
Once the input is preprocessed, a workflow in atlatl can be used to map the reads to the preprocessed reference.
`atlatl workflow_map --workdir tests/data/ --constructs tests/data/prep/constructs.fasta --in_reads tests/data/input/reads.fastq --reference tests/data/prep/ref_aug.fasta --annotations tests/data/prep/annotations.bed --cores 2`

#### Breakend Calculation
`atlatl breakends --alignment tests/data/align/alignment.bam --bedfile tests/data/prep/annotations.bed --out_dir tests/data/breakends --min_breakends 5`

#### Assebmle and Visualize Haplotypes
`atlatl group_assemble_and_visualize --breakends tests/data/breakends/construct_a.tsv --annotations tests/data/prep/annotations.bed --reads tests/data/input/reads.fastq --reference tests/data/prep/ref_aug.fasta --outdir tests/data/assembly --threads 4`

## Adjustment for Your Usage
You can have a look at the files in `tests/data/input/` and take these as an example to work with your own data. The aim of the tests is to build two constructs `construct_a` and `construct_b` which are defined by both .bed files with their exact name. The construct_x.bed files only define the intervals of source sequences from which the constructs are built. The sequences adressed in the construct_x.bed files need to be located in any of the given .fasta files (`ref.fasta` and `sequences.fasta`). The `annotations.bed` file add annotations to the sequences. You can and you should for example define the locus on the reference genome where the constructs will be inserted. This is necessary for atlatl to calculate the correct on-target efficiency.

Once you have defined the input files correctly, you can run all downstream steps easily without any additional data manipulation.

![Continuous Integration Status](https://github.com/bihealth/clear-CNV/workflows/CI/badge.svg)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/2eaafb57fbb74a46b918e9f58142c880)](https://www.codacy.com?utm_source=github.com&utm_medium=referral&utm_content=bihealth/clear-CNV&utm_campaign=Badge_Coverage)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/2eaafb57fbb74a46b918e9f58142c880)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=bihealth/clear-CNV&amp;utm_campaign=Badge_Grade)
[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
