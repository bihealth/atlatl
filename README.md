`atlatl workflow_prep --workdir tests/data/ --alleles tests/data/input/construct_a.bed tests/data/input/construct_b.bed --fastas tests/data/input/sequences.fasta --reference tests/data/input/ref.fasta --add_annotations tests/data/input/annotations.bed --cores 2`

`atlatl workflow_map_ont --workdir tests/data/ --constructs tests/data/prep/constructs.fasta --in_reads tests/data/input/reads.fastq --reference tests/data/prep/ref_aug.fasta --annotations tests/data/prep/annotations.bed --cores 2`

`atlatl breakends --alignment tests/data/align/alignment.bam --bedfile tests/data/prep/annotations.bed --out_dir tests/data/breakends --min_alignments 5`

`atlatl group_assemble_and_visualize --breakends tests/data/breakends/construct_a.tsv --annotations tests/data/prep/annotations.bed --reads tests/data/input/reads.fastq --reference tests/data/prep/ref_aug.fasta --outdir tests/data/assembly --vizdir tests/data/visualizations --threads 4`

![Continuous Integration Status](https://github.com/bihealth/clear-CNV/workflows/CI/badge.svg)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/2eaafb57fbb74a46b918e9f58142c880)](https://www.codacy.com?utm_source=github.com&utm_medium=referral&utm_content=bihealth/clear-CNV&utm_campaign=Badge_Coverage)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/2eaafb57fbb74a46b918e9f58142c880)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=bihealth/clear-CNV&amp;utm_campaign=Badge_Grade)
[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
