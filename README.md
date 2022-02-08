`atlatl workflow_prep --workdir . --alleles input/construct.bed --fastas input/sequences.fasta --reference input/ref.fasta --add_annotations input/annotations.bed --cores 2`

`atlatl workflow_map_ont --workdir . --constructs prep/constructs.fasta --in_reads input/reads.fastq --reference prep/ref_aug.fasta --annotations prep/annotations.bed --cores 2`



![Continuous Integration Status](https://github.com/bihealth/clear-CNV/workflows/CI/badge.svg)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/2eaafb57fbb74a46b918e9f58142c880)](https://www.codacy.com?utm_source=github.com&utm_medium=referral&utm_content=bihealth/clear-CNV&utm_campaign=Badge_Coverage)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/2eaafb57fbb74a46b918e9f58142c880)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=bihealth/clear-CNV&amp;utm_campaign=Badge_Grade)
[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
