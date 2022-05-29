#!/usr/bin/env python

import argparse
import os
from logzero import logger

from . import __version__
from . import tools
from . import workflows

#: The executables required for running atlatl.

def get_parser():
    """Return argparse command line parser."""
    parser = argparse.ArgumentParser(
        description="atlatl is a toolkit for targeted long read based genotyping for structural variations at candidate loci.",
    )
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers()

    # =========================================================================
    #  fasta_to_bed
    # =========================================================================
    parser_fasta_to_bed = subparsers.add_parser("fasta_to_bed",description="Saves a BED-file with a row for each sequences like this: sequence_name 0   len(sequence).")
    parser_fasta_to_bed.add_argument("-i","--input", type=os.path.abspath,required=True , help="input: path to input .fasta file.")
    parser_fasta_to_bed.add_argument("-o","--output", type=os.path.abspath,required=True , help="output: path to output .bed file.")
    parser_fasta_to_bed.set_defaults(func=tools.fasta_to_bed)
    
    # =========================================================================
    #  construct_and_save_sequences
    # =========================================================================
    parser_construct_and_save_sequences = subparsers.add_parser("construct_and_save_sequences",description="This script assembles allels according to given .bed-formatted files from a given set of sequences (fasta files). Creates indexes if necessary.")
    parser_construct_and_save_sequences.add_argument("-s","--sequences", type=os.path.abspath,required=True , nargs='+', help="paths to sequence files in fasta format.")
    parser_construct_and_save_sequences.add_argument("-b","--bedfiles", type=os.path.abspath,required=True , nargs='+', help="paths to .bed-files, each defining a final allele.")
    parser_construct_and_save_sequences.add_argument("-o","--outfasta", type=os.path.abspath,required=True , help="path to output file in fasta format.")
    parser_construct_and_save_sequences.add_argument("-e","--outbed", type=os.path.abspath, default=None, help="Optional: path to output file in bed format.")
    parser_construct_and_save_sequences.set_defaults(func=tools.construct_and_save_sequences)

    # =========================================================================
    #  write_fastq_of_transitive_reads
    # =========================================================================
    parser_write_fastq_of_transitive_reads = subparsers.add_parser("write_fastq_of_transitive_reads",description="Generates a fastq file of reads that align transitively (to reads that align) to a target region.")
    parser_write_fastq_of_transitive_reads.add_argument("-a","--alleles", type=os.path.abspath,required=True , help="input: path to the reference. fasta file with all regions of interest.")
    parser_write_fastq_of_transitive_reads.add_argument("-r","--reads", type=os.path.abspath,required=True , help="input: path to fastq-file with all ONT reads.")
    parser_write_fastq_of_transitive_reads.add_argument("-o","--output", type=os.path.abspath,required=True , help="output: path to output fastq file.")
    parser_write_fastq_of_transitive_reads.add_argument("-c","--cores", type=int, default=8, help="param: number of cores used by minimap2. Default is 8")
    parser_write_fastq_of_transitive_reads.add_argument("-d","--depth", type=int, default=0, help="param: maximum recursive depth in all-vs-all mapping and read selection. 0 is default.")
    parser_write_fastq_of_transitive_reads.set_defaults(func=tools.write_fastq_of_transitive_reads)

    # =========================================================================
    #  group_assemble_and_visualize
    # =========================================================================
    #parser_group_assemble_and_visualize = subparsers.add_parser("group_assemble_and_visualize",description="Generates a fastq file of reads that align transitively (to reads that align) to a target region.")
    #parser_group_assemble_and_visualize.add_argument("-a","--alleles", type=os.path.abspath,required=True , help="input: path to the reference. fasta file with all regions of interest.")
    #parser_group_assemble_and_visualize.add_argument("-r","--reads", type=os.path.abspath,required=True , help="input: path to fastq-file with all ONT reads.")
    #parser_group_assemble_and_visualize.add_argument("-o","--output", type=os.path.abspath,required=True , help="output: path to output fastq file.")
    #parser_group_assemble_and_visualize.add_argument("-c","--cores", type=int, default=8, help="param: number of cores used by minimap2. Default is 8")
    #parser_group_assemble_and_visualize.add_argument("-d","--depth", type=int, default=0, help="param: maximum recursive depth in all-vs-all mapping and read selection. 0 is default.")
    #parser_group_assemble_and_visualize.set_defaults(func=tools.write_fastq_of_transitive_reads)

    # =========================================================================
    #  ontarget_QC
    # =========================================================================
    parser_ontarget_QC = subparsers.add_parser("ontarget_QC",description="This functon calculates some descriptive statistics about the aligned reads. A table is printed to file and a scatter plot is saved.")
    parser_ontarget_QC.add_argument("-a","--alignments_path", type=os.path.abspath,required=True , help="path to e.g. alns.bam")
    parser_ontarget_QC.add_argument("-n","--annotations_path", type=os.path.abspath,required=True , help="path to bed-file that stores all annotations on all relevant targets")
    parser_ontarget_QC.add_argument("-o","--out_table", type=os.path.abspath,required=True , help="path to output table file (tsv)")
    parser_ontarget_QC.add_argument("-f","--out_fig", type=os.path.abspath,required=True , help="path to output figure file (html)")
    parser_ontarget_QC.set_defaults(func=tools.ontarget_QC)

    # =========================================================================
    #  add_constructs_ro_ref
    # =========================================================================
    parser_add_constructs_ro_ref = subparsers.add_parser("add_constructs_ro_ref",description="Equivalent to cat [input ...] | seqkit seq -w [linewidth] -o [output]")
    parser_add_constructs_ro_ref.add_argument("-i","--input", type=os.path.abspath,required=True , nargs='+', help="paths to files to cat.")
    parser_add_constructs_ro_ref.add_argument("-o","--output", type=os.path.abspath,required=True , help="path to output file.")
    parser_add_constructs_ro_ref.add_argument("-w","--linewidth", type=int, default=60, help="Linewidth of output .fasta file. Default is 60.")
    parser_add_constructs_ro_ref.set_defaults(func=tools.add_constructs_ro_ref)


    # =========================================================================
    #  breakends
    # =========================================================================
    parser_breakends = subparsers.add_parser("breakends",description="Significant breakends are calculated from given alignments per chromosome. A dataframe is printed to file which shows read names vs (breakends, segments). Segments are taken from a provided bed-file which describes the region of interest per chromosome.")
    parser_breakends.add_argument("-a","--alignment", type=os.path.abspath,required=True , help="input: path to indexed alignment file in .bam format.")
    parser_breakends.add_argument("-b","--bedfile", type=os.path.abspath,required=True , help="input: path to .bed file which specifies all regions of interest. IMPORTANT: only alignments covering the regions defined in this file are considered!")
    parser_breakends.add_argument("-o","--out_dir", type=os.path.abspath,required=True , help="output: path to the directory where the final dataframes are written to files.")
    parser_breakends.add_argument("--prefix", type=str, default="", help="prefix of output tables' names. The naming scheme looks like: [out_dir]/[prefix.][chr].tsv")
    parser_breakends.add_argument("--alpha", type=float, default=0.001, help="threshold for approximated significance of breakends. Default is 0.001. Has no effect if min_breakends > 0.")
    parser_breakends.add_argument("--radius", type=int, default=30, help="radius for breakpoints to account for significance. Default is 30.")
    parser_breakends.add_argument("--thresh_covered", type=float, default=0.5, help="threshold for coverage of an alignment per region to be considered covered. Should not be smaller than 0.8. Default is 0.5.")
    parser_breakends.add_argument("--cutoff_clipped", type=int, default=30, help="number of bases that are clipped but are not considered to indicate a break end. Default is 30.")
    parser_breakends.add_argument("--min_breakends", type=int, default=0, help="minimum number of break ends for a single break point to be considered and reported. min_breakends <= 0 means only significant breakends are considered and reported.")
    parser_breakends.add_argument("--all_ends", action='store_true', default=False, help="If given, all alignemnt ends are considered and reported, not only if they are breakends.")
    parser_breakends.add_argument("--no_breakends", action='store_true', default=False, help="If given, report only alignment ends if they are not breakends. Has no effect if all_ends is True.")
    parser_breakends.add_argument("--print_overlaps", action='store_true', default=False, help="If given, write overlapped regions to the output table showing the number and direction of overlaps.")
    parser_breakends.set_defaults(func=tools.generate_breakend_matrices)

    # =========================================================================
    #  visualize_assembly
    # =========================================================================

    parser_visualize_assembly = subparsers.add_parser("visualize_assembly",description="Visualizes alignments of a consensus sequence based on the read as the reference.")
    parser_visualize_assembly.add_argument("--alignments", type=os.path.abspath,required=True , help="input: path to a .bam file containing the mapped assembly")
    parser_visualize_assembly.add_argument("--annotations", type=os.path.abspath,required=True , help="input: path to bed file which annotates all relevant regions.")
    parser_visualize_assembly.add_argument("--outhtml", type=os.path.abspath,required=True , help="output: path to .html file with visualizations")
    parser_visualize_assembly.add_argument("--chrs", type=str,required=False , nargs='*',default=[], help="param: list of chromosome names to be visualized. If left blank, all chromosomes to which an assembly was aligned to are printed.")
    parser_visualize_assembly.add_argument("--thickness", type=int, default=5, help="param: in visualizations - thickness of objects.")
    parser_visualize_assembly.set_defaults(func=tools.visualize_assembly)

    # =========================================================================
    #  assemble_and_visualize
    # =========================================================================

    parser_assemble_and_visualize = subparsers.add_parser("assemble_and_visualize",description="Given a file with one read name per line, a file with reads and a reference, this function computes a consensus of the given reads and alignes it to the reference. Helpful for genotyping of a mixed read origins.")
    parser_assemble_and_visualize.add_argument("--readgroup", type=os.path.abspath,required=True , help="input: path to a text-file containing one read name per line.")
    parser_assemble_and_visualize.add_argument("--annotations", type=os.path.abspath,required=True , help="input:  path to bed file which annotates all relevant regions. IMPORTANT: only reads that align to the reference in the annotated regions are considered.")
    parser_assemble_and_visualize.add_argument("--reads", type=os.path.abspath,required=True , help="input: path to fastq(.gz) file containing all reads from which the named reads are selected.")
    parser_assemble_and_visualize.add_argument("--reference", type=os.path.abspath,required=True , help="input: path to the indexed reference file in fasta(.gz) format.")
    parser_assemble_and_visualize.add_argument("--outfiles_prefix", type=os.path.abspath,required=True , help="output: prefix of all output files.")
    parser_assemble_and_visualize.add_argument("--chrs", type=str,required=False , nargs='*',default=[], help="param: list of chromosome names to be visualized. If left blank, all chromosomes to which an assembly was aligned to are printed.")
    parser_assemble_and_visualize.add_argument("--threads", type=int, default=3, help="param: number of threads used by lamassemble. Default is 3.")
    parser_assemble_and_visualize.add_argument("--fasta_name", type=str,required=False, default='consensus' , help="param: Name of the sequence in the fasta file.")
    parser_assemble_and_visualize.add_argument("--technology", type=str,required=False, default='ont' , help="param: Name of the used technology. Choose from ONT: 'ont', PacBio: 'pb'. Default is 'ont'.")
    parser_assemble_and_visualize.add_argument("--image_format", type=str,required=False, default='html' , help="param: Image format name. Choose from [html,pdf,png]. Default is 'html'.")
    parser_assemble_and_visualize.add_argument("--thickness", type=int, default=5, help="param: in visualizations - thickness of objects.")
    parser_assemble_and_visualize.set_defaults(func=tools.assemble_and_visualize)

    # =========================================================================
    #  group_assemble_and_visualize
    # =========================================================================
    parser_group_assemble_and_visualize = subparsers.add_parser("group_assemble_and_visualize",description="Given the breakends table, this function extracts all informative read groups, assembles, and visualizes each of them.")
    parser_group_assemble_and_visualize.add_argument("--breakends", type=os.path.abspath,required=True , help="path to the .tsv file containing all breakends. This can be computed using atlatl generate_breakend_matrices.")
    parser_group_assemble_and_visualize.add_argument("--annotations", type=os.path.abspath,required=True , help="input:  path to bed file which annotates all relevant regions. IMPORTANT: only reads that align to the reference in the annotated regions are considered.")
    parser_group_assemble_and_visualize.add_argument("--reads", type=os.path.abspath,required=True , help="input: path to fastq(.gz) file containing all reads from which the named reads are selected.")
    parser_group_assemble_and_visualize.add_argument("--reference", type=os.path.abspath,required=True , help="input: path to the indexed reference file in fasta(.gz) format.")
    parser_group_assemble_and_visualize.add_argument("--outdir", type=os.path.abspath,required=True , help="output: directory to where all assembly and visualizations related data are saved.")
    parser_group_assemble_and_visualize.add_argument("--chrs", type=str,required=False , nargs='*',default=[], help="param: list of chromosome names to be visualized. If left blank, all chromosomes to which an assembly was aligned to are printed.")
    parser_group_assemble_and_visualize.add_argument("--threads", type=int, default=3, help="param: number of threads used by lamassemble. Default is 3.")
    parser_group_assemble_and_visualize.add_argument("--prefix", type=str,required=False, default="", help="param: prefix of all generated files.")
    parser_group_assemble_and_visualize.add_argument("--minimum_reads", type=int, default=5, help="param: minimum number of reads that are required to form an informative group to be assembled and visualized. Default is 5.")
    parser_group_assemble_and_visualize.add_argument("--fasta_name", type=str,required=False, default='consensus' , help="param: Name of the sequence in the fasta file.")
    parser_group_assemble_and_visualize.add_argument("--technology", type=str,required=False, default='ont' , help="param: Name of the used technology. Choose from ONT: 'ont', PacBio: 'pb'. Default is 'ont'.")
    parser_group_assemble_and_visualize.add_argument("--thickness", type=int, default=5, help="param: in visualizations - thickness of objects. Default is 5.")
    parser_group_assemble_and_visualize.add_argument("--img_format", type=str, default='html', help="param: output image in one of the formats: html, pdf, png")
    parser_group_assemble_and_visualize.set_defaults(func=tools.group_assemble_and_visualize)

    # =========================================================================
    #  WORKFLOW workflow_prep
    # =========================================================================
    workflow_prep_ont = subparsers.add_parser("workflow_prep",description="Prepares successive read mapping given artificial alleles and the reference.")
    workflow_prep_ont.add_argument("-w","--workdir", type=os.path.abspath,required=True , help="workdir of snakemake.")
    workflow_prep_ont.add_argument("-a","--alleles", type=os.path.abspath,required=True , nargs='+', help="list of paths to bed-files that define all additional alleles based on concatenations of sequences given in the provided 'fastas'.")
    workflow_prep_ont.add_argument("-s","--sequences", type=os.path.abspath,required=True , nargs='+', help="list of fasta-files that contain the sequences from which alleles defined in 'alleles' are constructed.")
    workflow_prep_ont.add_argument("-r","--reference", type=os.path.abspath,required=True , help="path to the reference which will be augmented by the newly constructed alleles and shall be used in downstream read mapping.")
    workflow_prep_ont.add_argument("-n","--add_annotations", type=os.path.abspath,required=True , help="path to a bedfile which contains additional annotations on the final alleles and / or the reference.")
    workflow_prep_ont.add_argument("-c","--cores", type=int, default=30, help="cores argument passed to snakemake.")
    workflow_prep_ont.set_defaults(func=workflows.workflow_prep)

    # =========================================================================
    #  WORKFLOW workflow_map
    # =========================================================================
    workflow_map = subparsers.add_parser("workflow_map",description="Creates a configfile and runs snakemake map reads to the augmented reference and produces feature counts for the given features and repeats. Run this after workflow_prep_ont.")
    workflow_map.add_argument("-w","--workdir",     type=os.path.abspath,   required=True, help="workdir of snakemake.")
    workflow_map.add_argument("-a","--constructs",  type=os.path.abspath,   required=True, help="fasta-file of constructed sequences of the final alleles.")
    workflow_map.add_argument("-i","--in_reads",    type=os.path.abspath,   required=True, help="fastq-file of reads.")
    workflow_map.add_argument("-r","--reference",   type=os.path.abspath,   required=True, help="path to augmented reference.")
    workflow_map.add_argument("-n","--annotations", type=os.path.abspath,   required=True, help="path to generated bed-file with annotations for the augmented reference.")
    workflow_map.add_argument("-c","--cores",       type=int,   default=30,   help="cores argument passed to snakemake.")
    workflow_map.add_argument("-d","--depth",       type=int,   default=0,    help="depth of recursion in read selection. Only useful if great re-arrangements of the knock-in are expected. Default is 0.")
    workflow_map.add_argument("--minimap_r",        type=int,   default=95,   help="-r parameter of minimap2. Deletions greater than this value are detected in breakend-analysis. Smaller deletions remain within single alignments. Default is 95. This is feasible for read error-rates of less than 5 percent. Use greater values for higher error-rates.")
    workflow_map.add_argument("--minimap_z",        type=int,   default=95,   help="-z parameter of minimap2. Inversions greater than this values can be detected. Default is 95. This is feasible for read error-rates of less than 5 percent. Use greater values for higher error-rates.")
    workflow_map.set_defaults(func=workflows.workflow_map)

#    workflow_map.add_argument("-p","--in_repeats", type=os.path.abspath,required=True , help="bed-file of reads for given reference.")

    # =========================================================================
    #  WORKFLOW workflow_map_is
    # =========================================================================
    #workflow_map_is = subparsers.add_parser("workflow_map_is",description="Creates a configfile and runs snakemake to map short reads to the augmented reference and produces feature counts for the given features and repeats. Run this after workflow_prep.")
    #workflow_map_is.add_argument("-w","--workdir", type=os.path.abspath,required=True , help="workdir of snakemake.")
    #workflow_map_is.add_argument("-i","--reads_table", type=os.path.abspath,required=True,  help="fastq-file of reads.")
    #workflow_map_is.add_argument("-r","--reference", type=os.path.abspath,required=True,  help="path to augmented reference.")
    #workflow_map_is.add_argument("-n","--annotations", type=os.path.abspath,required=True,  help="path to generated bed-file with annotations for the augmented reference.")
    #workflow_map_is.add_argument("-c","--cores", type=int, default=30, help="cores argument passed to snakemake.")
    #workflow_map_is.set_defaults(func=workflows.workflow_map_is)

    # =========================================================================
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
    else:
        logger.info("Starting analysis...")
        args.func(**vars(args))
        logger.info("All done. Have a nice day!")


if __name__ == "__main__":
    main()
