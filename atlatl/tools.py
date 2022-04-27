import sys
import argparse
import pathlib
import pandas as pd
import numpy as np
import shlex
import tempfile
import subprocess
from logzero import logger
from Bio import SeqIO

from . import helpers

def fasta_to_bed(
        input:str,
        output:str,
        **kwargs):
    """
    Saves a BED-file with a row for each sequence like this: seq_name 0   len(sequence)   seq_name.

    input: path to input .fasta file
    output: path to output .bed file
    """
    with open(output,'wt') as outf:
        with open(input,'r') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                print('\t'.join(map(str,[record.id,0,len(record.seq),record.id])),file=outf)

def construct_and_save_sequences(
        bedfiles:list,
        sequences:list,
        outfasta:str,
        outbed:str=None,
        **kwargs):
    """
    Receives several (a list of) .fasta files and .bed (tab separated) files, where the .bed files
    define the concatenation of substrings of the given sequences in the .fasta files.
    Genomic coordinates where start > end are reversed. Coordinates are 1-based.
    The output file is a .fasta file with all assembled alleles.
    If defined, a bed-file will be written that annotates all underlying sub-sequences.
    """
    sequences,annotations = helpers.assemble_alleles_from_bed_and_fasta(bedfiles, sequences)
    SeqIO.write(sequences, pathlib.Path(outfasta), "fasta")
    if outbed:
        annotations.to_csv(outbed,sep='\t',header=False,index=False)

def write_fastq_of_transitive_reads(
        reads:str,
        alleles:str,
        output:str,
        cores:int=4,
        depth:int=0,
        **kwargs):
    """
    Writes a fastq-file with all reads that are in a transitive relation of alignment with a set
    of reads that can be aligned to a given set of alleles. Returns a tuple:
    1. file handle of a temporary sam file.
    2. list of read names that are in transitive relation to the original target region of the reference. 
    
    reads: path to fastq-formatted reads 
    alleles: path to fasta-formatted file with expected subsequences (alleles) 
    output: path to putput fastq file.
    cores: number of threads used by minimap2 
    depth: max number of recursions in the ava-read mapping expansion 
    """
    with open(output,'w+') as outf:
        tmp_selected_reads_fastq = helpers.list_to_fastq(reads,helpers.select_and_expand_reads(reads,alleles,cores,depth))
        for line in tmp_selected_reads_fastq:
            outf.write(line)

def add_constructs_ro_ref(
        input:list,
        output:str,
        linewidth:int,
        **kwargs):
    """
    Equivalent to cat [input ...] | seqkit seq -w [linewidth] -o [output]

    input: list of paths to fasta-files
    ouput: path to output fasta-file
    linewidth: line width of output fasta file
    """
    with tempfile.NamedTemporaryFile(mode='w+') as fp:
        helpers.cat_files(input,fp.name)
        helpers.format_fasta(fp.name,output,linewidth)

def generate_breakend_matrices(
        alignment:str,
        bedfile:str,
        out_dir:str,
        prefix:str="",
        alpha:float=0.01,
        min_breakends:int=0,
        radius:int=20,
        thresh_covered:float=0.5,
        cutoff_clipped:int=30,
        all_ends:bool=False,
        no_breakends=False,
        print_overlaps=False,
        **kwargs):
    """
    Significant breakends are calculated from given alignments per chromosome.
    A dataframe is printed to file which shows read names vs (breakends, segments).
    Segments are taken from a provided bed-file which describes the region of interest per chromosome.

    alignment: path to indexed alignment file in .bam format
    bedfile: path to .bed file which specifies all regions of interest. IMPORTANT: only alignments covering the regions defined in this file are considered!
    out_dir: path to the directory where the final dataframes are written to files
    alpha: threshold for significance of breakends
    
    radius: radius for breakpoints to account for significance
    thresh_covered: threshold for coverage of an alignment per region to be considered covered. Should not be smaller than 0.8
    cutoff_clipped: number of bases that are clipped but are not considered to indicate a break end
    all_ends: if True, get all ends of alignments, not only breakends.
    no_breakends: report only alignment ends if they are not breakends. Has no effect if all_ends is True."""
    helpers.print_breakends_and_overlaps(alignment_path=alignment,
                                        bed_path=bedfile,
                                        save_dir=out_dir,
                                        prefix=prefix,
                                        alpha=alpha,
                                        min_breakends=min_breakends,
                                        r=radius,thresh=thresh_covered,
                                        cutoff=cutoff_clipped,
                                        all_ends=all_ends,
                                        no_be=no_breakends,
                                        print_overlaps=print_overlaps)

def visualize_assembly(alignments:str,annotations:str,outhtml:str,chrs:list,thickness:int=5,**kwargs):
    helpers.visualize_assembly(alignments_path=alignments,
                            annotations_path=annotations,
                            outpath=outhtml,
                            chrs=chrs,
                            thickness=thickness)

def assemble_and_visualize(readgroup:str,
                            annotations:str,
                            reads:str,
                            reference:str,
                            outfiles_prefix:str,
                            chrs:str,
                            prefix:str="",
                            image_format="html",
                            threads:int=3,
                            fasta_name:str="",
                            technology:str='ont',
                            thickness:int=5,
                            **kwargs):
    """
    Given a file with one read name per line, a file with reads and a reference, this function computes a consensus of the given reads
    and alignes it to the reference. Helpful for genotyping of a mixed read origins.

    readgroup: path to a text-file containing one read name per line
    annotations: path to bed file which annotates all relevant regions. IMPORTANT: only reads that align to the reference in the annotated regions are considered.
    reads: path to fastq(.gz) file containing all reads from which the named reads are selected
    reference: path to the indexed reference file in fasta(.gz) format.
    name_forward: name of final html output of forward assembly.
    name_reverse: name of final html output of reverse complement assembly.
    fastaout: path to the file to which the final consensus will be written.
    fastaout_reversed: same as [fastaout] but for the reverse complement sequence.
    bamout: path to the final compressed and indexed alignment file (/path/to/file.bam).
    bamout_reversed: same as [bamout] but for the reverse complement sequence.
    threads: number of threads used by lamasemble.
    fasta_name: Name of the sequence in the fasta file.
    technology: name of the used technology. Choose from 'ont' or 'pb'.
    space_read_to_ref: in visualizations: space between read and reference. Please increase if overlapping visuals occur.
    space_per_chr: in visualizations: space between each chromosme. Sould be greater than [space_read_to_ref]. Please increase if overlapping visuals occur.
    thickness: in visualizations: thickness of objects.
    """
    helpers.assemble_and_visualize(
                        readgroup=readgroup,
                        annotations=annotations,
                        reads=reads,
                        reference=reference,
                        outfiles_prefix=outfiles_prefix,
                        chrs=chrs,
                        threads=threads,
                        fasta_name=fasta_name,
                        technology=technology,
                        thickness=thickness,
                        image_format=image_format)

def group_assemble_and_visualize(
                            breakends:str,
                            annotations:str,
                            reads:str,
                            reference:str,
                            outdir:str,
                            chrs:list,
                            threads:int=3,
                            prefix:str="",
                            minimum_reads:int=5,
                            fasta_name:str="",
                            technology:str='ont',
                            thickness:int=5,
                            img_format='html',
                            **kwargs):
    """
    Given the breakends table, this function extracts all informative read groups, assembles, and visualizes each of them.

    breakends: path to the .tsv file containing all breakends. This can be computed using atlatl generate_breakend_matrices.
    annotations: path to bed file which annotates all relevant regions. IMPORTANT: only reads that align to the reference in the annotated regions are considered.
    reads: path to fastq(.gz) file containing all reads from which the named reads are selected
    reference: path to the indexed reference file in fasta(.gz) format.
    outdir: directory to where all assembly related data are saved.
    vizdir: directory to where all visualizations are saved.
    threads: number of threads used by lamasemble.
    prefix: prefix of all generated files.
    minimum_reads: minimum number of reads that are required to form an informative group to be assembled and visualized.
    fasta_name: Name of the sequence in the fasta file.
    technology: name of the used technology. Choose from 'ont' or 'pb'.
    space_read_to_ref: in visualizations: space between read and reference. Please increase if overlapping visuals occur.
    space_per_chr: in visualizations: space between each chromosme. Sould be greater than [space_read_to_ref]. Please increase if overlapping visuals occur.
    thickness: in visualizations: thickness of objects.
    """
    helpers.group_assemble_and_visualize(                            
                            breakends=breakends,
                            annotations=annotations,
                            reads=reads,
                            reference=reference,
                            outdir=outdir,
                            chrs=chrs,
                            threads=threads,
                            prefix=prefix,
                            minimum_reads=minimum_reads,
                            fasta_name=fasta_name,
                            technology=technology,
                            thickness=thickness,
                            img_format=img_format)
