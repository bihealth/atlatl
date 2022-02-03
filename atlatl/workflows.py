#%%
import subprocess
import pathlib
import yaml
from logzero import logger

#%%

#%%
def workflow_prep(
        workdir:str,
        alleles:list,
        fastas:list,
        reference:str,
        add_annotations:str=None,
        cores:int=4,
        **kwargs):
    """
    Creates a configfile and runs snakemake to prepare read mapping of an augmented reference with additional
    alleles which are created from given subsequences (fasta-files) and bed-files.
    
    workdir: workdimap_ontr of snakemake
    alleles: list of paths to bed-files that define all additional alleles based on concatenations of sequences given in the provided "fastas"
    fastas: list of fasta-files that contain the sequences from which alleles defined in "alleles" are constructed
    reference: path to the reference which will be augmented by the newly constructed alleles and shall be used in downstream read mapping
    add_annotations: path to a bedfile which contains additional annotations on the final alleles and / or the reference
    cores: cores argument passed to snakemake
    """
    # --- CREATE CONFIG --- #
    d = {
        'workdir': workdir,
        'alleles': alleles,
        'fastas':  fastas,
        'reference': reference,
        'add_annotations': add_annotations
    }
    # --- generate path --- #
    configpath = (
        pathlib.Path(__file__).absolute().parent
        / pathlib.Path("workflow")
        / pathlib.Path("config_prep.yml")
    )
    
    # --- save config --- #
    with open(configpath, "wt") as f:
        yaml.safe_dump(d,stream=f,default_flow_style=False)
    
    # --- run snakemake with config --- #
    arguments = [
        "snakemake",
        "--snakefile",
        str(
            pathlib.Path(__file__).absolute().parent
            / pathlib.Path("workflow")
            / pathlib.Path("Snakefile_prep")
        ),
        "--configfile",
        str(configpath),
        "--cores",
        str(cores),
    ]
    logger.info("calling snakemake via \n"+' '.join(arguments))
    subprocess.check_call(arguments)

def workflow_map_ont(
        workdir:str,
        in_reads:str,
        constructs: str,
        reference:str,
        annotations:str,
        depth:int=0,
        cores:int=4,
        minimap_r:int=95,
        minimap_z:int=95,
        **kwargs):
    """
    Creates a configfile and runs snakemake to map long reads to the augmented reference and produces feature counts for the given features.
    Run this after workflow_prep.
    
    workdir: workdir of snakemake
    in_reads: fastq-file of reads
    reference: path to augmented reference
    annotations: path to generated bed-file with annotations for the augmented reference
    cores: cores argument passed to snakemake
    """
    # --- CREATE CONFIG --- #
    d = {
        'workdir': workdir,
        'constructs': constructs,
        'in_reads': in_reads,
        'reference': reference,
        'annotations': annotations,
        'depth': depth,
        'minimap_r': minimap_r,
        'minimap_z': minimap_z
    }
    # --- generate path --- #
    configpath = (
        pathlib.Path(__file__).absolute().parent
        / pathlib.Path("workflow")
        / pathlib.Path("config_map_ont.yml")
    )
    
    # --- save config --- #
    with open(configpath, "wt") as f:
        yaml.safe_dump(d,stream=f,default_flow_style=False)
    
    # --- run snakemake with config --- #
    arguments = [
        "snakemake",
        "--snakefile",
        str(
            pathlib.Path(__file__).absolute().parent
            / pathlib.Path("workflow")
            / pathlib.Path("Snakefile_map_ont")
        ),
        "--configfile",
        str(configpath),
        "--cores",
        str(cores),
    ]
    logger.info("calling snakemake via \n"+' '.join(arguments))
    subprocess.check_call(arguments)




def workflow_map_is(
        workdir:str,
        reads_table:str,
        reference:str,
        annotations:str,
        cores:int=4,
        **kwargs):
    """
    Creates a configfile and runs snakemake to map short reads to the augmented reference and produces feature counts for the given features.
    Run this after workflow_prep.
    
    workdir: workdir of snakemake
    reads_table: path to a file containing rows: name\\tpath_to_R1.fastq\\tpath_to_R2.fastq
    reference: path to augmented reference
    annotations: path to generated bed-file with annotations for the augmented reference
    cores: cores argument passed to snakemake
    """
    # --- CREATE CONFIG --- #
    d = {
        'workdir': workdir,
        'reads_table': reads_table,
        'reference': reference,
        'annotations': annotations
    }
    # --- generate path --- #
    configpath = (
        pathlib.Path(__file__).absolute().parent
        / pathlib.Path("workflow")
        / pathlib.Path("config_map_is.yml")
    )
    
    # --- save config --- #
    with open(configpath, "wt") as f:
        yaml.safe_dump(d,stream=f,default_flow_style=False)
    
    # --- run snakemake with config --- #
    arguments = [
        "snakemake",
        "--snakefile",
        str(
            pathlib.Path(__file__).absolute().parent
            / pathlib.Path("workflow")
            / pathlib.Path("Snakefile_map_is")
        ),
        "--configfile",
        str(configpath),
        "--cores",
        str(cores),
    ]
    logger.info("calling snakemake via \n"+' '.join(arguments))
    subprocess.check_call(arguments)
