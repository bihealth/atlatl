
import pathlib
import subprocess
import shlex
from collections import defaultdict
from typing import Tuple
import pandas as pd
import numpy as np
import tempfile
import io
import os
from scipy.stats import binom
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from logzero import logger
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam
from tqdm import tqdm

# =============================================================================
#  construct_and_save_sequences

def create_index(
        path:str):
    """Checks if a .fai or .bai file exists else creates one with samtools.
    path
    """
    if pathlib.Path(path).suffix == '.fa' or pathlib.Path(path).suffix == '.fasta':
        if not pathlib.Path(path+'.fai').is_file():
            subprocess.check_call(['samtools','faidx',path])
        return path+'.fai'
    if pathlib.Path(path).suffix == '.bam':
        if not pathlib.Path(path+'.bai').is_file():
            subprocess.check_call(['samtools','index',path])
        return path+'.bai'
    raise IOError(f"File {path} is not .fa,.fasta or .bam.")

def create_chr_dict(
        fasta_paths:list):
    """Creates a dictionary \{chr:path_to_fasta\} for each given fasta file and indexed chromosome in it.
    Indexes fasta files automatically if necessary. Calls create_index internally.
    fasta_paths: list of paths to fasta-files
    """
    d = defaultdict(str)
    for file in fasta_paths:
        for line in open(create_index(file)):
            d[line.split('\t')[0]] = file
    return d

def bed_to_features_bed(
        input:str,
        output:str=None,
        **kwargs) -> pd.DataFrame:
    """
    Creates a bed file describing the subsequences of a given bedfile after the corresponding fasta-file was assembled (concatenated)
    Returns a pandas DataFrame of the new coordinates.
    
    input: path to bed-file of subsequences that form the allele if concatenated
    ouput: path to the bedfile annotating the assembled fasta-file
    """
    bed = pd.read_csv(input,sep='\t',header=None)
    allele_name = pathlib.Path(input).stem

    end = (bed.iloc[:,2]-bed.iloc[:,1]).cumsum().to_numpy()
    start = np.array([0,*end[:-1]])

    allel = pd.DataFrame([allele_name]*bed.shape[0],columns=["chr"])
    allel.loc[:,"start"] = start
    allel.loc[:,"end"] = end
    allel.loc[:,'id'] = bed.iloc[:,3]
    if output:
        allel.loc[:,['chr','start','end','id']].to_csv(pathlib.Path(output),sep='\t',header=False,index=False)
    return allel.loc[:,['chr','start','end','id']]

def assemble_alleles_from_bed_and_fasta(
        bedpaths:list,
        fasta_paths:list) -> tuple:
    """
    Concatenates each target defined in one bed-file to a final allele.
    Each region must exist in the given fasta-files. The final output is a list of Bio.SeqRecord objects.
    Creates fasta indexes via samtools if necessary.
    """
    C = create_chr_dict(fasta_paths)
    R = []
    D = pd.DataFrame()
    for bed in bedpaths:
        allele_name = pathlib.Path(bed).stem
        S = Seq('')
        D = pd.concat([D,bed_to_features_bed(bed)])
        bdf = pd.read_csv(pathlib.Path(bed),sep='\t',header=None)
        for i in range(bdf.shape[0]):
            s = Seq('')
            if int(bdf.iloc[i,1]) <= int(bdf.iloc[i,2]):
                locus = str(bdf.iloc[i,0])+':'+str(bdf.iloc[i,1])+'-'+str(bdf.iloc[i,2])
                proc = subprocess.Popen(['samtools','faidx',C[bdf.iloc[i,0]],locus], stdout=subprocess.PIPE)
                read = False
                for line in io.TextIOWrapper(proc.stdout, encoding="utf-8"):
                    if read:
                        s += Seq(line.rstrip())
                    read = True
                S += s
            else:
                locus = str(bdf.iloc[i,0])+':'+str(bdf.iloc[i,2])+'-'+str(bdf.iloc[i,1])
                proc = subprocess.Popen(['samtools','faidx',C[bdf.iloc[i,0]],locus], stdout=subprocess.PIPE)
                read = False
                for line in io.TextIOWrapper(proc.stdout, encoding="utf-8"):
                    if read:
                        s += Seq(line.rstrip())
                    read = True
                S += s[::-1]
        R.append(
            SeqRecord(
                S,
                id=allele_name,
                name=allele_name,
                description='Artificial chromosome composed of '+','.join(bdf.iloc[:,0])
                )
            )
    return R,D

# =============================================================================
#  write_fastq_of_transitive_reads(

def list_to_fastq(
        reads:str,
        slist:set) -> tempfile.NamedTemporaryFile:
    """
    Writes records from a reads containing file (fastq) specified in slist (python listof seq names) to a temporary file (output: fastq).
    
    reads: Path to fastq-file with all reads.
    slist: set of read names
    """
    tmp_selected_reads_list  = tempfile.NamedTemporaryFile(mode='w+')
    tmp_selected_reads_fastq = tempfile.NamedTemporaryFile(mode='w+')
    # create list of read names
    with open(tmp_selected_reads_list.name,'wt') as f:
        for line in sorted(slist):
            print(line,file=f)
    # call seqtk to filter fastq for slist
    c1 = shlex.split(f"seqtk subseq {reads} {tmp_selected_reads_list.name}")
    p1 = subprocess.Popen(
        c1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8"
    )
    output, error = p1.communicate()
    with open(tmp_selected_reads_fastq.name,'wt') as f:
        print(output,file=f)
    tmp_selected_reads_list.close()
    return tmp_selected_reads_fastq

def map_ava(
    reads_query:str,
    reads_ref:str,
    cores:int=4):
    """
    Maps reads to sequences in .fastq format using minimap2 all-vs-all alignment. The output is written to a tmp file hand the file handle returned.
    
    reads_query: path to fastq-formatted reads. 
    reads_ref: path to fastq-formatted reads. 
    cores: number of threads used by minimap2. 
    """
    tmp_alignment = tempfile.NamedTemporaryFile(mode='w+')

    # map selected reads on reads loci (alleles.fasta)
    c0 = shlex.split(f"minimap2 -t {cores} -a -x ava-ont {reads_ref} {reads_query} -o {tmp_alignment.name}")
    # extract unique reads
    p0 = subprocess.Popen(
        c0, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8"
    )
    output, error = p0.communicate()
    logger.info(error)
    return tmp_alignment

def map_on_alleles(
    alleles:str,
    reads:str,
    cores:int=4):
    """
    Maps reads to sequences in .fasta format using minimap2. The output is written to a temp file handle and returned.
    
    alleles: path to fasta-formatted file with expected subsequences (alleles) 
    reads: path to fastq-formatted reads 
    cores: number of threads used by minimap2 
    """
    tmp_alignment = tempfile.NamedTemporaryFile(mode='w+')

    # map reads on extracted loci (alleles.fasta)
    c0 = shlex.split(f"minimap2 -t {cores} -a -x map-ont {alleles} {reads} -o {tmp_alignment.name}")
    # extract unique reads
    p0 = subprocess.Popen(
        c0, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8"
    )
    output, error = p0.communicate()
    return tmp_alignment

def select_reads_on_alignment(
        align:str) -> set:
    """
    Creates a set of read names from a given alignment.
    
    align: path to .sam-formatted file containing the final alignment.

    returns a set of read names.
    """
    tmp_alignment = tempfile.NamedTemporaryFile(mode='w+')
    c1 = shlex.split(f"samtools view -F 4 {align} -o {tmp_alignment.name}")
    p1 = subprocess.Popen(
        c1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8"
    )
    output, error = p1.communicate()
    return set(pd.read_csv(tmp_alignment.name,usecols=[0],sep='\t',header=None).loc[:,0])

def get_names_of_fastq(reads:str) -> set:
    """
    Returns a list of strings which are the names of the records in a fastq-file.
    
    reads: Path to fastq-file.
    """
    c1 = shlex.split(f"seqkit seq {reads} -i -n")
    p1 = subprocess.Popen(
        c1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8"
    )
    output, error = p1.communicate()
    return set(output.split('\n'))

def select_and_expand_reads(
    reads:str,
    alleles:str,
    cores:int=4,
    depth:int=0) -> set:
    """
    Selects all reads that are in a transitive relation of alignment with a set
    of reads that can be aligned to a given set of alleles. Returns a tuple:
    1. file handle of a temporary sam file.
    2. list of read names that are in transitive relation to the original target region of the reference. 
    
    reads: path to fastq-formatted reads 
    alleles: path to fasta-formatted file with expected subsequences (alleles) 
    cores: number of threads used by minimap2 
    depth: max number of recursions in the ava-read mapping expansion 
    """
    C = get_names_of_fastq(reads)
    alignment = map_on_alleles(alleles,reads,cores)
    B = select_reads_on_alignment(alignment.name)
    alignment.close()
    C = C.difference(B)
    A = B.copy()
    # recursive step
    counter = 0
    if depth > 0:
        while True:
            logger.info("recurse: A=%d, B=%d, C=%d"%(len(A),len(B),len(C)))
            Bfastq, Cfastq = list_to_fastq(reads,B),list_to_fastq(reads,C)
            alignment = map_ava(Cfastq.name,Bfastq.name,cores)
            Bfastq.close()
            Cfastq.close()
            B = select_reads_on_alignment(alignment.name)
            alignment.close()
            C = C.difference(B)
            A = A.union(B)
            counter += 1
            if len(B) == 0 or len(C) == 0 or counter >= depth:
                return A
    else:
        logger.info("no recursion initiated. Using reads aligned directly to locus.")
        return A

# =============================================================================
#  add_alleles_ro_ref

def cat_files(infiles,outfile):
    with open(outfile, 'wt') as outfile:
        for fname in infiles:
            with open(fname) as f:
                for line in f:
                    outfile.write(line)

def format_fasta(SOURCE,FILE,linewidth=60):
    c0 = shlex.split(f"seqkit seq -w {linewidth} {SOURCE} -o {FILE}")
    p0 = subprocess.Popen(
        c0, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8"
    )
    output, error = p0.communicate()

# =============================================================================
#  break end extraction

def breakends_of_aligned_fragment(al_frag,cutoff=30,all_ends=False,no_be=False):
    # cutoff - number of bases that are clipped but not considered as a break end
    break_points=[]
    s,e = al_frag.cigartuples[0], al_frag.cigartuples[-1]
    bs,be = False,False
    if all_ends:
        bs,be = True,True
    else:
        if 4 <= s[0] <= 5 and s[1] > cutoff:
            bs = True
        if 4 <= e[0] <= 5 and e[1] > cutoff:
            be = True
        if no_be:
            bs = not bs
            be = not be
    if bs:
        break_points.append(al_frag.reference_start)
    if be:
        break_points.append(al_frag.reference_end)
    return np.array(break_points,dtype=int)

def bp_is_significant(alpha,break_point_candidates,num_breakpoints,all_BEs):
    return 1-binom.cdf(n=all_BEs,k=num_breakpoints,p=1/break_point_candidates) < alpha

def al_has_be(al,be,r,cutoff,all_ends,no_be):
    for be_al in breakends_of_aligned_fragment(al_frag=al,cutoff=cutoff,all_ends=all_ends,no_be=no_be):
        if be_al <= be + r and be_al >= be - r:
            return True
    return False

def find_significant_BEs(alignments,reflen,reads,r,alpha,min_breakends,cutoff,all_ends,no_be):
    arrs = [breakends_of_aligned_fragment(al_frag=al,cutoff=cutoff,all_ends=all_ends,no_be=no_be) for al in alignments]
    if len(arrs) == 0:
        return []
    x = np.sort(np.concatenate(arrs))
    y = np.zeros(len(x),dtype=int)
    while not (x==y).all():
        y = x
        x = np.apply_along_axis(np.round,0,[np.median(x[(x >= b -r) & (x <= b +r)]) for b in x])
    BEs = set(x.astype(int))
    bp_support = {b:0 for b in BEs}
    for al in alignments:
        for b in bp_support.keys():
            bp_support[b] += al_has_be(al,b,r,cutoff,all_ends=all_ends,no_be=no_be)
    significant_pbs = []
    for b in bp_support.keys():
        if min_breakends > 0:
            if bp_support[b] >= min_breakends:
                significant_pbs.append(b)
        else:
            if bp_is_significant(alpha=alpha,num_breakpoints=bp_support[b],break_point_candidates=len(BEs),all_BEs=len(x)):
                significant_pbs.append(b)
    return significant_pbs

# unused
def get_profiles_of_reads(alignments,significant_pbs,readnames,radius=20,cthresh=0.66,cutoff=30):
    profiles = {n:[] for n in readnames}
    for be in significant_pbs:
        for al in alignments:
            char = 'x' if al.get_overlap(be-radius,be+radius) >= radius * 2 * cthresh else '-'
            if al_has_be(al,be,radius,cutoff,all_ends=False,no_be=False):
                char = '|'
            profiles[al.query_name].append(char)
    return profiles
# unused end

def print_breakends_and_overlaps(alignment_path,bed_path,save_dir,prefix="",chrs=[],alpha=0.01,min_breakends=0,r=30,thresh=0.5,cutoff=30,all_ends=False,no_be=False,print_overlaps=False):
    """
    Significant breakends are calculated from given alignments per chromosome.
    A dataframe is printed to file which shows read names vs (breakends, segments).
    Segments are taken from a provided bed-file which describes the region of interest per chromosome.

    alignment_path: path to indexed alignment file in .bam format
    bed_path: path to .bed file which specifies all regions of interest. IMPORTANT: only alignments covering the regions defined in this file are considered!
    save_dir: path to the directory where the final dataframes are written to files
    prefix: prefix of output file name.
    chrs: list of names of chromosomes if only specific chromosomes are considered.
    alpha: threshold for significance of breakends
    r: radius for breakpoints to account for significance
    thresh: threshold for coverage of an alignment per region to be considered covered. Should not be smaller than 0.8
    cutoff: number of bases that are clipped but are not considered to indicate a break end
    all_ends: print all ends of alignments, not only breakends
    no_be: report only alignment ends if they are not breakends
    print_overlaps: print overlapping regions with direction of alignment
    """
    annotations = pd.read_csv(bed_path, sep='\t',header=None)
    if len(chrs) == 0:
        chrs = sorted(set(annotations.iloc[:,0]))
    else:
        for chr in chrs:
            if not chr in sorted(set(annotations.iloc[:,0])):
                raise "All specified chromosomes must be defined in the given annotations (bed) file"
    for chr in chrs:
        annot = annotations[annotations.iloc[:,0] == chr]
        with tempfile.TemporaryDirectory() as tmpdir:
            aln_file = pathlib.Path(tmpdir) / "aln.bam"
            # CRITICAL: {annot[1].min()}-{annot[2].max()} makes it ncessary to define the whole region in the bed file! Not very good.
            cmd_view = shlex.split(f"samtools view -b {str(alignment_path)} \"{chr}:{annot[1].min()}-{annot[2].max()}\" -o {str(aln_file)}")
            logger.info(' '.join(cmd_view))
            subprocess.call(cmd_view)
            cmd_index = shlex.split(f"samtools index {str(aln_file)}")
            subprocess.call(cmd_index)
            A = np.array([a for a in pysam.AlignmentFile(aln_file, "rb")])
            if len(A) == 0:
                logger.info(f"no reads aligned in region {chr}:{annot[1].min()}-{annot[2].max()}")
                continue
            reflen = max([a.reference_end for a in A]) - min([a.reference_end for a in A])
            reads = np.unique([al.query_name for al in A])
            significant_pbs = find_significant_BEs(alignments=A,reflen=reflen,reads=reads,r=r,alpha=alpha,min_breakends=min_breakends,cutoff=cutoff,all_ends=all_ends,no_be=no_be)
            # create 2 tables:
            named_regions = ['-'.join(map(str,annot.loc[i,[3,1,2]])) for i in annot.index]
            annot.index = named_regions
            df_re = pd.DataFrame([['']*len(named_regions)]*len(reads),index=reads if len(reads) else None,columns=named_regions,dtype=str)
            df_be = pd.DataFrame(np.zeros((len(reads),len(significant_pbs))),index=reads,columns=sorted(significant_pbs),dtype=int)
            for al in tqdm(A):
                for be in significant_pbs:
                    if al_has_be(al=al,be=be,r=r,cutoff=cutoff,all_ends=all_ends,no_be=no_be):
                        df_be.loc[al.query_name,be] += 1
                for region in named_regions:
                    if al.get_overlap(*annot.loc[region,[1,2]]) > thresh*(annot.loc[region,2] - annot.loc[region,1]):
                        df_re.loc[al.query_name,region] += '<' if al.is_reverse else '>'
            df_re.replace('^$','-',regex=True,inplace=True)
            df_be.columns = ['be_'+str(c) for c in df_be.columns]
            logger.info(f"found {len(significant_pbs)} significant breakends.")
            logger.info("printing dataframe to file: "+str(pathlib.Path(save_dir) / ((prefix+'.' if prefix != "" else prefix)+chr+'.tsv')))
            pathlib.Path(save_dir).mkdir(parents=True, exist_ok=True)
            if print_overlaps:
                pd.concat([df_be, df_re], axis=1).to_csv(pathlib.Path(save_dir) / ((prefix+'.' if prefix != "" else prefix)+chr+'.tsv'),sep='\t',header=True,index=True)
            else:
                df_be.to_csv(pathlib.Path(save_dir) / ((prefix+'.' if prefix != "" else prefix)+chr+'.tsv'),sep='\t',header=True,index=True)

# =============================================================================
#  QC after mapping

def ontarget_QC(alignments_path, annotations_path, out_table, out_fig):
    """
    This functon calculates some descriptive statistics about the aligned reads. A table is printed to file and a scatter plot is saved.

    alignments_path: path to e.g. alns.bam
    annotations_path: path to bed-file that stores all annotations on all relevant targets
    out_table: path to output table (tsv)
    out_fig: path to output fig (html)
    """
    out_table = pathlib.Path(out_table)
    out_fig = pathlib.Path(out_fig)
    cmd = f"bedtools merge -d 10 -i {annotations_path}"
    p1 = subprocess.Popen(
            shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8"
        )
    output, _ = p1.communicate()
    targets = list(map(lambda x: x[0]+':'+x[1]+'-'+x[2],[line.split('\t') for line in output.rstrip().split('\n')]))
    D = pd.DataFrame()
    for target in targets:
        selected_alns = tempfile.NamedTemporaryFile(mode='w+',suffix='.bam')
        cmd = f"samtools view -b {alignments_path} {target} -o {selected_alns.name}"
        subprocess.check_call(shlex.split(cmd))
        for a in pysam.AlignmentFile(selected_alns.name, "rb"):
            D = pd.concat([D,pd.DataFrame([a.query_name, target, a.query_length, a.query_alignment_length, a.reference_length]).T])
    D.columns = ["read","target","read length","aligned bp","ref bp"]
    D.index = range(D.shape[0])
    elems = set(D.loc[:,["read","target"]].T.apply(lambda x: '_'.join(x)))
    # --- reduce DF --- #
    D_final = pd.DataFrame(index = elems, columns=["read","target","read length","aligned bp","ref bp"])
    D_final["read length"] = np.zeros(D_final.shape[0])
    D_final["aligned bp"] = np.zeros(D_final.shape[0])
    D_final["ref bp"] = np.zeros(D_final.shape[0])
    D_final = D_final.astype({"read length":"int32","aligned bp":"int32","ref bp":"int32"})
    for i in D.index:
        r,t = D.loc[i,["read","target"]]
        id = '_'.join([r,t])
        D_final.loc[id,"read length"] += int(D.loc[i,"read length"])
        D_final.loc[id,"read"] = r
        D_final.loc[id,"target"] = t
        D_final.loc[id,"aligned bp"] += int(D.loc[i,"aligned bp"])
        D_final.loc[id,"ref bp"] += int(D.loc[i,"ref bp"])
    # --- finalize and save DF --- #
    D_final["target size"] = [int(x.split('-')[1]) - int(x.split('-')[0].split(':')[1]) for x in D_final["target"]]
    D_final["read alignment percentage"] = D_final["aligned bp"] / D_final["read length"]
    D_final["ref percentage"] = D_final["ref bp"] / D_final["target size"]
    out_table.parent.mkdir(parents=True, exist_ok=True)
    if not str(out_table).endswith('.tsv'):
        out_table = str(out_table) + '.tsv'
        logger.warn(f"renamed input file to {out_table}")
    D_final.to_csv(out_table, sep='\t',index=False)
    # --- print figure --- #
    fig_lp = px.scatter(D_final, x="read length", y="read alignment percentage",
        opacity=0.5, color="target", marginal_x="histogram", marginal_y="histogram",
        hover_data=D_final.columns)
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    if not str(out_fig).endswith('.html'):
        out_fig = str(out_fig) + '.html'
        logger.warn(f"renamed input file to {out_fig}")
    fig_lp.write_html(out_fig)

# =============================================================================
#  selective_assembly

def on_target_QC(alignments_path, annotations_path, output):
    # pick each consecutive region R in annotations
    annotations = pd.read_csv(annotations_path,sep='\t')
    annotations.astype({'city':'string', 'num_candidates':'int32'})
    # for r in R:
    #   select all alignments
    #   table.loc[readname,aligned] += aligned bases
    #   table.loc[readname,length]  = read length

# gets a path to breakends.tsv and an out_dir name to isolate all informative readgroups and
# returns the generated filepaths.
def readgroups_to_files(path_bes,out_dir,prefix,minimum_reads=5):
    out_dir = pathlib.Path(out_dir)
    path_bes = pathlib.Path(path_bes)
    # load breakends
    BEs = pd.read_csv(path_bes,sep='\t',index_col=0)
    BEs.columns = BEs.columns[[c.startswith("be_") and c[3:].isdigit() for c in BEs.columns]]
    # isolate informative reads
    informative_BEs = BEs[BEs.sum(axis=1) > 0]
    combinations = set([tuple(informative_BEs.loc[i,:]) for i in informative_BEs.index])
    combo_counter = {c:[0,list()] for c in combinations}
    for i in informative_BEs.index:
        combo_counter[tuple(informative_BEs.loc[i,:])][0] += 1
        combo_counter[tuple(informative_BEs.loc[i,:])][1].append(i)
    if prefix == "":
        prefix = path_bes.stem
    # write read groups to files and return their names
    filenames = []
    # calculate the groups to consider
    selected_groups = []
    for key in combo_counter:
        if combo_counter[key][0] >= minimum_reads:
            selected_groups.append(key)
    for i,key in enumerate(selected_groups):
        wpath = out_dir / (prefix+'.'+str(i).zfill(len(str(len(selected_groups))))+'.txt')
        wpath.parent.mkdir(exist_ok=True,parents=True)
        with open(wpath,'w') as f:
            for line in combo_counter[key][1]:
                print(line,file=f)
        filenames.append(wpath)
    return filenames

def align_assembly(reference, fastain, bamout, threads):
    cmd_align = shlex.split(f"minimap2 -a -x asm10 -t {threads} -r95,95 -z95,95 {reference} {fastain}")
    cmd_compress = shlex.split("samtools view -b")
    cmd_sort = shlex.split(f"samtools sort -o {bamout}")
    cmd_index = shlex.split(f"samtools index {bamout}")
    p_align = subprocess.Popen(cmd_align, stdout=subprocess.PIPE)
    p_compress= subprocess.Popen(cmd_compress,stdin=p_align.stdout, stdout=subprocess.PIPE)
    p_align.stdout.close()
    p_sort= subprocess.Popen(cmd_sort,stdin=p_compress.stdout)
    p_compress.stdout.close()
    logger.info(f"aligning reads, compressing and sorting alignment..")
    logger.info(f"calling: {' '.join(cmd_align)}")
    p_sort.communicate()
    subprocess.call(cmd_index)

def consensus_from_given_read_names(readnames,reads,reference,fastaout,fastaout_reversed,bamout,bamout_reversed,threads,fasta_name,technology):
    mats = {'ont':"promethion.mat",
            'pb':'sequel-II-CLR.mat'}
    if not technology in mats.keys():
        logger.error(f"technology must be chosen from: {mats.keys()}")
        return
    cmd_whitelines = shlex.split(f"sed -i '/^$/d' {readnames}")
    logger.info(f"removing empty lines from readnames: {''.join(cmd_whitelines)}")
    p_whitelines = subprocess.check_call(cmd_whitelines)
    # zcat reads and grep them, write to tmpfile (readnames,reads,tmpfile)
    cmd_read = shlex.split(f"zcat -f {reads}")
    cmd_grep = shlex.split(f"grep --no-group-separator -A 3 -f {readnames}") # pipe from read
    selected_reads = tempfile.NamedTemporaryFile(suffix='.fastq') # stdout of grep to this file
    p_read = subprocess.Popen(cmd_read, stdout=subprocess.PIPE)
    p_grep = subprocess.Popen(cmd_grep, stdin=p_read.stdout,stdout=selected_reads)
    logger.info(f"selecting reads from {readnames} in {reads}..")
    p_grep.communicate()
    # consensus from tmpfile to fastaout (tmpfile,fastaout)

    # subprocess.check_call(shlex.split(f"cp {selected_reads.name} /fast/work/users/mayv_c/development/atlatl/tests/data/assembly/selected_reads.fastq"))

    matfile = os.path.join(os.path.dirname(__file__), "lamassemble", mats[technology])
    cmd_smbl = shlex.split(f"lamassemble -P {threads} -n {fasta_name} {matfile} {selected_reads.name}")
    with open(fastaout, 'wb') as f:
        logger.info(f"calling: {' '.join(cmd_smbl)}")
        subprocess.call(cmd_smbl,stdout=f)
    align_assembly(reference, fastaout, bamout, threads)
    with open(fastaout_reversed,'wb') as f:
        cmd_reverse = shlex.split(f"seqtk seq -r {fastaout}")
        subprocess.check_call(cmd_reverse,stdout=f)
    align_assembly(reference, fastaout_reversed, bamout_reversed, threads)

# =============================================================================
#  visualization

# overlapt of triples (start,end,height)
def overlap(a, b):
    return (a[2] == b[2]) * max(0, min(a[1], b[1]) - max(a[0], b[0]))

# a collection of triples (start,end,height)
def get_heights(tiles):
    L = np.array([(*t,0) for t in tiles])
    # build conflict graph
    G = pd.DataFrame(np.zeros((len(tiles),len(tiles))))
    while True:
        for a in G.index:
            for b in G.columns:
                if a != b:
                    G.loc[a,b] = overlap(L[a], L[b])
        conflicts = G.astype(bool).sum()
        if G.astype(bool).sum().sum() == 0:
            break
        lift_tile = conflicts[conflicts == conflicts.max()].index[-1]
        L[lift_tile][2] += 1
    return np.array([l[2] for l in L])

# function to return start and end of alignment on query
def start_end_query_on_read(query,reverse=False):
    hard_clipped_start = query.cigartuples[0][1] * (query.cigartuples[0][0] == 5)
    s = hard_clipped_start + query.query_alignment_start
    e = s + query.query_alignment_length
    if query.is_reverse != reverse:
        l = query.infer_read_length()
        e,s = l-s, l-e
    return s,e

def visualize_assembly(alignments_path:pathlib.PosixPath,annotations_path:pathlib.PosixPath,outpath:pathlib.PosixPath,chrs:list,thickness:int=5):
    aln_path   = pathlib.Path(alignments_path)
    # load all alignments to dict, where each key corresponds to a chr
    alignments = defaultdict(list)
    # pick only the first read in bam file
    samplename = ""
    for a in pysam.AlignmentFile(aln_path, "rb"):
        if not samplename:
            samplename = a.query_name
        if a.query_name == samplename:
            alignments[a.reference_name].append(a)
    if len(chrs) > 0:
        chromosomes = list(set(chrs).intersection(set(np.unique(np.array(list(alignments.keys()))))))
        if len(chromosomes) == 0:
            logger.warn("No chromosomes selected to print. check spelling or leave chromosmes blank")
            return
    else:
        chromosomes = np.unique(np.array(list(alignments.keys())))
    fig = make_subplots(rows=len(chromosomes), cols=1)
    for chri,chr in enumerate(chromosomes):
        alns = np.array(list(alignments[chr]))
        refstart, refend = min([a.reference_start for a in alns]), max([a.reference_end for a in alns])
        readstart = min([start_end_query_on_read(a)[0] for a in alns])
        # select and trim annotations to fit them into the selected interval
        tmp_region_file = tempfile.NamedTemporaryFile(suffix='.bed')
        tmp_annot_file = tempfile.NamedTemporaryFile(suffix='.bed')
        with open(tmp_region_file.name,'w') as f:
            print('\t'.join(list(map(str,[chr,refstart,refend]))),file=f)
        cmd_intersect = shlex.split(f"bedtools intersect -nonamecheck -a {annotations_path} -b {tmp_region_file.name}")
        with open(tmp_annot_file.name,'w') as g:
            subprocess.check_call(cmd_intersect,stdout=g)
            if pathlib.Path(tmp_annot_file.name).stat().st_size:
                features = pd.read_csv(tmp_annot_file.name,sep='\t',header=None)
            else:
                logger.warn(f"No annotations selected for {chr} in interval {str(refstart)},{str(refend)}")
                continue
        features.columns = ["chr","start","end","name"]
        features.loc[:,['start','end']] -= refstart
        features = features.sort_values(by='start')
        features.index = range(features.shape[0])
        a = alns[0]
        sum_rev = sum([int(a.is_reverse) * (a.reference_end - a.reference_start) for a in alns])
        sum_fwd = sum([int(not a.is_reverse) * (a.reference_end - a.reference_start) for a in alns])
        alns_reversed = sum_fwd >= sum_rev
        # define colors
        red,blue,dgrey,lgrey = 'rgb(201, 55, 44)','rgb(44, 96, 201)','rgb(100,100,100)','rgb(200,200,200)'
        if alns_reversed:
            red,blue = blue,red
        annotation_color = 'rgb(150,150,220)'
        # calc height for alignments on read
        alns_read_heights = get_heights([start_end_query_on_read(a) for a in alns])
        alns_ref_heights = get_heights([(a.reference_start,a.reference_end) for a in alns])
        space_read_to_ref = thickness * (alns_read_heights.max() + alns_ref_heights.max() + 10)
        # --- draw read --- #
        fig.add_trace(
            go.Scatter(
                x=np.array([0, 0, a.infer_read_length(),a.infer_read_length(),0]) - readstart,
                y=np.array([0, thickness, thickness, 0, 0]) + space_read_to_ref,
                fill='toself',
                fillcolor= lgrey,
                hoveron = 'fills', # select where hover is active
                line_color='grey',
                text=f"Read, len:{'{:,}'.format(a.infer_read_length())}",
                hoverinfo = 'text',
                showlegend=False),row=chri+1, col=1)
        # --- draw ref --- #
        fig.add_trace(
            go.Scatter(
                x=np.array([0, 0, refend - refstart, refend - refstart, 0]),
                y=np.array([0, thickness, thickness, 0, 0]),
                fill='toself',
                fillcolor= dgrey,
                hoveron = 'fills', # select where hover is active
                line_color='grey',
                text=f"Ref= {chr}:{'{:,}'.format(refstart)}-{'{:,}'.format(refend)}",
                hoverinfo = 'text',
                showlegend=False),row=chri+1, col=1)
        # --- draw alignments --- #
        for i,a in enumerate(alns):
            col = blue if a.is_reverse else red
            h = alns_read_heights[i]
        # --- alignments on read --- #
            qstart,qend = start_end_query_on_read(a,reverse=False)
            x_aread=np.array([
                        qstart,
                        qstart,
                        qend,
                        qend,
                        qstart]) - readstart
            y_aread=np.array([0, thickness, thickness, 0, 0]) + space_read_to_ref - thickness - thickness*h
            aln_text = f"alignment: {'{:,}'.format(qstart)} {'{:,}'.format(qend)}; ref: {'{:,}'.format(a.reference_start)} {'{:,}'.format(a.reference_end)}"
            #positions[a] = (x[[1,2]],y[0]) # use y[0] to connect to bottom
            fig.add_trace(
                go.Scatter(
                    x = x_aread,
                    y = y_aread,
                    fill='toself',
                    fillcolor=col,
                    hoveron = 'fills', # select where hover is active
                    line_color=col,
                    text=aln_text,
                    hoverinfo = 'text',
                    name=f"alignment {i+1}",
                    showlegend=False
                    ),row=chri+1, col=1)
        # --- alignments on reference --- #
            rstart,rend = a.reference_start, a.reference_end
            h = alns_ref_heights[i]
            x_aref=np.array([
                        rstart,
                        rstart,
                        rend,
                        rend,
                        rstart]) - refstart
            y_aref=np.array([0, thickness, thickness, 0, 0]) + thickness + thickness*h
            fig.add_trace(
                go.Scatter(
                    x = x_aref,
                    y = y_aref,
                    fill='toself',
                    fillcolor=col,
                    hoveron = 'fills', # select where hover is active
                    line_color=col,
                    text=aln_text,
                    hoverinfo = 'text',
                    name=f"alignment {i+1}",
                    showlegend=False
                    ),row=chri+1, col=1)
            # --- draw ref to read lines --- #
            fig.add_trace(
                go.Scatter(
                    x = np.array([x_aref[1],x_aread[0],x_aread[3],x_aref[2],x_aref[1]]),
                    y = np.array([y_aref[1],y_aread[0],y_aread[3],y_aref[2],y_aref[1]]),
                    fill='toself',
                    hoveron = 'fills', # select where hover is active
                    line_color=col,
                    text=aln_text,
                    hoverinfo = 'text',
                    name=f"alignment {i+1}",
                    showlegend=False
                    ),row=chri+1, col=1)
        # --- draw annotations --- #
        feature_heights = get_heights([tuple(features.loc[j,['start','end']]) for j in features.index])
        for i in features.index:
            fs,fe,ft = features.loc[i,['start','end','name']]
            fh = feature_heights[i]
            fig.add_trace(
                go.Scatter(
                    x=np.array([fs, fs, fe, fe, fs]),
                    y=np.array([0, thickness, thickness, 0, 0]) - fh * thickness,
                    fill='toself',
                    fillcolor=annotation_color,
                    hoveron = 'fills', # select where hover is active
                    line_color='grey',
                    text=f"Ref= {chr} {ft}:{'{:,}'.format(refstart+fs)}-{'{:,}'.format(refstart+fe)}", # maybe add coordinates
                    hoverinfo = 'text',
                    showlegend=False),row=chri+1, col=1)
        # save fig and return
    fig.update_layout(showlegend=False)
    pathlib.Path(outpath).parent.mkdir(parents=True,exist_ok=True)
    if outpath.suffix == '.html':
        fig.write_html(outpath)
    else:
        fig.write_image(outpath)
    #return fig

# =============================================================================
#  selective_assembly & visualization

def group_assemble_and_visualize(
                            breakends:str,
                            chrs:list,
                            annotations:str,
                            reads:str,
                            reference:str,
                            outdir:str,
                            threads:int=3,
                            prefix:str="",
                            minimum_reads:int=5,
                            fasta_name:str="",
                            technology:str='ont',
                            thickness:int=5,
                            img_format='html'):
    """
    Generates visualizations of all informative read groups found in a table of breakends.
    """
    fmts = ['html','png','pdf']
    if not img_format in fmts:
        logger.error("Must choose file format from ")
    filenames = readgroups_to_files(breakends,outdir,prefix,minimum_reads)
    for i,readgroup in enumerate(filenames):
        strkey = str(i).zfill(len(str(len(filenames))))
        used_prefix = strkey if prefix == "" else prefix + '.' + strkey
        outfiles_prefix = pathlib.Path(outdir) / used_prefix
        assemble_and_visualize(
                            readgroup=readgroup,
                            annotations=annotations,
                            reads=reads,
                            reference=reference,
                            outfiles_prefix=outfiles_prefix,
                            threads=threads,
                            chrs=chrs,
                            fasta_name=fasta_name,
                            technology=technology,
                            thickness=thickness)

def assemble_and_visualize(readgroup:str,
                            annotations:str,
                            reads:str,
                            reference:str,
                            outfiles_prefix:str,
                            chrs:list,
                            threads:int=3,
                            fasta_name:str="",
                            technology:str='ont',
                            image_format:str='html',
                            thickness:int=5):
    """
    Generates visualizations of assemblies of given read names ('readgroup').

    readgroup: .txt with one read name per line.
    """
    pn = pathlib.Path(outfiles_prefix).parent
    fn = pathlib.Path(outfiles_prefix).stem
    consensus_from_given_read_names(readgroup,
                                    reads,
                                    reference,
                                    pn / '.'.join([fn,"fwd","fasta"]),
                                    pn / '.'.join([fn,"rvs","fasta"]),
                                    pn / '.'.join([fn,"fwd","bam"]),
                                    pn / '.'.join([fn,"rvs","bam"]),
                                    threads,
                                    fasta_name,
                                    technology)
    visualize_assembly(pn / '.'.join([fn,"fwd","bam"]),
                        annotations,
                        pn / '.'.join([fn,"fwd",image_format]),
                        chrs,
                        thickness=thickness)
    visualize_assembly(pn / '.'.join([fn,"rvs","bam"]),
                        annotations,
                        pn / '.'.join([fn,"rvs",image_format]),
                        chrs,
                        thickness=thickness)
