import os
import sys
import glob
from tqdm import tqdm
import subprocess as sp
from Bio.Nexus import Nexus
from Bio import SeqIO, AlignIO, Alphabet
from multiprocessing.dummy import Pool as ThreadPool

"""
Builds a species tree using a core genome alignment
Steps
    1. Do core genome alignment of all taxa with parsnp
    2. Split .xmfa file into seperate .fna, 1 per aligned region, excluding all non-informational alignments*
    3. Convert each .fna file into a nexus file with Biopython AlignIO
    4. Concatenate all the nexus files into a signle nexus alignment file
    5. Run MrBayes on the single alignment file to build a species tree

* XMFA files contain a bunch of seperate alignments, each dilineated with a =. non-informational alignments meaning that all sequences in the alignment are the same, so they provide no information during the tree building process.

Input data is the files containing the genomes of the taxa you want to build a species tree for (_genomic.fna files from NCBI)
File chain
    1. ./genome/*_genomic.fna
    2. ./species_tree_WGS/WGS.xmfa
    3. ./species_tree_WGS/fasta/WGS_1.aln_fna
                               /WGS_2.aln_fna
                               ...
    4. ./species_tree_WGS/nexus/WGS_1.nex
                               /WGS_2.nex
                               ...
    5. ./species_tree_WGS/WGS.nex
    6. ./species_tree_WGS/WGS_species_tree.con.tree
                          + other MrBayes files

Output is all the intermediary alignment files and the output of MrBayes (i.e the species tree WGS_species_tree.con.tree)
"""
# Directories
base_dir = 'species_tree_WGS'
fasta_files = '{}/fasta_alignments'.format(base_dir)
nexus_files = '{}/fasta_alignments'.format(base_dir)

# Alignment Files
mauve_path = '/home/sid/thesis_SidReed/apps/mauve_snapshot_2015-02-13/linux-x64/progressiveMauve'
parsnp_path = '/usr/local/bin/parsnp'
mbscriptname = 'WGS_mrbayes_script.txt'
logfilename = 'WGS_mrbayes_log.txt'


def getBase(filename, splitchar='.'):
    base = splitchar.join(os.path.basename(filename).split(splitchar)[0:-1])
    return base


def multiAlignGenomes(genomde_dir, threads, method='parsnp'):
    # align core genomes of all taxa using parsnp, creates a multi-alingment including all taxas used to build a species tree
    # uses more information thatn 16S genes only and works even when there are
    # not enough annotated 16S genes
    outfile = '{}/WGS.xmfa'.format(base_dir)
    if os.path.isfile(outfile):
        return outfile
    if method == 'parsnp':
        # -c uses all files in the target dir
        # -r ! randomly choose a reference fasta
        # multi core genome alignment for all members of a genus
        parsnpcmd = '{} -c -r ! -p {} -d {} -o {} '.format(
            parsnp_path, threads, genome_dir, base_dir)
        os.system(parsnpcmd)
        os.system('mv {}/parsnp.xmfa {}'.format(base_dir, outfile))
        return outfile
    elif method == 'mauve':
        cmd = '{} --output={} {}/*.fna'.format(mauve_path, outfile, genome_dir)
        os.system(cmd)
        return outfile
    else:
        raise ValueError('invalid alinger choose {parsnp|mauve}')


def getHeaders(xmfa):
    # get list of input genome headers (names) from xmfa file
    cmd = ['grep "^##SequenceHeader"', xmfa, '| cut  -d" "  -f2']
    headers = sp.getoutput(' '.join(cmd)).split('\n')
    headers = [x.strip('\'\"') for x in headers]
    return headers


def trimXMFA(xmfa):
    # remove all comment, = lines from xmfa so it looks like a fasta
    temp_aln = getBase(xmfa)
    temp_aln = '{}/{}.aln_fna'.format(base_dir, temp_aln)
    os.system('grep "^[>agctAGCT-]" {} >| {}'.format(xmfa, temp_aln))
    return temp_aln


def convertXMFAToFastas(xmfa):
    # take xmfa and convert each section into a seperate fasta aln file
    # ignore biopython comparison warning
    # warnings.filterwarnings("ignore", category=BiopythonWarning)
    headers = getHeaders(xmfa)
    trimmed_xmfa = trimXMFA(xmfa)
    records = list(SeqIO.parse(open(trimmed_xmfa), 'fasta'))
    n_headers = len(headers)
    for ri in range(0, len(records) // n_headers):
        base_name = '.'.join(os.path.basename(xmfa).split('.')[0:-1])
        output_file = '{}/fasta/{}_{}.aln_fna'.format(base_dir, base_name, ri)
        if os.path.isfile(output_file):
            next
        aln = records[n_headers * ri:n_headers * (ri + 1)]
        seqs = [a.seq.upper() for a in aln]
        if len(
                set(seqs)) == 1:  # if all seqs are the same, no point in including this genome section
            for i in range(0, n_headers):
                aln[i].id = headers[i].strip('\'\">')
                aln[i].seq = aln[i].seq.upper()
            with open(output_file, 'w') as output:
                SeqIO.write(aln, output, 'fasta')
    return None


def convertFastaToNexus(fasta_file):
    # convert a fasta aln file to a nexus file
    fasta_base = getBase(fasta_file)
    nexus_file = '{}/nexus/{}.nex'.format(base_dir, fasta_base)
    AlignIO.convert(fasta_file, 'fasta',
                    nexus_file, 'nexus',
                    alphabet=Alphabet.generic_dna)
    return None


def concatNexusAlignments(processes):
    # take the list of fasta alingments and convert each to a nexus file and
    # concat all the nexus files into 1 alingment
    pool = ThreadPool(processes)
    already_done = [x.split('.')[0]
                    for x in os.listdir('{}/nexus'.format(base_dir))]
    fastas = ['{}/fasta/{}'.format(base_dir, file) for
              file in os.listdir('{}/fasta'.format(base_dir))
              if file.split('.')[0] not in already_done]
    list(tqdm(pool.imap(convertFastaToNexus, fastas),
              total=len(fastas), desc='Fastas to Nexus...'))
    combined_nexus = '{}/WGS.nex'.format(base_dir)
    if os.path.isfile(combined_nexus):
        return combined_nexus
    nexus = ['{}/nexus/{}'.format(base_dir, file) for
             file in os.listdir('{}/nexus'.format(base_dir))]
    nexus = [(filename, Nexus.Nexus(filename)) for filename in nexus]
    combined = Nexus.combine(nexus)
    combined.write_nexus_data(filename=open(combined_nexus, 'w'))
    return combined_nexus


def buildSpeciesTree(nexus_alignment):
    # run mrbayes to produce a species tree from a parsnp core genome alignment
    mbf =\
        """set autoclose=yes nowarn=yes
execute {}
lset nst=6 rates=gamma
mcmc ngen=10000 savebrlens=yes file=WGS_species_tree
sump burnin=250
sumt burnin=250
quit""".format(nexus_alignment)
    open(mbscriptname, 'w').write(mbf)  # create mr bayes script
    with open(mbscriptname, 'rb', 0) as mbscript, open(logfilename, 'wb', 0) as logfile:
        logtxt = sp.run(['mb'],
                        stdin=mbscript,
                        stdout=logfile,
                        check=True)
    for fp in glob.glob('./WGS*'):
        os.rename(fp, os.path.join(base_dir, fp))
    return None


def main(genus, genome_dir, threads):
    os.system('mkdir -p {}'.format(base_dir))
    os.system('mkdir -p {}/fasta'.format(base_dir))
    os.system('mkdir -p {}/nexus'.format(base_dir))
    xmfa = multiAlignGenomes(genome_dir, threads)
    convertXMFAToFastas(xmfa)
    nexus = concatNexusAlignments(threads)
    if os.path.isfile('{}/WGS_species_tree.con.tree'.format(base_dir)):
        return None
    buildSpeciesTree(nexus)
    return None


if __name__ == '__main__':
    if len(sys.argv) > 1:
        genus = sys.argv[1]
    else:
        genus = os.path.basename(os.getcwd())
    if len(sys.argv) > 2:
        genome_dir = sys.argv[2]
    else:
        genome_dir = './genome'
    if len(sys.argv) > 3:
        threads = int(sys.argv[3])
    else:
        threads = 1
    main(genus, genome_dir, threads)
