import os
import sys
import json

def parse(filepath):
    lines = [x.strip() for x in open(filepath).readlines()]
    ratelist = lines[lines.index('$rates')+4].split(' ')[1:]
    ratelist = [x for x in ratelist if x != '']
    crate,ncrate = [float(x) for x in ratelist]
    semlist = lines[lines.index('$se$rates')+4].split(' ')[1:]
    semlist = [x for x in semlist if x != '']
    csem,ncsem = [float(x) for x in semlist]
    return {'crispr_indel_rate':crate,
            'non-crispr_indel_rate':ncrate,
            'crispr_indel_sem':csem,
            'non-crispr_indel_sem':ncsem}

if __name__ == '__main__':
    basedir = sys.argv[1]
    genus = os.path.basename(basedir)
    outd = parse(os.path.join(basedir,'markophylo_results.txt'))
    outpath = os.path.join(basedir,'{}_markophylo_results.json'.format(genus))
    json.dump(outd,open(outpath,'w'))

#for dir in `ls | grep -v -E "missing*|notEn*|\.|noCRISPR*" `; do python ../../scripts/parse_markophylo_output.py "$dir"; done
#find genus_data -maxdepth 2 -name "*markophylo_results.json" -exec cp {} marko_reports/ \;
