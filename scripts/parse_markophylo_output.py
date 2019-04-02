import os
import sys
import json

def parse(filepath,genus):
    lines = [x.strip() for x in open(filepath).readlines()]
    ratelist = lines[lines.index('$rates')+4].split(' ')[1:]
    ratelist = [x for x in ratelist if x != '']
    crate,ncrate = [float(x) for x in ratelist]
    semlist = lines[lines.index('$se$rates')+4].split(' ')[1:]
    semlist = [x for x in semlist if x != '']
    csem,ncsem = [float(x) for x in semlist]
    return {'genus':genus,
            'crispr_indel_rate':crate,
            'non-crispr_indel_rate':ncrate,
            'crispr_indel_sem':csem,
            'non-crispr_indel_sem':ncsem}

def parse2(filepath,genus):
    lines = [x.strip() for x in open(filepath).readlines()]
    #rates
    ratelist = lines[lines.index('$rates')+4].split(' ')[1:]
    ratelist = [x for x in ratelist if x != '']
    rates = [float(x) for x in ratelist]
    if len(rates) == 2:
        crate,ncrate = rates
    else:
        crate, ncrate = 0,rates[0]
    #std errs
    semlist = lines[lines.index('$se$rates')+4].split(' ')[1:]
    semlist = [x for x in semlist if x != '']
    if len(semlist) == 2:
        csem,ncsem = semlist
    else:
        csem,ncsem = 0,semlist[0]
    return {'genus':genus,
            'crispr_indel_rate':crate,
            'non-crispr_indel_rate':ncrate,
            'crispr_indel_sem':csem,
            'non-crispr_indel_sem':ncsem}

if __name__ == '__main__':
    if len(sys.argv) > 2:
        basedir = sys.argv[1]
        genus = os.path.basename(basedir)
        inpath = os.path.join(os.getcwd(),basedir,'markophylo_results.txt')
        outd = parse2(inpath,genus)
        outpath = os.path.join(basedir,'{}_markophylo_results.json'.format(genus))
        json.dump(outd,open(outpath,'w'))
    else:
        basedir = sys.argv[1]
        genus = os.path.basename(basedir)
        inpath = os.path.join(os.getcwd(),basedir,'markophylo_results.txt')
        outd = parse(inpath,genus)
        outpath = os.path.join(basedir,'{}_markophylo_results.json'.format(genus))
        json.dump(outd,open(outpath,'w'))

#for dir in `ls | grep -v -E "missing*|notEn*|\.|noCRISPR*" `; do python ../../scripts/parse_markophylo_output.py "$dir"; done
#find genus_data -maxdepth 2 -name "*markophylo_results.json" -exec cp {} marko_reports/ \;
