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
    if len(rates) != 3:
        rates = [0].extend(rates)
        print(rates)
    crate,ncrate,intrate = rates
    try:
        #std errs
        semlist = lines[lines.index('$se$rates')+4].split(' ')[1:]
        semlist = [x for x in semlist if x != '']
        if len(semlist) != 3:
            semlist = [0].extend(semlist)
            print(semlist)
        csem,ncsem,intsem = semlist
        return {'genus':genus,
                'crispr_indel_rate':crate,
                'non-crispr_indel_rate':ncrate,
                'internal_indel_rate':intrate,
                'crispr_indel_sem':csem,
                'non-crispr_indel_sem':ncsem,
                'internal_indel_sem':ncsem}
    except ValueError:
        return {'genus':genus,
                'crispr_indel_rate':crate,
                'non-crispr_indel_rate':ncrate,
                'internal_indel_rate':intrate}


if __name__ == '__main__':
    basedir = sys.argv[1]
    genus = os.path.basename(basedir)
    inpath = os.path.join(os.getcwd(),basedir,'markophylo_results.txt')
    outd = parse2(inpath,genus)
    outpath = os.path.join(basedir,'{}_markophylo_results.json'.format(genus))
    json.dump(outd,open(outpath,'w'))

#for dir in `ls | grep -v -E "missing*|notEn*|\.|noCRISPR*" `; do python ../../scripts/parse_markophylo_output.py "$dir"; done
#find genus_data -maxdepth 2 -name "*markophylo_results.json" -exec cp {} marko_reports/ \;
