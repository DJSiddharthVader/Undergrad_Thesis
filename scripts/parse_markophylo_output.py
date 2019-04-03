import os
import sys
import json
from functools import reduce

def hascat(lines,cat):
    boollist = [True if cat in l else False for l in lines ]
    return reduce(lambda x,y:x or y,boollist)

def parse2(filepath,genus):
    lines = [x.strip() for x in open(filepath).readlines()]
    #parse rates
    rates = lines[lines.index('$rates')+4].split(' ')[1:]
    rates = [float(x) for x in rates if x != '']
    #parse categories
    has_c = hascat(lines,'$crispr')
    has_nc = hascat(lines,'$non_crispr')
    if has_c and has_nc:
        crate,ncrate,intrate = rates
    elif has_c:
        crate,intrate = rates
        ncrate = 0
    elif has_nc:
        ncrate,intrate = rates
        crate = 0
    else:
        raise ValueError('Only has internal partition?')
    #parse std errs
    try:
        sems = lines[lines.index('$se$rates')+4].split(' ')[1:]
        sems = [float(x) for x in sems if x != '']
        if has_c and has_nc:
            csem,ncsem,intsem = sems
        elif has_c:
            csem,intsem = sems
            ncsem = 0
        elif has_nc:
            ncsem,intsem = sems
            csem = 0
        else:
            raise ValueError('Only has internal partition?')
        return {'genus':genus,
                'crispr_indel_rate':crate,
                'non-crispr_indel_rate':ncrate,
                'internal_indel_rate':intrate,
                'crispr_indel_sem':csem,
                'non-crispr_indel_sem':ncsem,
                'internal_indel_sem':intsem}
    except ValueError:
        return {'genus':genus,
                'crispr_indel_rate':crate,
                'non-crispr_indel_rate':ncrate,
                'internal_indel_rate':intrate,
                'crispr_indel_sem':-1,
                'non-crispr_indel_sem':-1,
                'internal_indel_sem':-1}

if __name__ == '__main__':
    basedir = sys.argv[1]
    genus = os.path.basename(basedir)
    inpath = os.path.join(os.getcwd(),basedir,'markophylo_results.txt')
    outd = parse2(inpath,genus)
    print(outd)
    outpath = os.path.join(basedir,'{}_markophylo_results.json'.format(genus))
    json.dump(outd,open(outpath,'w'))

#for dir in `ls | grep -v -E "missing*|notEn*|\.|noCRISPR*" `; do python ../../scripts/parse_markophylo_output.py "$dir"; done
#find genus_data -maxdepth 2 -name "*markophylo_results.json" -exec cp {} marko_reports/ \;

#for g in `ls`; do cd "$g"; Rscript ~/thesis_SidReed/scripts/markophylo_indel_estiamtes.R; cd -; python ~/thesis_SidReed/scripts/parse_markophylo_output.py "$g"; done
#run on all genera dir in a dir

#DEPRECIATED
##rates
#ratelist = lines[lines.index('$rates')+4].split(' ')[1:]
#ratelist = [x for x in ratelist if x != '']
#rates = [float(x) for x in ratelist]
#if len(rates) != 3:
#    none = [0]
#    none.extend(rates)
#    rates = none
#crate,ncrate,intrate = rates
#try:
#    #std errs
#    semlist = lines[lines.index('$se$rates')+4].split(' ')[1:]
#    semlist = [x for x in semlist if x != '']
#    if len(semlist) != 3:
#        none = [0]
#        none.extend(semlist)
#        semlist = none
#    csem,ncsem,intsem = semlist
#    return {'genus':genus,
#            'crispr_indel_rate':crate,
#            'non-crispr_indel_rate':ncrate,
#            'internal_indel_rate':intrate,
#            'crispr_indel_sem':csem,
#            'non-crispr_indel_sem':ncsem,
#            'internal_indel_sem':ncsem}
#except ValueError:
#    return {'genus':genus,
#            'crispr_indel_rate':crate,
#            'non-crispr_indel_rate':ncrate,
#            'internal_indel_rate':intrate}
