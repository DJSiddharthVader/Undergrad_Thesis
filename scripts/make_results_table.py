import os
import sys
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
import scipy.stats as sst
from functools import reduce

missing_val = -1
marko_filename = 'markophylo_results.txt'
marko_file_errors = ['[1] "plotting tree"','"Error']
report_filename = 'markophylo_results.txt'

def hasRate(lines,cat):
    boollist = [True if cat in l else False for l in lines ]
    return reduce(lambda x,y: x or y, boollist)

def getRates(lines,pattern,has_c,has_nc):
    if pattern not in lines:
        return missing_val,missing_val,missing_val
    rates = lines[lines.index(pattern)+4].split(' ')[1:]
    if set(marko_file_errors).intersection(set(rates)):
        return missing_val,missing_val,missing_val
    rates = [float(x) for x in rates if x != '']
    if has_c and has_nc:
        crate,ncrate,intrate = rates
    elif has_c:
        crate,intrate = rates
        ncrate = missing_val
    elif has_nc:
        ncrate,intrate = rates
        crate = missing_val
    else:
        raise ValueError('Only has internal partition?')
    return crate,ncrate,intrate

def parseMarko(basedir,genus):
    filepath = os.path.join(basedir,genus,marko_filename)
    if not(os.path.exists(filepath)):
        return None
    lines = [x.strip() for x in open(filepath).readlines()]
    if len(lines) == 0:
        return None
    has_c = hasRate(lines,'$crispr')
    has_nc = hasRate(lines,'$non_crispr')
    crate,ncrate,intrate = getRates(lines,'$rates',has_c,has_nc)
    csem,ncsem,intsem = getRates(lines,'$se$rates',has_c,has_nc)
    return {'genus':genus,
            'crispr_indel_rate':crate,
            'non-crispr_indel_rate':ncrate,
            'internal_indel_rate':intrate,
            'crispr_indel_sem':csem,
            'non-crispr_indel_sem':ncsem,
            'internal_indel_sem':intsem}

def makeMarkoDF(basedir,inputs):
    marko_data=[parseMarko(basedir,genus) for genus in tqdm(inputs)]
    marko_df = pd.DataFrame(list(filter(None,marko_data)))
    return marko_df

def reportToDF(report_file):
    if not(os.path.exists(report_file)):
        return pd.DataFrame()
    report = json.load(open(report_file))
    report['modularity.mean'] = np.mean(report['modularity'])
    report['modularity.sem'] = sst.sem(report['modularity'])
    report['assortativity.mean'] = np.mean(report['assortativity'])
    report['assortativity.sem'] = sst.sem(report['assortativity'])
    report.pop('modularity')
    report.pop('assortativity')
    return pd.io.json.json_normalize(report)

def makeReportDF(basedir,genera,smode):
    filepaths = [os.path.join(basedir,genus,'{}_{}_stat_report.json'.format(genus,smode)) for genus in genera]
    dfs = [reportToDF(report_file) for report_file in tqdm(filepaths)]
    report_df = pd.concat([x for x in dfs if not(x.empty)])
    return report_df

def main(basedir,marko_outpath,report_outpath,smode):
    inputs = [genus for genus in os.listdir(basedir) \
              if os.path.isdir(os.path.join(basedir,genus))]
    marko_df = makeMarkoDF(basedir,inputs)
    marko_df.to_csv(marko_outpath,index=False)
    report_df = makeReportDF(basedir,inputs,smode)
    report_df.to_csv(report_outpath,index=False)
    combined_df = marko_df.merge(report_df,how='outer',on='genus')
    combined_df.to_csv('all_{}_results.csv'.format(smode),index=False)
    return None

if __name__ == '__main__':
    if len(sys.argv) > 1:
        basedir = sys.argv[1]
    else:
        basedir =  os.getcwd()
    if len(sys.argv) > 2:
        marko_outpath = sys.argv[2]
    else:
        marko_outpath = 'all_markophylo_results.csv'
    if len(sys.argv) > 3:
        smode = sys.argv[3]
    else:
        smode = 'WGS'
    if len(sys.argv) > 4:
        report_outpath= sys.argv[4]
    else:
        report_outpath = 'all_{}_report_results.csv'.format(smode)
    print(basedir,marko_outpath,smode,report_outpath)
    main(basedir,marko_outpath,report_outpath,smode)

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
#def mainSingle(basedir,outpath):
#    genus = os.path.basename(basedir)
#    inpath = os.path.join(basedir,'markophylo_results.txt')
#    outdata = parse2(inpath,genus)
#    json.dump(outdata,open(outpath,'w'))

