import re
import sys
import json
import requests
import pandas as pd
from bs4 import BeautifulSoup
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

base = 'https://www.ncbi.nlm.nih.gov/assembly/'

def geturllist(reffile):
    gcas = list(pd.DataFrame(json.load(open(reffile)))['GCA'])
    return gcas

def getGCANameInfo(gca):
    gcatrim = '_'.join(gca.split('_')[:2])
    url = base + gcatrim
    data = requests.get(url)
    soup = BeautifulSoup(data.text,'html.parser')
    categories =  [x.text[:-2] for x in soup.find_all('dt')]
    values = [x.text for x in soup.find_all('dd')]
    ddata = {c:v for c,v in zip(categories,values)}
    ddata.update({'GCA':gca})
    return ddata

def main(reffile,outf,processes):
    pool = ThreadPool(processes)
    urllist = geturllist(reffile)
    results = pool.map(getGCANameInfo, geturllist(reffile))
    json.dump(results,open(outf,'w'))

if __name__ == '__main__':
    reffile = '/home/sid/thesis_SidReed/CRISPRone_files/mp_CRISPRoneAnnotations_01_11_18SR.json'
    if len(sys.argv) == 3:
        processes = int(sys.argv[2])
    else:
        processes = 8
    main(reffile,sys.argv[1],processes)
