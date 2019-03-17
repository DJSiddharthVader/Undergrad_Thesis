import re
import sys
import json
import requests
import pandas as pd
from tqdm import tqdm
from bs4 import BeautifulSoup
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool


#def updateVersion(gca,newv):
#    try:
#        gcan,asm = gca.split('_')[1:]
#    except ValueError:
#        print(gca)
#        print(gca.split('_'))
#        sys.exit()
#    gcan,vn = gcan.split('.')
#    gcan = '.'.join([gcan,str(newv)])
#    ngca = '_'.join(['GCA',gcan,asm])
#    return(ngca)

#def makeSoup(gca):
#    gcatrim = '_'.join(['GCF',gca.split('_')[1]])
#    url = baseurl + gcatrim
#    data = requests.get(url)
#    soup = BeautifulSoup(data.text,'html.parser')
#    try:
#        soup = makeSoup(gca)
#        table = soup.find_all('table')[3]
#        correctVersion = True
#    except IndexError:
#        version += 1
#        print('updating version')
#        print(gca)
#    return(soup)

def getgcalist(reffile):
    gcas = list(pd.DataFrame(json.load(open(reffile)))['GCA'])
    gcas = [x for x in gcas if 'GCA_' in x]
    return gcas

def getAccession(soup,gca):
    try:
        table = soup.find_all('table')[3]
    except IndexError:
        sys.exit()
    row = table.tbody.find_all('tr')[0]
    acc = row.find_all('td')[-1].text
    return(acc)

def makeRequest(gca):
    baseurl = 'https://www.ncbi.nlm.nih.gov/assembly/'
    gcan = '_'.join(gca.split('_')[:2]) #GCF_01010101.1
    url = baseurl + gcan
    req = requests.get(url)
    return(req)

def correctURL(request):
    if str(request) != '<Response [404]>':
        rev = [x for x in BeautifulSoup(request.text,'html.parser').find_all('a') if 'revision' in x.text]
        if len(rev) != 0:
            return True
        else:
            return False
    else:
        return False

def makeSoup(gca):
    req = makeRequest(gca)
    isCorrect = correctURL(req)
    if isCorrect:
        soup = BeautifulSoup(req.text,'html.parser')
        return soup,gca
    try:
        base = gca.split('_')[1].split('.')[0] #GCA_010101.1_ASM11v1 -> 010101
    except IndexError:
        print(gca)
    bv = int(gca.split('_')[1].split('.')[1])
    for version in range(bv,5):
        for prefix in ['GCA','GCF']:
            ngca = '{}_{}.{}'.format(prefix,base,version)
            req = makeRequest(ngca)
            isCorrect = correctURL(req)
            if isCorrect:
                soup = BeautifulSoup(req.text,'html.parser')
                return soup,ngca
    return False,gca

def getGCANameInfo(gca):
    soup,gca = makeSoup(gca)
    if soup == False:
        return {'GCA':gca,'Accession Number':'none'}
    categories =  [x.text[:-2] for x in soup.find_all('dt')]
    values = [x.text for x in soup.find_all('dd')]
    ddata = {c:v for c,v in zip(categories,values)}
    ddata.update({'GCA':gca})
    ddata.update({'Accession Number':getAccession(soup,gca)})
    return ddata

def main(reffile,outf,processes):
    pool = ThreadPool(processes)
    gcalist = getgcalist(reffile)#[9600:]
    #results = pool.map(getGCANameInfo,urllist)
    results = list(tqdm(pool.imap(getGCANameInfo,gcalist),total=len(gcalist),desc='scrape'))
    results = [x for x in results if x['Accession Number'] != 'none']
    json.dump(results,open(outf,'w'))

if __name__ == '__main__':
    reffile = '/home/sid/thesis_SidReed/data/allCRISPRAnnotationData/CRISPRone_files/mp_CRISPRoneAnnotations_01_11_18SR.json'
    if len(sys.argv) == 3:
        processes = int(sys.argv[2])
    else:
        processes = 8
    main(reffile,sys.argv[1],processes)
