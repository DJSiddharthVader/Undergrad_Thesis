import sys
import json
import requests
import itertools
import pandas as pd
from tqdm import tqdm
from bs4 import BeautifulSoup
from functools import partial
from multiprocessing.dummy import Pool as ThreadPool

baseurl = 'https://www.ncbi.nlm.nih.gov/assembly/'
prefixes = ['GCA','GCF']
"""
Will return an error if 8 >= processes, not with 4
    requests.exceptions.ConnectionError: ('Connection aborted.', RemoteDisconnected('Remote end closed connection without response'))
think it has to do with  making too many connections at a time so NCBI just denies the request
Maybe adding a sleep after every requests will resolve error, but may end up being slower anyways
Takes ~1.5-2 hours looking at 11091 GCAs using 8 processes

Failed Accessions:
  - GCA_000090405.1_ASM9040v2     (3648)
  - GCA_000009225.1_ASM922v1      (3924)
  - GCA_000967155.1_HUSEC2011CHR1 (9199)
"""

def getGCAList(reffile):
    #get GCA accessions for organisms, use these to map bac. strains to fasta files
    gcas = list(pd.DataFrame(json.load(open(reffile)))['GCA'])
    gcas = [x for x in gcas if 'GCA_' in x]
    gcas.sort() #for error tracing if there are specific problem gcas
    return gcas

def getAccession(soup):
    #soup is the python object of a web page for parsing
    #given the NCBI web page soup object, extract the bac. straina accession
    try:
        table = soup.find_all('table')[3]
    except IndexError:
        sys.exit()
    row = table.tbody.find_all('tr')[0]
    acc = row.find_all('td')[-1].text
    return(acc)

def makeRequest(gca):
    #fetch NCBI webpage for a given GCA accession
    gcan = '_'.join(gca.split('_')[:2]) #GCF_01010101.1, accession
    url = baseurl + gcan
    req = requests.get(url)
    return(req)

def checkURL(request):
    #checks if url is not dead, else returns false
    if str(request) != '<Response [404]>':
        rev = [x for x in BeautifulSoup(request.text,'html.parser').find_all('a') if 'revision' in x.text]
        if len(rev) != 0:
            return True
        else:
            return False
    else:
        return False

def makeGCAList(gca,max_version):
    try:
        base = gca.split('_')[1].split('.')[0] #GCA_010101.1_ASM11v1 -> 010101
    except IndexError:
        print(gca)
        return False
    base_version = int(gca.split('_')[1].split('.')[1]) #get version from accession number
    url_list = ['{}_{}.{}'.format(prefix,base,version) for prefix,version in itertools.product(prefixes,list(range(base_version,max_version)))]
    return url_list

def makeSoup(gca,max_version=5):
    #Gets NCBI page for a gca number, checking multiple version numbers, returing a soup object
    ugca = gca
    req = makeRequest(gca)
    if checkURL(req): #if the url returns a webpage, return the soup object
        soup = BeautifulSoup(req.text,'html.parser')
    else: #else  try the same url with the version number +1 until it works or reaches maxversion
        url_list = makeGCAList(gca,max_version)
        if url_list == False: #no base version number
            soup = False
        for ugca in url_list:
            req = makeRequest(ugca)
            if checkURL(req):
                soup = BeautifulSoup(req.text,'html.parser')
        soup = False
    return soup,ugca

def getGCAInfo(gca):
    #For a gca number extract information about the organism from the NCBI page and save as a dictionary
    try:
        soup,gca = makeSoup(gca)
    except requests.exceptions.ConnectionError:
        print('Could not fetch {}'.format(gca))
        return {'GCA':gca,'Accession Number':'Failed Connection'}
    if soup == False:
        return {'GCA':gca,'Accession Number':'none'}
    categories =  [x.text[:-2] for x in soup.find_all('dt')]
    values = [x.text for x in soup.find_all('dd')]
    ddata = {c:v for c,v in zip(categories,values)}
    ddata.update({'GCA':gca})
    ddata.update({'Accession Number':getAccession(soup)})
    return ddata

def main(reffile,outf,processes):
    pool = ThreadPool(processes)
    gcalist = getGCAList(reffile)#[9600:]
    #session = requests.Session()
#    scrape_fnc = partial(getGCAInfo,session=session)
#    results = list(tqdm(pool.imap(scrape_fnc,gcalist),total=len(gcalist),desc='scrape'))
    results = list(tqdm(pool.imap(getGCAInfo,gcalist),
                        total=len(gcalist),
                        desc='scrape'))
    results = [x for x in results if x['Accession Number'] != 'none']
    json.dump(results,open(outf,'w'))


if __name__ == '__main__':
    reffile = sys.argv[1] #'/home/sid/thesis_SidReed/data/CRISPROne.json'
    outfile = sys.argv[2] #/home/sid/thesis_SidReed/data/NCBI_data.json
    processes = int(sys.argv[3])
    main(reffile,outfile,processes)

#DEPRECIATED
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


