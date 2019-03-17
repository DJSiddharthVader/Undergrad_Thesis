import os
import re
import sys
import math
import json
import time
import requests
from timeit import timeit
from bs4 import BeautifulSoup
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

#reffile os allCRISPRonelinks.html
base = 'http://omics.informatics.indiana.edu/CRISPRone/' #base for getting reports

def getGCA_Tax_link(lientry):
    reportlink = lientry.a.get('href')
    text = lientry.get_text().split(' | ')
    name = text[0]
    gca = text[2]
    return {'GCA':gca,'Name':name,'Report_URL':reportlink}

def geturllist(reffile):
    allsoup = BeautifulSoup(''.join(open(reffile).readlines()),'html.parser')
    lientries = allsoup.find_all('li')[6:] #skip first 5 entries, headers
    return lientries #list of all report links for each genome

def getCRISPRInfo(reporturl):
    data = requests.get(reporturl)
    dsoup = BeautifulSoup(data.text,'html.parser')
    nothascrispr = list(set(re.findall(r'',dsoup.find_all('h3')[0].text))) #empty if has crispr
    if nothascrispr == ['']: #if has crispr
        #crisprreportlink = [x.get('href') for x in dsoup.find_all('a') if 'what=summary' in x.get('href')][0]
        crisprreportlinks = [x.get('href') for x in dsoup.find_all('a')]
        crl = 'nolinkfound'
        for x in crisprreportlinks:
            if 'what=summary' in x:
                crl = x
        if crl == 'nolinkfound':
            crispr_arrs = -1
            casprots = ['No System']
            numcasprots = -1
            types = 'No System'
        else:
            cdata = requests.get(base + crl)
            crispr_arrs = int(re.search(r'array=(\d+?);',cdata.text).groups(1)[0])
            casprots = re.findall(r'what=cas;des=.*?:(.*?);',cdata.text)
            numcasprots = len(casprots)
            #types = re.search(r'type=(.*?)$',cdata.text).groups(1)
            #the above does not work for some reason, so use grep below
            with open('tmp','w') as gp:
                gp.write(cdata.text)
            types = os.popen("grep -oP 'type=.*$' tmp | tail -1 | cut -d'=' -f2").read().strip('\n')
    else: #no crispr detected
        crispr_arrs = 0
        casprots = ['None']
        numcasprots = 0
        types = ['None']
    return {'No. Arrays (CRISPRone)':crispr_arrs,
            'No. Cas Proteins (CRISPRone)':numcasprots,
            'Cas Proteins (CRISPRone)':casprots,
            'System Types (CRISPRone)':types
           }

def addCRISPRToGCA(gcataxlinkdict):
    crisprannotation = getCRISPRInfo(gcataxlinkdict['Report_URL'])
    gcataxlinkdict.update(crisprannotation)
    return gcataxlinkdict

def main(reffile,outfile):
    linkinfo = [getGCA_Tax_link(x) for x in geturllist(reffile)]
    allinfo = [addCRISPRToGCA(x) for x in linkinfo]
    json.dump(allinfo,open(outfile,'w'))

def main_mp(reffile,outfile,processes):
    pool = ThreadPool(processes)
    #linkinfo = pool.map(getGCA_Tax_link,geturllist(reffile))
    #allinfo = pool.map(addCRISPRToGCA,linkinfo)
    urllist = geturllist(reffile)
    linkinfo = list(tqdm(pool.imap(getGCA_Tax_link,urllist),total=len(urllist),desc='scrape'))
    allinfo = list(tqdm(pool.imap(addCRISPRToGCA,linkinfo),total=len(linkinfo),desc='extracting'))
    json.dump(allinfo,open(outfile,'w'))


if __name__ == '__main__':
    reffile = '/home/sid/thesis_SidReed/data/CRISPRone_files/allCRISPRonelinks.html'
    outfile = str(sys.argv[1])
    if len(sys.argv) == 3:
        processes = int(sys.argv[2])
    else:
        processes = 8
    main_mp(reffile,'mp_' + outfile,processes)


