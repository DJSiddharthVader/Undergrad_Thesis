import os
import re
import sys
import json
import requests
from tqdm import tqdm
from bs4 import BeautifulSoup
from multiprocessing.dummy import Pool as ThreadPool

baseURL = 'http://omics.informatics.indiana.edu/CRISPRone/' #base for getting reports
#takes about ~15 minutes working on 11102 taxa using 4 cores (processes)

def getGCATaxLink(entry):
    reportlink = entry.a.get('href')
    info = entry.get_text().split(' | ')
    return {'GCA':info[2],'Name':info[0],'Report_URL':reportlink}

def getCRISPRInfo(reporturl):
    #for a link to a crispr_one page (1 per strain) scrape info about the crispr genes in the strain
    dsoup = BeautifulSoup(requests.get(reporturl).text,'html.parser')
    has_crispr = True if list(set(re.findall(r'',dsoup.find_all('h3')[0].text)))  == [''] else False #check if report says taxa has a crispr system  or not
    if has_crispr:
        reportlinks = [x.get('href') for x in dsoup.find_all('a')]
        has_info = [x if 'what=summary' in x else False for x in reportlinks]
        reportlink = list(filter(lambda a:a, has_info))[-1] if any(x for x in has_info) else 'nolinkfound' #get link to info about crisprs found in taxa
        if reportlink == 'nolinkfound':
            # no crispr information found
            crispr_arrs = -1
            casprots = ['No System']
            numcasprots = -1
            types = 'No System'
        else:
            #get crispr information found
            cdata = requests.get(baseURL + reportlink)
            crispr_arrs = int(re.search(r'array=(\d+?);',cdata.text).groups(1)[0])
            casprots = re.findall(r'what=cas;des=.*?:(.*?);',cdata.text)
            numcasprots = len(casprots)
            with open('tmp','w') as gp:
                gp.write(cdata.text)
            types = os.popen("grep -oP 'type=.*$' tmp | tail -1 | cut -d'=' -f2").read().strip('\n') # grep for the types of crispr systems detected
            #types = re.search(r'type=(.*?)$',cdata.text).groups(1)
            #the above does not work for some reason, so use grep
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

def addCRISPRToGCA(gcaDict):
    crisprannotation = getCRISPRInfo(gcaDict['Report_URL'])
    gcaDict.update(crisprannotation)
    return gcaDict

def getURLList(reffile):
    soup = BeautifulSoup(open(reffile),'html.parser')
    entries = soup.find_all('li')[6:] #skip first 5 entries, headers
    return entries #list of all report links for each genome

def main(reffile,outfile,processes=1):
    if processes > 1:
        pool = ThreadPool(processes)
        urllist = getURLList(reffile)
        linkinfo = list(tqdm(pool.imap(getGCATaxLink,urllist),total=len(urllist),desc='scrape'))
        allinfo = list(tqdm(pool.imap(addCRISPRToGCA,linkinfo),total=len(linkinfo),desc='extracting'))
    else:
        linkinfo = [getGCATaxLink(x) for x in getURLList(reffile)]
        allinfo = [addCRISPRToGCA(x) for x in linkinfo]
    json.dump(allinfo,open(outfile,'w'))


if __name__ == '__main__':
    reffile = sys.argv[1]#/home/sid/thesis_SidReed/data/CRISPROne.html
    outfile = sys.argv[2]#/home/sid/thesis_SidReed/data/CRISPROne.html.json
    processes = int(sys.argv[3])
    main(reffile,outfile,processes)


