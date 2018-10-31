import sys
import re
import requests
from requests_html import HTMLSession

base = 'http://omics.informatics.indiana.edu/CRISPRone/'
testurl='http://omics.informatics.indiana.edu/CRISPRone/check.php?id=GCA_000006605.1_ASM660v1&col=complete'

session = HTMLSession()
r = session.get(testurl)
data = 'No CRISPRS'
for x in r.html.links:
    if 'summary' in x:
        data = requests.get(base + x)
        for line in data.text:
            if 'what' in line:
                print
