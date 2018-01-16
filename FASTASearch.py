# -*- coding: utf-8 -*

###
#Author Ryoga Misu
#iGEM group of Tokyo Institute of Technology
#
#Usage : python FASTAsearch.py -g/--gene genename -k/--kind kind
# ex) python FASTAsearch.py -g 'CD27' -k 'human'
###

from bs4 import BeautifulSoup
import re
import urllib.request
import sys
import json
import os
import requests
from selenium import webdriver
import argparse
import time

def argparser():
    parser = argparse.ArgumentParser(add_help=True,
                                    prog='GeneINFO.py', # プログラム名
                                    usage='python GeneINFO.py -g/--gene genename -k/--kind kind')
    parser.add_argument("--gene",  "-g", required=True)
    parser.add_argument('--kind', '-k', required=True)
    args = parser.parse_args()
    gene = args.gene
    kind = args.kind.replace(' ', '+')
    return (gene, kind)

def makesoup(url):
    with urllib.request.urlopen(url) as response:
        html = response.read()
        time.sleep(1)
        try:
            soup = BeautifulSoup(html, "lxml")
        except:
            soup = BeautifulSoup(html, "html5lib")
    return (soup)

def makesoup_java(url):
    driver = webdriver.PhantomJS(service_log_path=os.path.devnull)
    driver.get(url)
    time.sleep(3)
    html = driver.page_source.encode('utf-8')  # more sophisticated methods may be available
    try:
        soup_java = BeautifulSoup(html, "lxml")
    except:
        soup_java = BeautifulSoup(html, "html5lib")
    return(soup_java)
    
def url_ncbi(genename, kind):
    url = 'https://www.ncbi.nlm.nih.gov/gene/?term=' + genename + '+' + kind
    soup = makesoup(url)
    if  'The following term was not found in Gene' in str(soup.findAll('span', class_='icon')):
        print('[ERROR] The following term was not found in Gene')
        sys.exit()

    elif 'No items found.' in str(soup.findAll('span', class_='icon')) :
        print('[ERROR] No items found.')
        sys.exit()
    else :
        pass

    if 'Showing Current items.' in str(soup.findAll('span', class_='icon')) and not '<h2>Search results</h2>' in str(soup.findAll('h2')):
        ncbiurl = url
    else:
        hrefname = re.compile('/gene/\w')
        atag = soup.findAll('a', href=hrefname)[2]
        href = str(atag).split(" ")[1]
        ncbiurl = 'https://www.ncbi.nlm.nih.gov' + href.split('"')[1]
    return(ncbiurl)


def ncbi_soup(ncbiurl):
    ncbisoup = makesoup(ncbiurl)
    return(ncbisoup)

def ncbi_checker(ncbisoup, genename):    
    Egenename = str(ncbisoup.find('dl', id="summaryDl").find('dd', class_='noline')).replace('<', '>').split('>')[2]    
    if Egenename.lower() == genename.lower():
        print('The input genename is same with Official Symbol')

    elif not Egenename == genename:
        for line in ncbisoup.find('dl', id="summaryDl").findAll('dd'):
            if '; ' in str(line):
                if genename.lower() or genename.upper() in str(line).replace('<', '>').replace('>', ' ').replace(';', '').split(' '):
                    print('The Official Symbol of the input Genename is ' + Egenename +', founded' )
    
def fasta_from_ncbi(ncbisoup):
    fastaurl = 'https://www.ncbi.nlm.nih.gov' + str(ncbisoup.findAll('a', title="Nucleotide FASTA report")[0]).split('"')[1]
    return(fastaurl)

def seq_fasta(fastaurl, genename, kind):
    try:
        os.mkdir('./' + genename + '_in_' + kind)
    except:
        pass
    fastasoup = makesoup_java(fastaurl)
    comment = '> '+ str(fastasoup.findAll('pre')).split(';')[1].split('\n')[0]
    f = open('./' + genename + '_in_' + kind  + '/' + genename + '_in_' + kind + '.fa',  'w')
    f.write(comment + '\n')
    for line in fastasoup.findAll('span', class_="ff_line" ):
        f.write(str(line).replace('<', '>').split('>')[2])
    f.close()

def Genbank_soup(ncbisoup):
    genbankurl = 'https://www.ncbi.nlm.nih.gov' + str(ncbisoup.findAll('a', title="Nucleotide GenBank report" )).split('"')[1]
    genbanksoup = makesoup_java(genbankurl)
    return(genbanksoup)

def Genbank_mRNA_CDS(Genbanksoup, genename, kind):
    mRNA = re.compile('\wmRNA_0')
    mRNAregion = str(str(Genbanksoup.findAll('span', id=mRNA, class_="feature")[0]).split('join')[1].split('/gene')[0]).replace('\n', '').replace(' ', '').replace('(', '').replace(')', '').split(',')
    CDS = re.compile('\wCDS_0')
    CDSregion = str(str(Genbanksoup.findAll('span', id=CDS, class_="feature")[0]).split('join')[1].split('/gene')[0]).replace('\n', '').replace(' ', '').replace('(', '').replace(')', '').split(',')
    fa = open('./'+ genename + '_in_' + kind  + '/' + genename + '_in_' + kind + '.fa', 'r')
    seq = fa.read().split('\n')[1]
    fa.close()
    mRNAout = open('./'+ genename + '_in_' + kind + '/' + genename + '_in_' + kind + '_mRNA.fa', 'w')
    for num in mRNAregion:
        first, second = num.split('..')
        start = int(first) - 1
        end = int(second)
        mRNAout.write(str(seq)[start:end].strip())
    mRNAout.close()

    CDSout = open('./'+ genename + '_in_' + kind + '/' + genename + '_in_' + kind + '_CDS.fa', 'w')
    for num in CDSregion:
        first, second = num.split('..')
        start = int(first) - 1
        end = int(second)
        CDSout.write(str(seq)[start:end].strip())
    CDSout.close()


if __name__ == '__main__' :

    genename, kind = argparser()

    ncbiurl = url_ncbi(genename, kind)

    ncbisoup = ncbi_soup(ncbiurl)

    ncbi_checker(ncbisoup, genename)

    fastaurl = fasta_from_ncbi(ncbisoup)

    seq_fasta(fastaurl, genename , kind)

    Genbanksoup = Genbank_soup(ncbisoup)

    Genbank_mRNA_CDS(Genbanksoup, genename, kind)

    print("Succeeded")
