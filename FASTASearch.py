# -*- coding: utf-8 -*

###
#Author Ryoga Misu
#iGEM group of Tokyo Institute of Technology
###

from bs4 import BeautifulSoup
import re
import urllib.request
import sys
import json
import os
import requests
from selenium import webdriver
import codecs
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
    kind = args.kind
    return (gene, kind)

def makesoup(url):
    with urllib.request.urlopen(url) as response:
        html = response.read()
        time.sleep(5)
        try:
            soup = BeautifulSoup(html, "lxml")
        except:
            soup = BeautifulSoup(html, "html5lib")
    return (soup)

def makesoup_java(url):
    driver = webdriver.PhantomJS(service_log_path=os.path.devnull)
    driver.get(url)
    time.sleep(5)
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

    if 'Showing Current items.' in str(soup.findAll('span', class_='icon')):
        ncbiurl = url
    else:
        hrefname = re.compile('/gene/\w')
        atag = soup.findAll('a', href=hrefname)[2]
        href = str(atag).split(" ")[1]
        ncbiurl = 'https://www.ncbi.nlm.nih.gov' + href.split('"')[1]
    return(ncbiurl)


def ncbi_soup(ncbiurl):
    ncbisoup = makesoup(ncbiurl)
    time.sleep(5)
    return(ncbisoup)

def ncbi_checker(ncbisoup, genename):    
    Egenename = str(ncbisoup.find('dl', id="summaryDl").find('dd', class_='noline')).replace('<', '>').split('>')[2]    
    if Egenename == genename:
        print('The input genename is same with Official Symbol name')

    if not Egenename == genename:
        for line in ncbisoup.find('dl', id="summaryDl").findAll('dd'):
            if '; ' in str(line):
                if genename in str(line).replace('<', '>').replace('>', ' ').replace(';', '').split(' '):
                    print('the Official Symbol of the input genename is ' + Egenename +', founded' )
    
def fasta_from_ncbi(ncbisoup):
    fastaurl = 'https://www.ncbi.nlm.nih.gov' + str(ncbisoup.findAll('a', title="Nucleotide FASTA report")[0]).split('"')[1]
    return(fastaurl)

def seq_fasta(fastaurl, genename, kind):
    fastasoup = makesoup_java(fastaurl)
    comment = '> '+ str(fastasoup.findAll('pre')).split(';')[1].split('\n')[0]
    f = open('output.txt',  'w')
    f.write(comment + '\n')
    for line in fastasoup.findAll('span', class_="ff_line" ):
        f.write(str(line).replace('<', '>').split('>')[2])
    f.close()

if __name__ == '__main__' :

    genename, kind = argparser()

    ncbiurl = url_ncbi(genename, kind)

    ncbisoup = ncbi_soup(ncbiurl)

    ncbi_checker(ncbisoup, genename)

    fastaurl = fasta_from_ncbi(ncbisoup)

    seq_fasta(fastaurl, genename , kind)

    print("Succeeded")
