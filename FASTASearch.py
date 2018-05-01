# coding=utf-8
import sys
import os
import argparse
from Bio import Entrez

Entrez.email = "igem2018tokyotech@gmail.com"

def argparser():
    parser = argparse.ArgumentParser(add_help=True,
                                    prog='GeneINFO.py', # プログラム名
                                    usage='python GeneINFO.py -g/--gene genename -k/--kind kind')
    parser.add_argument("--gene",  "-g", nargs='+', required=True)
    parser.add_argument('--kind', '-k', nargs='+', required=True)
    args = parser.parse_args()
    gene = '+'.join(args.gene)
    kind = '+'.join(args.kind)
    return (gene, kind)

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes" : True, "y" : True, "ye": True, "no" : False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer:" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '': 
            return valid[default]
            break
        elif choice in valid:
            return valid[choice]
            break
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")

def ncbi_for_genbankid(gene, kind):
    print('Fetching Gene Information...\n')
    handle = Entrez.esearch(db="gene", term= gene + '[gene] AND ' + kind  + '[Orgn]')
    record = Entrez.read(handle)
    handle.close()
    geneidlist = record['IdList']
    if len(geneidlist) == 0:
        sys.stdout.write('Bad request. Make sure that all words are spelled correctly or Try different keywords.')
        sys.exit()
        
    else:
        geneidnum = len(geneidlist)
        count = 1
        for geneid in geneidlist :
            handle = Entrez.efetch(db="gene", id=geneid, rettype="gb", retmode="text")
            geneinfo = handle.read()
            answer = query_yes_no("[Information]" + geneinfo + '\nAre you looking for this Gene? Please answer!')
            if answer:
                sys.stdout.write('\nAll right!\n')
                break
            elif answer == False:
                sys.stdout.write('\nOK, searching another choice...')
                count += 1
        if int(geneidnum) < int(count):
            sys.stdout.write("\n\nSorry, can't find what you are looking for.\nMake sure that all words are spelled correctly or Try different keywords.")
            sys.exit()
    
    print('Fetching Genbank ID and Information...\n')
    ncbidict={}
    count = 0
    for info in str(geneinfo).split('\n'):
        count += 1
        if count == 2:
            genename = info.replace('1. ', '')  
        elif 'Official Symbol' in info:
            OfficialSymbol, Name = info.split('and')
            key,value = OfficialSymbol.split(':')
            ncbidict[key] = value
            key, value = Name.split(':')
            ncbidict[key] = value
        else:
            try:
                key,value = info.split(':')
                ncbidict[key] = value
            except ValueError:
                pass

    for anno in ncbidict['Annotation'].split(' '):
        if 'NC_' in anno:
            genbankid = anno
        elif '..' in anno:
            region = anno.replace(',', '').replace('(', '').replace(')', '')
            start, end = region.split('..')
    handle.close()
    return(genename, genbankid, start, end)  

def DefineKind(genename, genbankid, start, end):
    handle = Entrez.efetch(db="nucleotide", id=genbankid, rettype="gb", retmode="text", seq_start = start, seq_stop = end)
    result = handle.read()
    for i in result.split('\n'):
        if '/organism' in i:
            kind = i.replace('/organism=', '').title().replace(' ', '').replace('"', '')
    return(kind)

def Gene_mRNA_CDS(genename, kind, genbankid, start, end):
    os.makedirs('./' + genename + '_in_' + kind, exist_ok = True)
    print('Fetching the DNA sequence of ' + genename +  '...')
    handle = Entrez.efetch(db="nucleotide", id=genbankid, rettype="gb", retmode="text", seq_start = start, seq_stop = end)
    gbresult = handle.read()
    sequence =  ''.join([i for i in gbresult.split('ORIGIN')[1] if not i.isdigit()]).replace('\n', '').replace(' ', '').replace('/', '')
    handle.close()
    f = open('./' + genename + '_in_' + kind  + '/' + genename + '_in_' + kind + '.fa',  'w')
    f.write(sequence)
    f.close()
    print('Fetching the CDS nucleotide of ' + genename +  '...')
    handle = Entrez.efetch(db="nuccore", id=genbankid, rettype="fasta_cds_na", retmode="text", seq_start = start, seq_stop = end)
    CDS = handle.read()
    count = 0
    try:
        for section in CDS.split('>'):
            list = section.split('\n')
            comment = list[0]
            if len(comment) == 0:
                continue
            seq = ''.join(list[1:])
            count += 1
            CDSout = open('./'+ genename + '_in_' + kind + '/' + genename + '_in_' + kind + '_CDS_nucleotide_' + str(count) + '.fa', 'w') 
            CDSout.write('>' + comment + '\n')
            CDSout.write(seq)

    except IndexError:
        print('No CDS_nucleotide founded')

    print('Fetching the CDS protein of ' + genename +  '...')
    handle = Entrez.efetch(db="nuccore", id=genbankid, rettype="fasta_cds_aa", retmode="text", seq_start = start, seq_stop = end)
    CDS_protein = handle.read()
    count = 0
    try:
        for section in CDS_protein.split('>'):
            list = section.split('\n')
            comment = list[0]
            if len(comment) == 0:
                continue
            seq = ''.join(list[1:])
            count += 1
            CDS_protein_out = open('./'+ genename + '_in_' + kind + '/' + genename + '_in_' + kind + '_CDS_protein_' + str(count) + '.fa', 'w')
            CDS_protein_out.write('>' + comment + '\n')
            CDS_protein_out.write(seq + '\n')
    
    except IndexError:
        print('No CDS_protein founded')


def main():
    inputgene, inputkind = argparser()

    genename, genbankid, start, end = ncbi_for_genbankid(inputgene, inputkind)
    
    kind = DefineKind(genename, genbankid, start, end)

    Gene_mRNA_CDS(genename, kind, genbankid, start, end)
    
    print("\nSucceeded")

if __name__ == '__main__':
    main()
