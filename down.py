#program to take the json created by getstats.py and download the files in three folders X,Y,Z

import json
from Bio import Entrez
from Bio import Medline
from bs4 import BeautifulSoup
import sys, os
import shutil

Entrez.email = 'hgfhgf@hgf.com'
reload(sys)
sys.setdefaultencoding('utf-8')

#converting pubmed ids to pubmed central ids as pubmed central has the text files
def idconv(idkeys):
    res=[]
    for i in range(0, len(idkeys), 600):
        for i in range(101):
            try:
                handle = Entrez.elink(db='pmc', dbfrom='pubmed', id=idkeys[i:i + 600])
                break
            except Exception as err:
                if i==100:
                    sys.exit(str(err))
                else:
                    pass
        for each in Entrez.read(handle):
            try:
                res.append(each['LinkSetDb'][0]['Link'][0]['Id'])
            except IndexError:
                pass
    return res

#writing content to file and saving it in respective folder
def filewrite(filename, cont, fold):
    if not os.path.exists(fold): os.mkdir(fold)
    with open(fold + '/' + filename + '.txt', 'wb') as outf:
        outf.write(cont)

#converting html content to text
def getcont(htm):
    bs = BeautifulSoup(htm, "html.parser")
    cont = '\n\n'.join([each.text for each in bs.findAll('p')])
    return cont

#download by pubmed id
def downtxt(pid):
    for i in range(101):
        try:
            handle = Entrez.efetch(db='pmc', id=pid[0], rettype='full', retmode='html')
            break
        except Exception as err:
            if i==100:
                sys.exit(str(err))
            else:
                pass
    htmlstr = handle.read()
    cont = pid[3] + ' by ' + ' ,'.join(pid[2]) + '\n' + getcont(htmlstr)
    return cont

#create folder if not existing
def cpyfil(src, dst):
    if not os.path.exists(dst): os.mkdir(dst)
    shutil.copy(src, dst)

#main function that takes json file as input and saves files by folder whose name is the same as the mesh term
def main():
    with open("res.json", "r")as inf:
        jsondat = inf.read()
    iddic = json.loads(jsondat)
    for term in iddic.keys():
        if not os.path.exists(term): os.makedirs(term)
        kdic = iddic[term]
        tmplst = []
        auts = []
        for k in kdic.keys():
            for kee in kdic[k].keys():
                for each in kdic[k][kee]:
                    tmplst.append([each['id'], each['authors'], each['title'], each['date'].replace(' ','_')])
                    if k == "id1a": auts += each['authors']
        auts = list(set(auts))
        pmids=[each[0] for each in tmplst]
        convids=idconv(pmids)
        ids = zip(convids,pmids,[each[1] for each in tmplst],[each[2] for each in tmplst],[each[3] for each in tmplst])
        for i in ids:
            print(i)
            if len(i[2])>1:
                nm=i[1]+'_collobarative_'+i[4]
            else:
                nm = i[1] +'_'+ i[4]
            filewrite(nm, downtxt(i), term)
            for a in auts:
                if a in i[2]:
                    cpyfil(term + '/' + nm + ".txt", term + '/' + a)
            os.remove(term + '/' + nm + ".txt")


main()
#print (downtxt(['22574783']))
