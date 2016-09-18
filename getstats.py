''' This script takes input of mesh terms from text file and retrieves information regarding number of files
and number of authors and also outputs a visualization the number of publications in three categories
X: Contains single-author publications where each author has three or more publications in his name.
Y: The author appearing in ‘X’ must have three or more collaborative publications.
Z: Single-author Publications authored by those unrelated to ‘X’ having a minimum of two publications(proper) each.
Author:seenu.andi-rajendran'''
i
from Bio import Entrez
from Bio import Medline
import plotly, json, sys, time, httplib, urllib2
import collections
from socket import error as socket_error
import plotly.graph_objs as go
import errno
from Queue import Queue
import threading

LIC = 'pmc cc license'
Entrez.email = "fgregbdbg@bobxx.com"
LIM = 3
batch_size = 10000

#searches for term in database and outputs list of records
def search(Term):
    handle = Entrez.esearch(db='pubmed', term=Term, retmax=100000000)
    record = Entrez.read(handle)
    idlist = record["IdList"]
    count = len(idlist)
    #print(count)
    records=[]
    for start in range(0, count, batch_size):
        #end = min(count, start + batch_size)
        for attempt in range(101):
            try:
                handle = Entrez.efetch(db='pubmed', id=idlist, rettype="medline", retmode="text", retstart=start, retmax=batch_size)
                break
            except Exception as err:
                if attempt==100:
                    sys.exit(str(err))
                else:
                    pass
        records += list(Medline.parse(handle))
    print(Term)
    #print(len(records))
    return records

#counting number of authors who have more than v publications
def filtdic(cntr, v):
    keys = list()
    for key, val in cntr.iteritems():
        if val >= v:
            keys.append(key)
    return keys

#output graph with counts of number of authors for each mesh term
def plt(ls):
    terms = ls.keys()
    id1a = [len(ls[term]['id1a']) for term in terms]
    id1b = [len(ls[term]['id1b']) for term in terms]
    ids2 = [len(ls[term]['id2']) for term in terms]
    trace1 = go.Bar(x=terms, y=id1a, name='Number of authors with ' + str(LIM) + ' or more single-author-publications')
    trace2 = go.Bar(x=terms, y=id1b, name='Number of authors with ' + str(LIM) + ' or more multi-author-publications')
    trace3 = go.Bar(x=terms, y=ids2, name='Number of articles with only one author unrelated to previous')
    data = [trace1, trace2, trace3]
    layout = go.Layout(showlegend=True, barmode='stack')
    fig = go.Figure(data=data, layout=layout)
    plot_url = plotly.offline.plot(fig, filename='graph.html')

#counting all authors in three categories
def getcounts(tm):
    term = '"' + tm + '"' + '[MeSH Terms]' + 'AND ("review"[Publication Type]) AND "journal article"[Publication Type] NOT "letter"[Publication Type] NOT "comment"[Publication Type]'# + ' AND ' + LIC + '[filter]'
    records = search(term)
    recs = [record for record in records if record.has_key("AU")]
    # Get list of authors who have a minimum of three single-author-publications each
    tmpauts = list()
    for rec in recs:
        if len(rec["AU"]) == 1:
            tmpauts.append(rec["AU"][0])
    auts1a = filtdic(collections.Counter(tmpauts), LIM)

    # Get number of authors from the above list having more than 3 colloborative publications
    ids1b = {}
    ids1a = {}

    for aut in auts1a:
        tmpids = []
        cnt = 0
        auttmpids = []
        for rec in recs:
            if len(rec['AU']) == 1 and aut in rec['AU']:
                auttmpids.append({'authors': rec['AU'], 'title': rec['TI'], 'id': rec['PMID'],'date':rec['DP']})
            if len(rec['AU']) > 1 and aut in rec['AU']:
                cnt += 1
                tmpids.append({'authors': rec['AU'], 'title': rec['TI'], 'id': rec['PMID'],'date':rec['DP']})
        if cnt >= LIM:
            ids1b[aut] = tmpids
            ids1a[aut] = auttmpids

    #Pubmed search query format
    term = '"' + tm + '"' + '[MeSH Terms]' + 'AND ("review"[Publication Type]) AND "journal article"[Publication Type] NOT "letter"[Publication Type] NOT "comment"[Publication Type]'#+ ' AND ' + LIC + '[filter]'
    records = search(term)
    recs = [record for record in records if record.has_key("AU")]

    auts2 = [x for x in tmpauts if x not in ids1a.keys() and x not in ids1b.keys()]

    id2 = []

    for aut in list(set(auts2)):
        tmpids = []
        cnt=0
        for rec in recs:
            if len(rec["AU"]) == 1 and rec["AU"][0] == aut:
                tmpids.append({'authors': rec['AU'], 'title': rec['TI'], 'id': rec['PMID'],'date':rec['DP']})
                cnt += 1
        if cnt >= 2:
            id2+=tmpids

    res = {'id1a': ids1a, 'id1b': ids1b,'id2':id2}
    return res

#main function where text file with mesh terms is input and fed to the above functions. Prints json file which is then used to download the files
def main():
    with open('mesh.txt', 'r') as inf:
        terms = [each[:-1] for each in inf.readlines()]
    resdic = dict()
    for term in terms:
        resdic[term] = getcounts(term)
    plt(resdic)
    with open('res.json', 'w') as outf:
        outf.write(json.dumps(resdic, ensure_ascii=False))

#experimental threading function to speed up
def threadmain():
    with open('mesh.txt', 'r') as inf:
        terms = [each[:-1] for each in inf.readlines()]
    q = Queue()
    resq = Queue()
    #print resq.qsize()
    def _worker():
        while True:
            item = q.get()
            resq.put([item, getcounts(item)])
            q.task_done()
    for j in range(1,25):
        t = threading.Thread(target=_worker)
        t.start()
    for term in terms:
        q.put(term)
    q.join()
    resdic={}
    while not resq.qsize()==0:
        oitem=resq.get()
        resdic[oitem[0]]=oitem[1]
    plt(resdic)
    with open('res.json', 'w') as outf:
        outf.write(json.dumps(resdic, ensure_ascii=False))


#mydic = search('"eye diseases"[MeSH Terms] AND pmc cc license[filter]AND ("review"[Publication Type]) AND "journal article"[Publication Type] NOT "letter"[Publication Type] NOT "comment"[Publication Type]')
#print(len(mydic))
main()
