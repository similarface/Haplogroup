__author__ = 'similarface'
import re

jrsidlist=[]
jchrpos=[]
for j in open('/Users/similarface/Documents/chrposrsidxinpian.data','r'):
    jlines=j.strip().split(' ')
    jrsid=jlines[0]
    jchr=jlines[1]
    jpos=jlines[2]
    jrsidlist.append(jrsid)
    jchrpos.append(jchr+':'+jpos)

match=0
matchlist=[]

nomatch=0
nomatchlist=[]

for i in open('/Users/similarface/Documents/20160331SecondPrimerSNP.bed','r'):
    lines=re.split('\s*',i.strip())
    try:
        rsid=lines[5]
        chr=lines[0]
        pos=lines[2]
        if rsid in jrsidlist:
            match=match+1
            matchlist.append(chr+':'+pos)
        elif chr+':'+pos in jchrpos:
            match=match+1
            matchlist.append(chr+':'+pos)
        else:
            nomatch=nomatch+1
            nomatchlist.append(chr+':'+pos)
    except IndexError,e:
        chr=lines[0]
        pos=lines[2]
        if chr+':'+pos in jchrpos:
            match=match+1
            matchlist.append(chr+':'+pos)
        else:
            nomatch=nomatch+1
            nomatchlist.append(chr+':'+pos)
print("match:"+str(match))
print("nomatch:"+str(nomatch))
print(matchlist)
print(nomatchlist)