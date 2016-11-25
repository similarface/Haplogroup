#coding:utf-8
__author__ = 'similarface'

file='/Users/similarface/Downloads/Root.txt'


infodict={}
ii=0
for line in open(file):
    headtag=0
    ii=ii+1
    if line.startswith(".") :
        for char in line:
            if char=='.':
                headtag=headtag+1
            elif char==' ':
                pass
            else:
                break
    infodict[ii]={'headtag':headtag,'data':line}

def getParent(infodict,col):
    for i in range(col,0,-1):
        if infodict[i]['headtag']==infodict[col]['headtag']-1:
            return i
import re
def gethaplop(line):
    if "~" in line:
        return line[:line.index("~")+1].strip()
    else:
        return re.split('\s*',line.strip())[0].strip()

def getdata(line):
    if "~" in line:
        return line[line.index("~"):].strip().replace('/',', ').replace('~','').strip()
    else:
        return  line[len(re.split('\s*',line)[0])+1:].strip().replace('/',', ').replace('~','').replace('\n','')

dictlen=infodict.viewkeys().__len__()

for jj in range(1,dictlen):
    print("insert into T_Y_TREE VALUES ("+str(jj)+","+str(getParent(infodict,jj))+","+str(infodict[jj]['headtag'])+",'"+gethaplop(infodict[jj]['data'][infodict[jj]['headtag']*2:])+"','"+getdata(infodict[jj]['data'][infodict[jj]['headtag']*2:])+"');")

'''
create table T_Y_TREE(
id int,
pid int,
levels int,
snp varchar(64),
markers VARCHAR (4000)
)
'''