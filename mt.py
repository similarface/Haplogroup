__author__ = 'similarface'
import json
from  collections import defaultdict
mttreejson = './mt-json-new'
with open(mttreejson, 'r') as mttreedata:
    mttree = json.loads(mttreedata.read())
resultdict={}
nolablename={}
def ScanTree(tree):
    name=tree['name']
    children=tree['children']
    rsid=tree['rsid']
    resultdict[name]={'rsid':rsid}
    if name.startswith('NoLabe'):
       nolablename[name]={'rsid':rsid}
    if len(children)!=0:
        for child in children:
            ScanTree(child)
ScanTree(mttree)

def getNolable(tree,lname,res):
    name=tree['name']
    children=tree['children']
    if name==lname:
        for childx in children:
            llname=childx['name']
            if llname.startswith('NoLabel'):
                return  res.append(llname)
    else:
        if len(children)!=0:
            for child in children:
                getNolable(child,lname,res)

print(nolablename)
poss={}
for i in open('/Users/similarface/Documents/yhapmap/mt-dna.txt','r'):
    if i.strip()!='':
        try:
            rsids=i.strip().split('\t')[1].split(' ')
            haplog=i.strip().split('\t')[0]
        except Exception,e:
            # k=i.strip().replace('!','').replace('(','').replace(')','')
            # if k[0] in ['A','T','G','C'] and k[-1] in ['A','T','G','C']:
            #     poss[k[1:-1]]=k[0]+"->"+k[-1]
            #print(i.strip())
            if i.strip().startswith('A') or i.strip().startswith('T') or i.strip().startswith('G') or i.strip().startswith('C') or i.strip().startswith('(') or ('.' in i.strip().split('  ')[0]):
                temprsids=[]
                temprsidssample=[]
                for trs in i.strip().replace(')','').replace('(','').replace('!','').split('  '):
                    if '.' not in trs:
                        temprsids.append(trs)
                        temprsidssample.append(trs[1:-1])
                for k,v in nolablename.items():

                    if len(v['rsid'])==len(temprsids):
                        flag=True
                        for trsid in temprsidssample:
                            if trsid not in v['rsid']:
                                flag=False
                    if flag:
                        for xxx in temprsids:
                            poss[k]={xxx[1:-1]:xxx[0]+'->'+xxx[-1]}
                            print(haplog,k,xxx[1:-1],xxx[0],xxx[-1],i.strip().replace('(','').replace(')','').replace('!',''))
                #print haplog,temprsids
            else:
                poss[i.strip()]={}
        dictrsid={}
        for rsid in rsids:
            if '.' in rsid:
                pass
            if rsid=='':
                pass
            else:
                rs=rsid.replace('!','').replace('(','').replace(')','')
                #poss[rs[1:-1]]=rs[0]+"->"+rs[-1]
                dictrsid[rs[1:-1]]=rs[0]+"->"+rs[-1]
        poss[haplog]=dictrsid

import json
a=json.dumps(poss)
print a

# import json
# mttreejson = './mt-json-new'
# with open(mttreejson, 'r') as mttreedata:
#     mttree = json.loads(mttreedata.read())
# resultdict={}
# def ScanTree(tree):
#     name=tree['name']
#     children=tree['children']
#     rsid=tree['rsid']
#     resultdict[name]={'rsid':rsid}
#     if len(children)!=0:
#         for child in children:
#             ScanTree(child)
# if __name__=='__main__':
#     ScanTree(mttree)
#     print(resultdict)
#
#
#
