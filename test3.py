#coding:utf-8
import json
mtjson = './mt-json'
with open(mtjson, 'r') as ytreedata:
    mttree = json.loads(ytreedata.read())

global ii
ii = 0
def addPid(mttree,a,b):
    global ii
    ii = ii + 1
    a=ii
    name=mttree['name']
    mttree['id']=a
    mttree['pid']=b
    children=mttree['children']
    print name,a,b
    for i in range(0,len(children)):
        addPid(children[i],a,a)

if __name__=='__main__':
    addPid(mttree,ii,0)

# hapdict={}
# for line in open('/Users/similarface/Downloads/test2.txt','r'):
#     lines=line.strip().split(' ')
#     hapdict[lines[0]]={"id":lines[1],"pid":lines[2]}
#
# newfile=open('/Users/similarface/Documents/yhapmap/y-json-new','w')
# for i in open('/Users/similarface/Documents/yhapmap/y-json','r'):
#     if i.strip().replace('"','').startswith("name"):
#         b=i.strip().split(':')[1].strip()
#         name=b[:-1].replace('"','')
#         try:
#             id=hapdict[name]['id']
#             pid=hapdict[name]['pid']
#             #print(i)
#             newfile.writelines(i)
#             newfile.writelines('"id":'+id+",")
#             newfile.writelines('"pid":'+pid+",")
#
#         except KeyError,e:
#             newfile.writelines('------->'+name)
#     else:
#         newfile.writelines(i.strip())
#
# newfile.close()