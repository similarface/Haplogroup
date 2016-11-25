#coding:utf-8

hapdict={}
for line in open('/Users/similarface/Downloads/test3.txt','r'):
    lines=line.strip().split(' ')
    hapdict[lines[0]]={"id":lines[1],"pid":lines[2]}

newfile=open('./mt-json-new','w')
for i in open('./mt-json','r'):
    if i.strip().replace('"','').startswith("name"):
        b=i.strip().split(':')[1].strip()
        name=b[:-1].replace('"','')
        try:
            id=hapdict[name]['id']
            pid=hapdict[name]['pid']
            #print(i)
            newfile.writelines(i)
            newfile.writelines('"id":'+id+",")
            newfile.writelines('"pid":'+pid+",")

        except KeyError,e:
            newfile.writelines('------->'+name)
    else:
        newfile.writelines(i.strip())

newfile.close()

