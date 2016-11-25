__author__ = 'similarface'
import json
data={}
for i in open('/Users/similarface/Documents/reslut1/ysnplistposmutaition','r'):
    ii=i.split('\t')
    data[ii[0]]=ii[1].strip()

data_string = json.dumps(data)

file='/Users/similarface/Documents/reslut1/ysnplistposmutaition.json'
a=open(file,'w')
a.writelines(data_string)
a.close()