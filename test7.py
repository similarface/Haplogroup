#coding:utf-8
__author__ = 'similarface'
from collections import defaultdict
import json
import math
def getmttreejson():
    '''
    根据文件获取mttree字典
    :param filename:
    :return:
    '''
    mttreejson = './mttreedb.json'
    with open(mttreejson, 'r') as mttreedata:
        mttree = json.loads(mttreedata.read())
    return mttree

def scanmttree(tree,re):
    '''
    遍历字典获取snp出现的频率数
    :param tree:
    :param re:
    :return:
    '''
    rsid=tree['rsid']
    for rs in rsid:
        re[rs]=re[rs]+1
    #print(rsid)
    children=tree['children']
    if len(children)!=0:
        for child in children:
            scanmttree(child,re)

def getMaxItem(re):
    '''
    获取最大的次数
    :param re:
    :return:
    '''
    maxcount=0
    for k,v in re.items():
        if v>maxcount:
            maxcount=v
    return maxcount

def getWeightDict():
    '''
    获取权重字典
    算法：F=max(v)
    for i in I put:w(i)=10-9lnf(i)/ln(F)
    :param re:
    :param maxcount:
    :return:
    '''
    re=defaultdict(int)
    mttree=getmttreejson()
    scanmttree(mttree,re)
    weightdict={}
    maxcount=getMaxItem(re)
    for k,v in re.items():
        weightdict[k]=10-round(9*math.log10(v)/math.log10(maxcount),1)
    return weightdict

if __name__ == '__main__':
    print(getWeightDict())