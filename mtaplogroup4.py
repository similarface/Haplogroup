# coding:utf-8
__author__ = 'similarface'
import snp
import json
import os
from collections import deque
import logging
import logging.handlers

LOG_FILE = 'mtaplogroup.log'
handler = logging.handlers.RotatingFileHandler(LOG_FILE, maxBytes=1024 * 10240, backupCount=5)  # 实例化handler
fmt = '%(asctime)s - %(filename)s:%(lineno)s - %(name)s - %(message)s'
formatter = logging.Formatter(fmt)  # 实例化formatter
handler.setFormatter(formatter)  # 为handler添加formatter

'''
读取23andme的芯片数据
'''

mttreejson = './mttreedb.json'

with open(mttreejson, 'r') as mttreedata:
    mttree = json.loads(mttreedata.read())
from mtjson import getDictFromFile
from test7 import getWeightDict
def getTotalWeightTree(sdata):
    weightTree=0
    weightdict=getWeightDict()
    redirecttree=getDictFromFile()
    for k,v in redirecttree.items():
        rsids=v['pos']
        if len(rsids)!=0:
            for rsid in rsids:
                try:
                    if rsid[-1] == sdata[rsid[1:-1]].genotype[0]:
                        weightTree=weightTree+weightdict[rsid]
                except Exception,e:
                    pass
    return weightTree

weightTree=getTotalWeightTree(snp.reader_dict_mtpos('/Users/similarface/wocao.txt'))

def getMthaologroup(mttree, sdata, result):
    '''
    :param ytree: Y染色体的haplogroup树
    :param ysnp: Y染色体突变字典
    :param sdata: 芯片数据
    :param result: 结果
    :return:
    '''

    id = mttree['id']
    pid = mttree['pid']
    name = mttree['name']
    children = mttree['children']
    rsids = mttree['rsid']
    # 匹配的rsid
    matchrsid = []
    # 不匹配的rsid
    nomathrsid = []
    weight = 0
    weightdict=getWeightDict()
    qfeziweight=0
    qfenmuweight=0
    lain=0

    zerolain=2.01
    for rsid in rsids:
        try:
            if rsid[-1] == sdata[rsid[1:-1]].genotype[0]:
                logger.debug("pos {0},mutation {1},name {2} ,match {3}".format(rsid,sdata[rsid[1:-1]].genotype[0], name,'YES'))
                matchrsid.append(rsid)
                weight = weight + 1
                qfeziweight=weightdict[rsid]+qfeziweight
                qfenmuweight=weightdict[rsid]+qfenmuweight
                lain=lain+1
                zerolain=0
            else:
                logger.debug("pos {0},mutation {1},name {2} ,match {3}".format(rsid,sdata[rsid[1:-1]].genotype[0], name,'NO'))
                #weight=weight-1
                qfenmuweight=weightdict[rsid]+qfenmuweight
                nomathrsid.append(rsid)
        except KeyError, e:
            qfenmuweight=weightdict[rsid]+qfenmuweight
            nomathrsid.append(rsid)
            try:
                logger.debug("pos {0},mutation {1},name {2} ,match {3}".format(rsid, '?', name,
                                                                                       'UNKOWNU'))
            except KeyError, e:
                try:
                    logger.debug(
                        "pos {0},mutation {1},name {2} ,match {3}".format(rsid, sdata[rsid[1:-1]].genotype[0],
                                                                                  name, 'UNKOWNU'))
                except KeyError, e:
                    logger.debug(
                        "pos {0},mutation {1},name {2} ,match {3}".format(rsid, '?', name, 'UNKOWNU'))

    #currentfenziweight=result[pid]['qfeziweight']+qfeziweight
    currenthoumianfenmuweight=result[pid]['qfenmuweight']+qfenmuweight
    #print(id,currentfenziweight,currentfenmuweight)
    if qfenmuweight!=0:
        result[id] = {'pid': pid, 'name': name, 'matchrsid': matchrsid, 'nomathrsid': nomathrsid,
                  'weight': (0.5*qfeziweight/qfenmuweight+0.5*qfeziweight/currenthoumianfenmuweight),
                      'lain':lain,"qfeziweight":qfeziweight,
                  "qfenmuweight":currenthoumianfenmuweight}
    else:
        result[id] = {'pid': pid, 'name': name, 'matchrsid': matchrsid, 'nomathrsid': nomathrsid,
                  'weight': 0+0.5*qfeziweight/currenthoumianfenmuweight,
                      'lain':lain,"qfeziweight":qfeziweight,
                  "qfenmuweight":currenthoumianfenmuweight}


    print(id,result[id]['weight'])
    if len(children) != 0:
        for itemdict in children:
            getMthaologroup(itemdict, sdata, result)
    else:
        pass


def getSnpCountSum(ytree):
    name=ytree['rsid']
    global allrsid
    allrsid=allrsid+name
    children=ytree['children']
    if len(children)!=0:
        for i in range(0,len(children)):
            getSnpCountSum(children[i])



def getPath(mttree, deques, id):
    '''
    根据ID 获取树的路径
    :param y-tree等级树
    :param deques: 是个双端队列 可以前面插入值
    :param id: 树中的某个节点
    :return:返回队列 deques
    '''
    # id=0表示到达树的顶层
    if id == "0":
        return ''
    # 左侧加入该节点的父亲节点
    try:
        deques.appendleft({mttree[id]['name']:{'matchrsid':mttree[id]['matchrsid'],'lain':mttree[id]['lain'],'nomatchrsid':mttree[id]['nomathrsid']}})
    except Exception,e:
        print(id)
    #deques.appendleft({mttree[id]['name']: mttree[id]['matchrsid']})

    # 获取ytree树的父节点的PID
    pid = mttree[id]['pid']
    getPath(mttree, deques, pid)


def getMaxWeightNode(re):
    # 权重值
    imax = -1000
    # 字典的key
    maxid = 0
    lre={}
    for k,v in re.items():
        lre[k]=v['weight']
    return sorted(lre.iteritems(),key=lambda asd:asd[1],reverse=True)


def getMthaologroupResult(slugfilename):
    # 芯片数据
    sdata = snp.reader_dict_mtpos(slugfilename)
    # 存储运算结果
    result = {"0": {'matchrsid': None, 'nomathrsid': None, 'weight': 0, 'name': None,'lain':0,'qfeziweight':0,'qfenmuweight':0}}
    #
    getMthaologroup(mttree, sdata, result)
    return result


def getMthaologroupPath(slugfilename):
    result = getMthaologroupResult(slugfilename)
    # 获取权重最大的节点
    alltreepath=[]
    i=0
    oldmaxid=-100000
    lenzeropath=-100000
    retreepath=None
    for node in getMaxWeightNode(result):
        i=i+1
        if i>10 :
            break
        maxid, imax = node
        if oldmaxid>imax:
            break
        oldmaxid=imax
        logger.info('maxid {0} ,imax {1}'.format(maxid, imax))
        treepath = deque([])
        # 获取权重最大的节点的路径
        #getPath(result, treepath, maxid)
        getPath(result, treepath, maxid)

        for item in treepath:
            for k,v in item.items():
                print(k,v['matchrsid'])

        logger.info('->'.join(getFormatYhaologroupPath(treepath)[0]))
        #logger.info('>'.join(getFormatYhaologroupPath(treepath)[1]))
        zeropath=getFormatYhaologroupPath(treepath)[1]
        #print '->'.join(getFormatYhaologroupPath(treepath)[0])
        print '->'.join(getFormatYhaologroupPath(treepath)[0])
        if lenzeropath<len(zeropath):
            lenzeropath=len(zeropath)
            retreepath=treepath

        print(getZeroCount(zeropath),getSumLain(zeropath),len(zeropath))
        #print(treepath)
        alltreepath.append(treepath)
    print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    print '->'.join(getFormatYhaologroupPath(treepath)[0])
    return retreepath


def getFormatYhaologroupPath(treepath):
    re = []
    lain=[]
    for item in treepath:
        for k,v in item.items():
            re.append(k)
            lain.append(v['lain'])
    return re,lain

def getZeroCount(lain):
    count=0
    for i in lain:
        if i==0:
            count=count+1
    return count

def getSumLain(lain):
    sums=0
    for i in lain:
        sums=sums+i
    return sums

def getYhaologroupJson():
    pass

from update23andmefilepos import update23andmeFilePos
import sys
if __name__ == "__main__":
    #---sample 01----无争议
    #N 60 H1e A2a5 [16189,7028] [deque([{u'MT-DNA': []}, {u"L1'2'3'4'5'6": []}, {u"L2'3'4'5'6": []}, {u"L2'3'4'6": [u'16189']}, {u"L3'4'6": []}, {u"L3'4": []}, {u'L3': []}, {u'N': []}, {u'R': []}, {u'R0': []}, {u'HV': []}, {u'H': [u'7028']}, {u'H1': []}, {u'H1e': []}]
    #filename = '/Users/similarface/Documents/user60_file27_yearofbirth_1982_sex_XY.23andme.txt'
    #filename='/Users/similarface/Documents/reslut1/genome_wang_mian_Full_20151021231213.txt'
    filename='/Users/similarface/wocao.txt'
    #filename='/Users/similarface/Documents/user10_file3_yearofbirth_1982_sex_XY.23andme.txt'

    #filename=sys.argv[1]
    #sample 02 -->争议
    #N 1785  H3     T2e1a1a
    #filename = '/Users/similarface/Documents/user1785_file1041_yearofbirth_unknown_sex_unknown.23andme.txt'
    #filename = '/Users/similarface/Downloads/11.23andme.176'

    #sample 03 -->
    #N 1875  HV9a   T2e1a1a
    #filename = '/Users/similarface/Documents/user1875_file1113_yearofbirth_unknown_sex_unknown.23andme.txt'


    loggername = os.path.basename(filename)[0:os.path.basename(filename).rindex('.')]
    logger = logging.getLogger(loggername)  # 获取名为tst的logger
    logger.addHandler(handler)  # 为logger添加handler
    logger.setLevel(logging.DEBUG)

    # for file in ['/Users/similarface/Documents/reslut1/genome_wang_mian_Full_20151021231213.txt',
    #              '/Users/similarface/Documents/reslut1/genome_wang_jun_chen_Full_20151021231123.txt',
    #              '/Users/similarface/Downloads/I2b1.txt',
    #              '/Users/similarface/Downloads/H6a1b2.txt',
    #              '/Users/similarface/Downloads/1718.23andme.993',
    #              '/Users/similarface/Downloads/1718.23andme.993',
    #              '/Users/similarface/Downloads/4763.23andme.3353',
    #              '/Users/similarface/Downloads/T',
    #              '/Users/similarface/Documents/wegeneRawData.txt']:
    #update23andmeFilePos(filename)
    getMthaologroupPath(filename)
    #sdata = snp.reader_dict_mtpos(filename)
    #print(getTotalWeightTree(sdata))
    # getSnpCountSum(ytree)
    # rsids=list(set(allrsid))
    # print(','.join(rsids))
    # print len(rsids)
