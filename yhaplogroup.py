# coding:utf-8
__author__ = 'similarface'
import snp
import json
import os
from collections import deque
import logging
import logging.handlers

LOG_FILE = 'yhaplogroup.log'
handler = logging.handlers.RotatingFileHandler(LOG_FILE, maxBytes=1024 * 1024, backupCount=5)  # 实例化handler
fmt = '%(asctime)s - %(filename)s:%(lineno)s - %(name)s - %(message)s'
formatter = logging.Formatter(fmt)  # 实例化formatter
handler.setFormatter(formatter)  # 为handler添加formatter
'''
读取23andme的芯片数据
'''

ytreejson = './y-json-new'
with open(ytreejson, 'r') as ytreedata:
    ytree = json.loads(ytreedata.read())
# Y染色体突变库
ysnpjson = './ysnp.json'
with open(ysnpjson, 'r') as ysnpjsondata:
    ysnp = json.loads(ysnpjsondata.read())

def getYhaologroup(ytree, ysnp, sdata, result):
    '''
    :param ytree: Y染色体的haplogroup树
    :param ysnp: Y染色体突变字典
    :param sdata: 芯片数据
    :param result: 结果
    :return:
    '''
    id = ytree['id']
    pid = ytree['pid']
    name = ytree['name']
    children = ytree['children']
    rsids = ytree['rsid']
    # 匹配的rsid
    matchrsid = []
    # 不匹配的rsid
    nomathrsid = []
    weight = 0
    selfweight=0
    for rsid in rsids:
        try:
            if ysnp[rsid][-1] == sdata[rsid].genotype[0]:
                logger.debug("pos {0},mutation {1},snp {2},name {3} ,match {4}".format(rsid, ysnp[rsid][-1],
                                                                                       sdata[rsid].genotype[0], name,
                                                                                       'YES'))
                matchrsid.append(rsid)
                weight = weight + 1
                selfweight=selfweight+1
            else:
                logger.debug("pos {0},mutation {1},snp {2},name {3} ,match {4}".format(rsid, ysnp[rsid][-1],
                                                                                         sdata[rsid].genotype[0], name,
                                                                                         'NO'))
                nomathrsid.append(rsid)
        except KeyError, e:
            nomathrsid.append(rsid)
            try:
                logger.debug("pos {0},mutation {1},snp {2},name {3} ,match {4}".format(rsid, ysnp[rsid][-1], '?', name,
                                                                                       'UNKOWNU'))
            except KeyError, e:
                try:
                    logger.debug(
                        "pos {0},mutation {1},snp {2},name {3} ,match {4}".format(rsid, '?', sdata[rsid].genotype[0],
                                                                                  name, 'UNKOWNU'))
                except KeyError, e:
                    logger.debug(
                        "pos {0},mutation {1},snp {2},name {3} ,match {4}".format(rsid, '?', '?', name, 'UNKOWNU'))

    result[id] = {'pid': pid, 'name': name, 'matchrsid': matchrsid, 'nomathrsid': nomathrsid,
                  'weight': result[pid]["weight"] + weight,'selfweight':selfweight}
    if len(children) != 0:
        for itemdict in children:
            getYhaologroup(itemdict, ysnp, sdata, result)
    else:
        pass

global allrsid
allrsid=[]

def getSnpCountSum(ytree):
    name=ytree['rsid']
    global allrsid
    allrsid=allrsid+name
    children=ytree['children']
    if len(children)!=0:
        for i in range(0,len(children)):
            getSnpCountSum(children[i])



def getPath(ytree, deques, id):
    '''
    根据ID 获取树的路径
    :param y-tree等级树
    :param deques: 是个双端队列 可以前面插入值
    :param id: 树中的某个节点
    :return:返回队列 deques
    '''
    # id=0表示到达树的顶层
    if id == 0:
        return ''
    # 左侧加入该节点的父亲节点
    deques.appendleft({ytree[id]['name']:{'matchrsid':ytree[id]['matchrsid'],'selfweight':ytree[id]['selfweight'],'id':id,'pid':ytree[id]['pid'],'nomatchrsid':ytree[id]['nomathrsid']}})
    #deques.appendleft({ytree[id]['name']: ytree[id]['matchrsid']})

    # 获取ytree树的父节点的PID
    pid = ytree[id]['pid']
    getPath(ytree, deques, pid)


def getMaxWeightNode(re):
    # 权重值
    imax = -10000
    # 字典的key
    maxid = 0
    for k, v in re.items():
        if v['weight'] > imax:
            imax = v['weight']
            maxid = k
    return maxid, imax


def getYhaologroupResult(slugfilename):
    # 芯片数据
    sdata = snp.reader_dict_ypos(slugfilename)
    # 存储运算结果
    result = {0: {'matchrsid': None, 'nomathrsid': None, 'weight': 0, 'name': None,'selfweight':0,'pid':0}}
    #
    getYhaologroup(ytree, ysnp, sdata, result)
    return result


def getYhaologroupPath(slugfilename):
    result = getYhaologroupResult(slugfilename)
    for kkk,vvv in result.items():
        print kkk,vvv['pid'],vvv['name'],vvv['selfweight']
    # 获取权重最大的节点
    maxid, imax = getMaxWeightNode(result)
    logger.info('maxid {} ,imax {}'.format(maxid, imax))
    treepath = deque([])
    # 获取权重最大的节点的路径
    getPath(result, treepath, maxid)
    # for que in treepath:
    #         for kkk,vvv in que.items():
    #             print vvv['id'],vvv['pid'],kkk,vvv['selfweight']
    logger.info(treepath)
    logger.info('->'.join(getFormatYhaologroupPath(treepath)))
    print '->'.join(getFormatYhaologroupPath(treepath))
    return treepath


def getFormatYhaologroupPath(treepath):
    re = []
    for item in treepath:
        re.append(item.keys()[0])
    return re


def getYhaologroupJson():
    pass


if __name__ == "__main__":
    #filename = "/Users/similarface/Documents/reslut1/genome_wang_mian_Full_20151021231213.txt"
    filename='/Users/similarface/Documents/FangGE.txt'

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
    getYhaologroupPath(filename)

    # getSnpCountSum(ytree)
    # rsids=list(set(allrsid))
    # print(','.join(rsids))
    # print len(rsids)
