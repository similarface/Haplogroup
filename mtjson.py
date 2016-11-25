#coding:utf-8
__author__ = 'similarface'
import re
import json
'''
根据http://www.phylotree.org/
获取json字符串
'''

class Node:
    def __init__(self,id, pid, name, rsid):
        self.id = id
        self.pid = pid
        self.name = name
        self.rsid = rsid
        self.children=[]
    def default(self, o):
        return o.__dict__


def getDictFromFile():
    '''
    根据文件获取Dict
    :return:
    '''
    filedict={}
    linenum=0
    for line in open('./mt-dna.txt.new','r'):
        linenum=linenum+1
        tabnum=0
        for i in line:
            tabnum=tabnum+1
            if i!='\t':
                break
        try:
            nolableline=line.replace('\t','')

            if nolableline.startswith(' '):
                nolablelinelist=[xx.replace('!','').replace(')','').replace('(','') for xx in re.split('\s*',nolableline.strip()) if "." not in xx and "-" not in xx and xx!="reserved"]
                filedict[linenum]={'tabmum':tabnum-1,'haplog':"NoLable"+str(linenum),"pos":nolablelinelist}
            else:
                lines=line.strip().split('\t')
                haplog=lines[0].strip()
                try:
                    pos=[yy.replace('!','').replace(')','').replace('(','') for yy in re.split('\s*',lines[1].strip()) if '.' not in yy and "-" not in yy and yy !="reserved"]
                except IndexError,e:
                    pos=[]
                filedict[linenum]={'tabmum':tabnum,'haplog':haplog,"pos":pos}
        except Exception,e:
            print(line)
    return filedict


def getResult():
    '''
    获取层级关系 原理：当前行的tab数目少1并且是当前行的最近的前面就是它的父亲
    :return:
    '''
    result={}
    a=getDictFromFile()
    for i in range(1,len(a.keys())):
        locala=a[i]
        if i==1:
            result[i]=Node(str(i),str(0),locala['haplog'],locala['pos'])
        else:
            pid=scanfileforwardgetpid(a,i,locala['tabmum'])
            result[i]=Node(str(i),str(pid),locala['haplog'],locala['pos'])
    return result


#####################Begin#######################
class NodeWeight:
    def __init__(self,id, pid, name, weight):
        self.id = id
        self.pid = pid
        self.name = name
        self.weight = weight
        self.children=[]
    def default(self, o):
        return o.__dict__

class NodeJSONWeightEncoder(json.JSONEncoder):
    def default(self, node):
        return {"children":node.children,"weight":node.weight,"name":node.name }

def getIdPidFileResult(filename):
    '''
    ID PID name weight
    :return:
    '''
    result={}
    for line in open(filename,'r'):
        lines=line.strip().split(' ')
        result[lines[0]]=NodeWeight(lines[0],lines[1],lines[2],lines[3])
    return result

def WeightTreeNode(result):
    node=None
    for k ,v in result.items():
        currentNode = v;
        parentId = currentNode.pid;
        if parentId!='0':
            parentNode = result[parentId]
            parentNode.children.append(currentNode)
        else:
            node = currentNode
    return NodeJSONWeightEncoder().encode(node)

#####################END#######################


def scanfileforwardgetpid(filedict,forwardnum,suojinshu):
    '''
    :param filedict: 自带呢
    :param forwardnum: 前置位号
    :param suojinshu: 缩进行数
    :return:
    找到当前行的父亲
    '''
    for ii in range(forwardnum,0,-1):
        if filedict[ii]['tabmum']==(suojinshu-1):
            return ii

class NodeJSONEncoder(json.JSONEncoder):
    def default(self, node):
        return {"children":node.children,"rsid":node.rsid,"pid":node.pid,"id":node.id, "name":node.name }

def assembleNode(result):
    node=None
    for k ,v in result.items():
        currentNode = v;
        parentId = currentNode.pid;
        if parentId!='0':
            parentNode = result[int(parentId)]
            parentNode.children.append(currentNode)
        else:
            node = currentNode
    return NodeJSONEncoder().encode(node)

def getTreeFromPid(redict,id):
    ids=[]
    for k,v in redict.items():
        if v['pid']==id:
            ids.append(v['id'])
    return ids

if __name__=='__main__':

    #rootNodes=assembleNode(getResult())
    rootNodes=WeightTreeNode(getIdPidFileResult('/Users/similarface/Documents/reslut1/demo/G-10.0/qilunYweight'))

    print rootNodes
    #print(getDictFromFile())
