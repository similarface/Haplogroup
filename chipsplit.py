#coding:utf-8
__author__ = 'similarface'
import snp
import MySQLdb
from collections import defaultdict
'''
读取23andme的芯片数据
'''
import csv
class SnpCompare(snp.snp):
    def __init__(self,name,chromosome,position,genotype,haplogroup,treesnp,mutation,flag,marker):
        snp.snp.__init__(self,name,chromosome,position,genotype)
        self.mutation=mutation
        self.haplogroup=haplogroup
        self.flag=flag
        self.treesnp=treesnp
        self.marker=marker
    def __str__(self):
        print(self.name+'\t'+self.chromosome+'\t'+self.position+'\t'+self.genotype+'\t'+self.haplogroup+'\t'+self.treesnp+'\t'+self.mutation+'\t'+self.flag)

def getsnpdb():
    '''
    获取突变数据库
    filename:突变数据库的信息 ep:A102,R,FGC8634; Y3519,15658416,G->A
    :param filename:
    :return:
    '''
    filename='/Users/similarface/Documents/SNPindexY.csv'
    ysnpdb={}
    reader = csv.reader(file(filename, 'rb'))
    for line in reader:
        #key:mark的值
        ysnpdb[line[0]]={'haplogroup':line[1],'othername':line[2],'ypos':line[3],'mutation':line[4]}
    return ysnpdb

class SlugY:
    def __init__(self,rsid,chromosome,position,genotype):
        self.rsid=rsid
        self.chromsome=chromosome
        self.position=position
        self.genotype=genotype

def getSlugYList(file):
    sluglist=[]
    with open(file) as slugfile:
        for line in slugfile.readlines():
            if line.startswith("#"):
                pass
            else:
                lines=line.split('\t')
                if lines[1]=='Y':
                    sluglist.append(SlugY(lines[0],lines[1],lines[2],lines[3].strip()))

def getslugdictFromFile(file):
    slugdict={}
    with open(file) as slugfile:
        for line in slugfile.readlines():
            if line.startswith("#"):
                pass
            else:
                lines=line.split('\t')
                if lines[1]=='Y':
                    slugdict[lines[2]]={'rsid':lines[0],'genotype':lines[3].strip()}
    return slugdict

def getConndb():
    return MySQLdb.connect(host='192.168.30.252',port=3306,db='mfsourcedb',user='dna',passwd='dna',charset='utf8').cursor()

# def gethaplogroup(startpos,ysnpdb,oneman,db,weight):
#     SQL="select id,pid,snp,levels,markers from T_Y_TREE where pid=%s"%(startpos)
#     result=db.execute(SQL)
#     if result==0:
#         return startpos
#     else:
#         #allredict=defaultdict(int)
#         for line in db.fetchall():
#             allow=True
#             id=line[0]
#             pid=line[1]
#             snp=line[2]
#             levels=line[3]
#             markers=[i.strip() for i in line[4].split(',')]
#             markerslen=len(markers)
#             nomatch=0
#             match=0
#             defmatch=0
#             if id==1925:
#                 print(id)
#             for marker in markers:
#                 try:
#                     ypos,mutation=ysnpdb[marker]['ypos'],ysnpdb[marker]['mutation']
#                     info=oneman[ypos]
#                     if info['genotype']==mutation[-1]:
#                         #allredict[line[2]]=allredict[line[2]]+1
#                         #print(line[2],marker,"YES",ypos)
#                         match=match+1
#                     elif info['genotype']!=mutation[-1] and info['genotype']!='--':
#                         #print(line[2],marker,'NO',ypos)
#                         nomatch=nomatch+1/(markerslen*0.01)
#                         if nomatch/markerslen > 0.2:
#                                 allow=False
#                                 break
#                 except Exception,e:
#                     if match>0:
#                         defmatch=defmatch+0.01
#                     else:
#                         defmatch=defmatch+0.5
#
#             weight[id]=match-nomatch-defmatch+weight[pid]+levels*0.7
#             if allow:
#                 #print(id,snp,levels)
#                 a=gethaplogroup(id,ysnpdb,oneman,db,weight)
#             else:
#                 print(id)

def gethaplogroup2(startpos,ysnpdb,oneman,db,weight):
    #查询Y-tree
    SQL="select id,pid,snp,levels,markers from T_Y_TREE where pid=%s"%(startpos)
    #受到印象的行数
    result=db.execute(SQL)
    if result==0:
        return startpos
    else:
        #遍历数据库查询
        for line in db.fetchall():
            allow=True
            id=line[0]
            pid=line[1]
            treesnp=line[2]
            levels=line[3]
            markers=[i.strip() for i in line[4].split(',')]
            #markers的长度
            markerslen=len(markers)
            #不匹配的marker
            nomatch=0
            #真实的不匹配数目
            relnomatch=0
            #匹配的marker
            match=0
            #匹配的markers
            matchmarks=[]
            #不匹配的markers
            nomatchmarks=[]
            #没有找到的marker
            defmatch=0
            #真实的没有找到的marker
            reldefmatch=0
            for marker in markers:
                try:
                    #获取突变数据库ypos的的突变
                    ypos,mutation=ysnpdb[marker]['ypos'],ysnpdb[marker]['mutation']
                    #获取对应ypos 芯片数据的snp
                    snp=oneman[ypos]
                    #如果snp突变了
                    if snp.genotype[0]==mutation[-1]:
                        #突变值＋1
                        match=match+1
                        matchmarks.append(SnpCompare(snp.name,snp.chromosome,snp.position,snp.genotype,ysnpdb[marker]['haplogroup'],treesnp,mutation,'YES',marker))

                    #如果对应ypos上的值不等于
                    elif snp.genotype[0]!=mutation[-1] and snp.genotype!='--':
                        nomatchmarks.append(SnpCompare(snp.name,snp.chromosome,snp.position,snp.genotype,ysnpdb[marker]['haplogroup'],treesnp,mutation,'NO',marker))
                        #如果不匹配 不匹配的marker数目递增
                        nomatch=nomatch+1/(markerslen*0.01)
                        #如果不匹配率大于20% 中断匹配
                        relnomatch=relnomatch+1
                        # if nomatch/markerslen > 0.9:
                        #         allow=False
                        #print(snp.position,snp.genotype[0],mutation,treesnp,id,pid,"NO")
                except Exception,e:
                    #没有找到的marker
                    #发现匹配的每次+0.01
                    reldefmatch=reldefmatch+1
                    if match>0:
                        defmatch=defmatch+0.01
                    else:
                        defmatch=defmatch+0.5

            #权重算法： 匹配的-不匹配的-默认不知道的+加上父级的权重+层次＊0.7
            weight[id]={'weight':match-nomatch-defmatch+weight[pid]['weight']+levels*0.7,'id':id,'pid':pid,'match':match,'nomatch':relnomatch,'dfmatch':reldefmatch,'marklen':markerslen,'matchmarks':matchmarks,'nomatchmarks':nomatchmarks,'snpontree':treesnp}
            if allow:
                #print(id,snp,levels)
                a=gethaplogroup2(id,ysnpdb,oneman,db,weight)
            else:
                print(id)

if __name__ == "__main__":
    #getSlugYList(file)
    #getsnpdb(snpydb)
    ysnpdb=getsnpdb()
    #oneman=snp.reader_dict_ypos('/Users/similarface/Documents/reslut1/genome_wang_mian_Full_20151021231213.txt')
    #oneman=snp.reader_dict_ypos('/Users/similarface/Documents/reslut1/genome_wang_jun_chen_Full_20151021231123.txt')
    #oneman=snp.reader_dict_ypos('/Users/similarface/Downloads/I2b1.txt')
    #oneman=snp.reader_dict_ypos('/Users/similarface/Downloads/H6a1b2.txt')
    #oneman=getslugdictFromFile('/Users/similarface/Downloads/1718.23andme.993')
    #oneman=snp.reader_dict_ypos('/Users/similarface/Downloads/1718.23andme.993')
    #oneman=snp.reader_dict_ypos('/Users/similarface/Downloads/4763.23andme.3353')
    oneman=snp.reader_dict_ypos('/Users/similarface/Documents/wegeneRawData.txt')


    #J2b2
    #
    #oneman=getslugdictFromFile('/Users/similarface/Downloads/R1b1b2a1a1.txt')
    db=getConndb()
    weight={1:{'weight':0,'id':1,'pid':0,'match':0,'nomatch':0,'dfmatch':0,'matchmarks':[],'nomatchmarks':[],'snpontree':''}}
    gethaplogroup2(1,ysnpdb,oneman,db,weight)
    db.close()
    id, max_weight = 100, -100000
    for k,v in weight.items():
        if v['weight'] > max_weight:
            max_weight = v['weight']
            id=k
    kk={}
    pid=id
    while True:
        try:

            print(pid,weight[pid]['weight'],str(weight[pid]['match']),str(weight[pid]['nomatch']),weight[pid]['match'],weight[pid]['snpontree'])
            # if weight[pid]['matchmarks']:
            #     print('-----matchmarks begin-------')
            #     for i in weight[pid]['matchmarks']:
            #         print(i.chromosome,i.flag,i.genotype,i.haplogroup,i.mutation,i.name,i.position,i.treesnp)
            #     print('-----matchmarks end-------')
            # if weight[pid]['nomatchmarks']:
            #     print('-----nomatchmarks begin-------')
            #     for j in weight[pid]['nomatchmarks']:
            #         print(i.chromosome,i.flag,i.genotype,i.haplogroup,i.mutation,i.name,i.position,i.treesnp)
            #     print('-----nomatchmarks end-------')
            pid=weight[pid]['pid']
            kk=weight[pid]
        except KeyError,e:
            break