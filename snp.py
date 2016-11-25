#coding:utf-8
__author__ = 'similarface'
import sys

def reader(infile):
    '''
    读取文件返回list
    :param infile: 
    :return:
    '''
    handle = open(infile,"r")
    data = []
    for x in handle:
        if x[0] != "#":
            single = x.split("\t")
            data.append(snp(single[0],single[1],single[2],single[3].rstrip()))
    handle.close()
    return data

def reader_dict(infile):
    '''
    读取文件获取字典
    :param infile: 
    :return:
    '''
    handle = open(infile,"r")
    data = {}
    for x in handle:
        if x[0] != "#":
            single = x.split("\t")
            #KEY=rsid
            # data[single[0]] = snp(single[0],single[1],single[2],single[3].rstrip())
    handle.close()
    return data

def reader_dict_pos(infile):
    '''
    读取文件获取字典 
    字典的key是染色体位置
    :param infile: 
    :return:
    '''
    handle = open(infile,"r")
    data = {}
    for x in handle:
        if x[0] != "#":
            single = x.split("\t")
            #KEY=rsid
            data[single[1]+":"+single[2]] = snp(single[0],single[1],single[2],single[3].rstrip())
    handle.close()
    return data

def reader_dict_ypos(infile):
    '''
    读取文件获取Y染色体字典 
    字典的key是Y染色体位置
    :param infile: 
    :return:
    '''
    handle = open(infile,"r")
    data = {}
    for x in handle:
        if x[0] != "#":
            single = x.split("\t")
            if single[1]=='Y':
                data[single[2]] = snp(single[0],single[1],single[2],single[3].rstrip())
    handle.close()
    return data

def reader_dict_xpos(infile):
    '''
    读取文件获取X染色体字典 
    字典的key是X染色体位置
    :param infile: 
    :return:
    '''
    handle = open(infile,"r")
    data = {}
    for x in handle:
        if x[0] != "#":
            single = x.split("\t")
            if single[1]=='X':
                data[single[2]] = snp(single[0],single[1],single[2],single[3].rstrip())

def reader_dict_mtpos(infile):
    '''
    读取文件获取MT染色体字典 
    字典的key是MT染色体位置
    :param infile: 
    :return:
    '''
    handle = open(infile,"r")
    data = {}
    for x in handle:
        if x[0] != "#":
            single = x.split("\t")
            try:
                if single[1]=='MT'  and single[3].strip()!='--':
                    data[single[2]] = snp(single[0],single[1],single[2],single[3].rstrip())
            except IndexError,e:
                print(x)
    return data



class snp():
    '''
    name,==rsid
    chromosome,==染色体
    position==位置
    genotype==基因型
    '''
    def __init__(self,name,chromosome,position,genotype):
        self.name = name
        self.chromosome = chromosome
        self.position = position
        self.genotype = genotype

    def __str__(self):
        print(self.name+'\t'+self.chromosome+'\t'+self.position+'\t'+self.genotype)

#reader(sys.argv[1])