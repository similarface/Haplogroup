#coding:utf-8
__author__ = 'similarface'

rsiddict={}
for i in open('/Users/similarface/Downloads/rsidref.tab','r'):
    lines=i.strip().split('\t')
    rsiddict[lines[0]]=[lines[1],lines[2]]

def phrse(gen,genotype):
    if gen=='A':
        if gen in genotype:
            return True
        elif genotype=='TT':
            return True
        else:
            return False
    if gen=='T':
        if gen in genotype:
            return True
        elif genotype=='AA':
            return True
        else:
            return False
    if gen=='G':
        if gen in genotype:
            return True
        elif genotype=='CC':
            return True
        else:
            return False
    if gen=='C':
        if gen in genotype:
            return True
        elif genotype=='GG':
            return True
        else:
            return False


#
# unfond=0
# found=0
# match=0
# unmatch=0
# rrssiidd=[]
# all=0
# for j in open('/Users/similarface/Downloads/cc','r'):
#     all=all+1
#     line=j.strip().split('\t')
#     try:
#         if phrse(rsiddict[line[2]][0],line[3]):
#             print(j.strip()+'\t'+rsiddict[line[2]][0]+'\t'+rsiddict[line[2]][1]+'\t'+'YES')
#             match=match+1
#         else:
#             print(j.strip()+'\t'+rsiddict[line[2]][0]+'\t'+rsiddict[line[2]][1]+'\t'+'NO')
#             unmatch=unmatch+1
#         found=found+1
#     except Exception,e:
#         print(j.strip()+'\t'+'?'+'\t'+'UNKOWN')
#         unfond=unfond+1
#         rrssiidd.append(line[2])
# print("总数:"+str(all)+'\t'+'找到:'+str(found)+'\t'+'未找到:'+str(unfond)+'\t'+'匹配:'+str(match)+'\t'+'不配:'+str(unmatch))
# print(rrssiidd)


nowdict={}
for  kk in open('/Users/similarface/Downloads/kk','r'):
    linesx=kk.strip().split('\t')
    rsidx=linesx[2]
    v1=linesx[3].split(':')[1].strip()
    nowdict[rsidx]=v1

for li in open('/Users/similarface/Downloads/ddd','r'):
    ref=li.strip().split('\t')[4]
    rsid=li.strip().split('\t')[2]
    try:
        if phrse(ref,nowdict[rsid]):
            print(li.strip()+'\t'+'YESYES')
        else:
            print(li.strip()+'\t'+'NONO')
    except Exception,e:
        print(li.strip()+'\t'+'UNKOWN')