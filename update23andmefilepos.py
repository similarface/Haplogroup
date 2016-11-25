#coding:utf-8
__author__ = 'similarface'
import os,sys
import tempfile
import shutil
import MySQLdb

def getConnection():
    '''
    连接数据库
    :return:
    '''
    try:
        db=MySQLdb.connect(host='192.168.30.252',port=3306,db='gendb',user='dna',passwd='dna',charset='utf8')
        cursor = db.cursor()
    except MySQLdb.Error, e:
        error_code = e.args[0]
        error_msg = 'MySQL error! ', e.args[0], e.args[1]
        print error_msg
    return cursor


def update23andmeFilePos(f_23andmefile):
    conn=getConnection()
    '''
    使用文件名称作为samplename的名称
    :param vcffile:
    :return:
    '''
    #临时文件
    temp_file=tempfile.mktemp(dir='/tmp')
    file=open(temp_file,'w+b')
    if os.path.exists(f_23andmefile):
        fopen=open(f_23andmefile,'r')
    else:
        print "file %s not found" % f_23andmefile
        sys.exit()
    for line in fopen:
        if line.startswith("#"):
            file.write(line)
        else:
            try:
                lines=line.split('\t')
                if lines[1] not  in ['Y','MT']:
                    pass
                else:
                    pos=getPosFromRsid(conn,lines[0])
                    if pos:
                        file.write('\t'.join([lines[0],lines[1],pos,lines[3]]))
                    else:
                        pos=getPosFromiRsid(conn,lines[0])
                        if pos:
                            file.write('\t'.join([lines[0],lines[1],pos,lines[3]]))
                        else:
                            file.write(line)
            except Exception,e:
                file.write(line)

    fopen.close()
    file.close()
    if os.path.exists(f_23andmefile):
        os.remove(f_23andmefile)
    shutil.copy(temp_file,f_23andmefile)
    try:
        os.remove(temp_file)
    except OSError:
        pass
    finally:
        conn.close()

def getPosFromRsid(conn,rsid):
    SQL="select pos from T_DBSNP_HG19_138 where rsid='%s' " %(rsid)
    conn.execute(SQL)
    try:
        results = conn.fetchone()
        if results!=None and results!="":
            return results[0]
        else:
            return None
    except Exception,e:
        return None

def getPosFromiRsid(conn,irsid):
    SQL="select position from t_23andme_addrefs where rsid='%s' " %(irsid)
    conn.execute(SQL)
    try:
        results = conn.fetchone()
        if results!=None and results!="":
            return results[0]
        else:
            return None
    except Exception,e:
        return None


if __name__=='__main__':
    update23andmeFilePos('/Users/similarface/Documents/wocao.txt')