__author__ = 'similarface'
from collections import defaultdict
import glob
ip = r"?P<ip>[\d.]*"
date = r"?P<date>\d+"
month = r"?P<month>\w+"
year = r"?P<year>\d+"
log_time = r"?P<time>\S+"
method = r"?P<method>\S+"
request = r"?P<request>\S+"
status = r"?P<status>\d+"
bodyBytesSent = r"?P<bodyBytesSent>\d+"
refer = r"""?P<refer>
         [^\"]*
         """
userAgent=r"""?P<userAgent>
            .*
           """
import re

import linecache
def readline(path):
    return linecache.getlines(path)

ipad=0
result=defaultdict(int)
p = re.compile(r"(%s)\ -\ -\ \[(%s)/(%s)/(%s)\:(%s)\ [\S]+\]\ \"(%s)?[\s]?(%s)?.*?\"\ (%s)\ (%s)\ \"(%s)\"\ \"(%s).*?\"" %( ip, date, month, year, log_time, method, request, status, bodyBytesSent, refer, userAgent ), re.VERBOSE)

acctype=set()
linuxcount=0
Androidcount=0
shoujidict=defaultdict(int)
for line in open('/data4/access.log','r'):
    m = re.findall(p, line)
    try:
        if len(m)!=0:
            #"Mozilla/5.0 (Linux; U; Android 4.4.4; zh-cn; SM-G5108Q Build/KTU84P)
            #(Linux; U; Android 4.4.2; zh-cn; HUAWEI MT7-CL00 Build/HuaweiMT7-CL00) AppleWebKit/533.1 (KHTML, like Gecko)Version/4.0 MQQBrowser/5.4 TBS/02544
            result[m[0][10]]=result[m[0][10]]+1

            try:
                content=m[0][10]
                begin=content.index('(')
                end=content.index(')')
                shortcontent=content[begin+1:end]
                acctype.add(shortcontent.split(';')[0])
                #print(shortcontent)
                if 'Linux' in shortcontent:
                    linuxcount=linuxcount+1
                    shouji=shortcontent.split(';')[4]
                    shoujidict[shouji]=shoujidict[shouji]+1
                if 'Windows' in shortcontent:
                    pass
                if 'iPad' == shortcontent.split(';')[0]:
                    ipad=ipad+1


            except Exception,e:
                pass

    except Exception,e:
        #print(m)
        pass

fopen=open('/tmp/shouji.txt','w')
for k,v in shoujidict.items():
    fopen.writelines(k+'\t'+str(v)+"\n")
fopen.close()

print(ipad)