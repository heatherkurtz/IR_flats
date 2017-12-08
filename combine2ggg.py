import combineggg as combine
import sys,os,tempfile
from mx import DateTime
from numpy import *
from threading import Thread

props = "11528,11624,12051,12005,11709,11738,11602,11663,12064,12197,12224,12203,12332,11696,12184,12099,12307,12329,12065,12061,12068,12286,12283,12167,11591.12328,11108,11142,11149,11153,11166,11189,11202,11208,11343,11359,11519,11520,11534,11541,11557,11563,11584,11597,11600,11644,11650,11666,11669,11694,11700,11702,11735,11838,11840,11587"
#props = "11108,11142,11149,11153,11166,11189,11202,11208,11343,11359,11519,11520,11534,11541,11557,11563,11584,11597,11600,11644,11650,11666,11669,11694,11700,11702,11735,11838,11840,11587"
props = props +",12108,12106,12616,12453,12076,12286,12070,12440,12460,12581"
props = ""
def irange(sequence):
# Returns a i,object list. i.e. for i,object in irange(list)
# produces both an iteration/rank number and the member object itself
      return map(None,range(len(sequence)),sequence)

def get_filters(start,end,exptime=300):
    import MySQLdb
    host = "asterix.stsci.edu"
    database = "wfc3"
    global props
    

    
    s = "select distinct(filter) from sky where sun_alt is not null and filter like \"f%s\" and proposid in (%s) and (sky-mode)/sky<0.05 and exptime>%f and start>=%f and end<%f" % ("%",props,exptime,start,end)
    db=MySQLdb.connect(host=host,db=database,user="npirzkal",passwd="medusa")
    c=db.cursor()
    c.execute(s)
    #print s
    lines = c.fetchall()
    c.close()
    db.close()
    
    filters = [l[0] for l in lines]
    return filters
    
    
def get_datasets(start,end,exptime,filter=None):
    import MySQLdb
    host = "asterix.stsci.edu"
    database = "wfc3"
    global props

    if props!="":
        ptxt = " and proposid in (%s) " % (props)
    else:
        ptxt = ""

 #   if filter!="ALL":
 #       s = "select dataset from sky where sun_alt is not null and filter=\"%s\" and proposid in (%s) and (sky-mode)/sky<0.05 and exptime>%f and start>=%f and end<%f" % (filter,props,exptime,start,end)
    if filter=="ALL":
        s = "select dataset from sky where sun_alt is not null and filter in (%s) %s and (sky-mode)/sky<0.05 and exptime>%f and start>=%f and end<%f" % ("\"F098M\",\"F105W\",\"F110W\",\"F125W\",\"F160W\"",ptxt,exptime,start,end)
    elif filter=="NOF160W":
        s = "select dataset from sky where sun_alt is not null and filter in (%s) %s and (sky-mode)/sky<0.05 and exptime>%f and start>=%f and end<%f" % ("\"F098M\",\"F105W\",\"F110W\",\"F125W\"",ptxt,exptime,start,end)
    else:
        s = "select dataset from sky where sun_alt is not null and filter=\"%s\" %s and (sky-mode)/sky<0.05 and exptime>%f and start>=%f and end<%f" % (filter,ptxt,exptime,start,end)

    print s
    
    db=MySQLdb.connect(host=host,db=database,user="npirzkal",passwd="medusa")
    c=db.cursor()
    c.execute(s)
    #print s
    lines = c.fetchall()
    c.close()
    db.close()
    datasets = [ll[0] for ll in lines]
    return datasets

def load_list2(f,datasets,start=0,end=1e32):
    lines = open(f).readlines()
    fs = []
    for l in lines:
        ws = l.split()
        expstart = float(ws[0])
        filename = ws[1]
        
        dataset = os.path.split(filename)[-1].split("_flt.fits")[0]
        if dataset not in datasets:
            continue
       # print "hellO"
        camera = ws[5]
        target = ws[4]
        filt = ws[8]
        proposid = ws[6]
        #print camera,target,filt,proposid
        if camera!="IR":
            continue
        if target!="NONE":
            continue
        if expstart<start or expstart>end:
            continue
        fs.append(filename.split("[")[0])
    return fs
    
def load_list(datasets,start=0,end=1e32):
    import sqlite3
    import os
    #conn = sqlite3.connect("/grp/hst/wfc3a/Database/ql_072111.db")
    conn = sqlite3.connect("/grp/hst/wfc3a/Database/ql.db")
    c = conn.cursor()
    s = "select m.dir,m.filename,ir.detector,ir.targname,ir.filter,ir.proposid,ir.expstart,ir.exptime from master as m, IR_FLT_0 as ir  where m.id=ir.id and ir.detector=\"IR\" "

    c.execute(s)
    fs = []
    starts = []
    for row in c:
        filename = os.path.join(row[0],row[1])
	dataset = row[1].split("_flt.fits")[0]
	if dataset not in datasets:
		continue
        camera = "IR"
        target = row[3]
        filt = row[4]
        proposeid = int(row[5])
        expstart = float(row[6])
        exptime = float(row[7])

        if camera!="IR":
            print "Not IR"
            continue
        if target=="NONE":
            print "No target"
            continue
        if expstart<start or expstart>end:
            print "No start"
            continue
        if filt[0]!="F":
            print "No F"
            continue
        fs.append(filename)
        starts.append(expstart)
    return fs,starts



def histo(filter,ndays=7,date0 = "2009-07-17"):
    start_date = map(int,date0.split("-"))
    print start_date
    d = DateTime.Date(start_date[0],start_date[1],start_date[2])
    start0 = d.mjd
    today = DateTime.today().mjd
    

    for d in range(start0,today,ndays):
        start = d
        end = start + ndays
        datasets = get_datasets(start,end,exptime=300.,filter=filter)
        e = DateTime.DateTimeFromMJD(end)
        end_date = "%s-%s-%s" % (e.year,e.month,e.day)
        e = DateTime.DateTimeFromMJD(start)
        start_date = "%s-%s-%s" % (e.year,e.month,e.day)
        print start_date,end_date,len(datasets)
        if len(datasets) < 10:
            continue
#        files = load_list("/tmp/all4.lst",datasets,start=start,end=end)
        files = load_list(datasets,start=start,end=end)

        res = []
        for f in files:
            data = combine.run_sextractor(f,badval = 32,start=start,end=end)
            res.append(data)
        print "Combining...",
        newdata,newdata2,newdata3 = combine.combine(res,badpix=32)
        print "Done."
        combine.write_fits("s_%s-%s.%s.fits" % (start_date,end_date,filter),newdata)

def quad_stats_all(fs):
    import glob,Statistics
    
    files = glob.glob(fs)
    skies = []
    for f in files:
        sky = quad_stats(f)
        skies.append(sky)
        
    skies = asarray(skies)
    ss = transpose(skies)

    
    print "Q1:",Statistics.average(ss[0]),Statistics.standardDeviation(ss[0])
    print "Q2:",Statistics.average(ss[1]),Statistics.standardDeviation(ss[1])
    print "Q3:",Statistics.average(ss[2]),Statistics.standardDeviation(ss[2])
    print "Q4:",Statistics.average(ss[3]),Statistics.standardDeviation(ss[3])
    return skies
    
def quad_stats(f):
    """Compute the mode and median of each quads in file f"""
    import pyfits,Statistics,combine
    from mx import DateTime
    
    d0 = f.split("-")
    d = DateTime.Date(int(d0[0][2:]),int(d0[1]),int(d0[2]))
    start0 = d.mjd
    
    fin = pyfits.open(f)
    pad = 10
    d1 = ravel(fin[0].data[507+pad:1015-pad,0+pad:508-pad])  
    d2 = ravel(fin[0].data[507+pad:1015-pad,507+pad:1015-pad])  
    d3 = ravel(fin[0].data[0+pad:508-pad,507+pad:1015-pad])  
    d4 = ravel(fin[0].data[0+pad:508-pad,0+pad:508-pad])  
    
    sky1 = Statistics.median(d1)
    mode1,m11,m12 = combine.mode_conf(d1)

    sky2 = Statistics.median(d2)
    mode2,m21,m22 = combine.mode_conf(d2)

    sky3 = Statistics.median(d3)
    mode3,m31,m32 = combine.mode_conf(d3)

    sky4 = Statistics.median(d4)
    mode4,m41,m42 = combine.mode_conf(d4)

    print f,start0,sky1,sky2,sky3,sky4 
    return sky1,sky2,sky3,sky4



def get_persistence_mask(f,maxval=30000.,numdayback=1.):
    import MySQLdb,pyfits
    
    print "Looking for persistence"
    
    fin = pyfits.open(f)
    data = fin[1].data
    exptime = fin[0].header["EXPTIME"]
    start = fin[0].header["EXPSTART"]
    fin.close()
    
    host = "asterix.stsci.edu"
    database = "wfc3"
    
    s = "select dataset from sky where start>=%f and end<%f" % (start-numdayback,start)
    db=MySQLdb.connect(host=host,db=database,user="npirzkal",passwd="medusa")
    c=db.cursor()
    c.execute(s)
    lines = c.fetchall()
    c.close()
    db.close()
    datasets = [ll[0] for ll in lines]
    
    print datasets
#    files2 = load_list("/tmp/all4.lst",datasets)
    files2,tt = load_list(datasets)

    print files2
    
    mask = zeros(shape(data),'b')
    for f in files2:
        print "Opening:",f
        fin2 = pyfits.open(f)
        data = fin2[1].data
        exptime = fin2[0].header["EXPTIME"]
        end = fin2[0].header["EXPEND"]
        filter = fin2[0].header["FILTER"]
        if filter=="Blank": continue
        vg = (data*exptime)>maxval
        print shape(data),shape(vg)
        mask[vg] = 1
        print "Adding ",len(ravel(mask[vg])), " pixels affected from",fin2[0].header["FILTER"],fin2[0].header["TARGNAME"]," Total:",len(ravel(mask[mask==1]))
        fin2.close()
    return mask
        
        
        
if __name__=="__main__":
    import sys,futils,pyfits
    import time

    cache = 1
    if len(sys.argv)==5:
        cache = 0
    
    
    if os.path.isfile(sys.argv[1]):
        res = open(sys.argv[1]).readlines()[0:250]
        end_date = ""
        filter = ""
    else:
        start_date = map(int,sys.argv[1].split("-"))
        print start_date
        d = DateTime.Date(start_date[0],start_date[1],start_date[2])
        start = d.mjd

        end_date = map(int,sys.argv[2].split("-"))
        print end_date
        d = DateTime.Date(end_date[0],end_date[1],end_date[2])
        end = d.mjd

        days = end - start

        #days = float(sys.argv[2])
        #end = start + days
        #ndays = days
        
        #e = DateTime.DateTimeFromMJD(end)
        #end_date = "%s-%s-%s" % (e.year,e.month,e.day)
        print start_date,end_date, days
        ndays = days
                
        filter = sys.argv[3]
            
        MINTIME = 600


        
        datasets = get_datasets(start,end,MINTIME,filter=filter)
        print "datasets:",len(datasets)
        #print datasets
#        files = load_list("/tmp/all4.lst",datasets)
        files,starts = load_list(datasets)

        print "files:",len(files)
        start = min(starts)
        e = DateTime.DateTimeFromMJD(start)
        start_date = "%s-%s-%s" % (e.year,e.month,e.day)

        end = max(starts)
        e = DateTime.DateTimeFromMJD(end)
        end_date = "%s-%s-%s" % (e.year,e.month,e.day)
        print start_date,end_date
        
        res = []
        for i,f in irange(files):
            dataset = os.path.split(f)[-1].split("_flt.fits")[0]
            if not os.path.isfile(f):
                print f,"not found"
                continue
            print i," of ",len(files),":",f,dataset,pyfits.open(f)[0].header["FILTER"],pyfits.open(f)[0].header["PFLTFILE"]
            already_done = combine.check_if_done(dataset,updated=55767.0)
            OLDPFILE = pyfits.open(f)[0].header["PFLTFILE"].split("iref$")[-1]
            OLDPFILE = os.path.join("/grp/hst/cdbs/iref/",OLDPFILE)
            already_done = [already_done[0],already_done[1],already_done[2],OLDPFILE]


            if already_done[0] != None:
                print dataset, "is cached"
                if cache==1:
                    res.append(already_done)
                    continue
                else:
                    print "But I do not care"
            try:
            #if 1:
                OLDPFILE = pyfits.open(f)[0].header["PFLTFILE"].split("iref$")[-1]
                OLDPFILE = os.path.join("/grp/hst/cdbs/iref/",OLDPFILE)

                persist_mask = get_persistence_mask(f,maxval=30000.,numdayback=3.)
                #persist_mask = get_persistence_mask(f,maxval=30000.,numdayback=0.)

                ##name = tempfile.mktemp()
                data,allsum = combine.run_sextractor(f,badval = 32,start=start,end=end, persist_mask=persist_mask,logit=1)
		if data == None:
			continue                
#                print "Before:",len(ravel(data[data==32]))
#                data[persist_mask>0] = 32
#                print "After:",len(ravel(data[data==32]))

                combine.write_fits("%s-%s.%s.allsum.2fff.fits" % (sys.argv[1],end_date,filter),allsum)


                ##combine.write_fits(name,data)

                CACHEPATH = "/Users/npirzkal/WFC3/combine/cache"
                oname = os.path.split(f)[-1]
                oname = os.path.join(CACHEPATH,oname)
                print oname

                if cache==1:
                    combine.write_fits(oname,data)
                    ngoods = len(ravel(data[data!=32]))
                    combine.update_cache_info(f,oname,ngoods)
                    
                already_done = combine.check_if_done(dataset,updated=0)
                already_done = [already_done[0],already_done[1],already_done[2],OLDPFILE]
                if already_done[0]!=None:
                    res.append(already_done)
            except IOError:
                print f," went away!"
                

#print "Done for now"
#   sys.exit(1)
    
    print res

    s = []
    for r in res:
        s.append("%s %f %f\n" % (r[0],r[1],r[2]))

    open("combine2fff.%s-%s.%s.lst" % (sys.argv[1],end_date,filter),"w").writelines(s)
    print "Combining...",
    
    #modedata,avgdata,ndata,stdevdata,meddata = combine.combine_pixel(res,badpix=32)
    cubename = "%s-%s.%d.%s.%d.3dcube.fits" % (sys.argv[1],end_date,ndays,filter,MINTIME)
    ndata,avgdata,stdevdata,meddata,allsum = combine.combine_pixel(res,badpix=32,cubename=cubename)

    print "Done."
    
    
    for f in res:
        if os.path.isfile(f[0]):
            print f[0]
            #print "Cleaning ",f
            #os.unlink(f)
                
    #combine.write_fits("%s-%s.%s.mode.2f.fits" % (sys.argv[1],end_date,filter),modedata)
    combine.write_fits("%s-%s.%d.%s.%d.allsum2.2fff.fits" % (sys.argv[1],end_date,ndays,filter,MINTIME),allsum)

    combine.write_fits("%s-%s.%d.%s.%d.avg.2fff.fits" % (sys.argv[1],end_date,ndays,filter,MINTIME),avgdata)
    combine.write_fits("%s-%s.%d.%s.%d.n.2fff.fits" % (sys.argv[1],end_date,ndays,filter,MINTIME),ndata)
    combine.write_fits("%s-%s.%d.%s.%d.stdev.2fff.fits" % (sys.argv[1],end_date,ndays,filter,MINTIME),stdevdata)
    combine.write_fits("%s-%s.%d.%s.%d.med.2fff.fits" % (sys.argv[1],end_date,ndays,filter,MINTIME),meddata)
    
