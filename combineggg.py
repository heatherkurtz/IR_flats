import pyfits,os,sys
import Statistics,Stats2
from numpy import *
import sextractor
from threading import Thread
import threading
from mx import DateTime

ALLSUM = None


#OLDPFLATS = {"F098M":"sca2025ti_pfl.fits","F105W":"sca20260i_pfl.fits","F110W":"sca20261i_pfl.fits","F125W":"sca20262i_pfl.fits","F140W":"sca20269i_pfl.fits","F160W":"sca2026bi_pfl.fits"}
NEWPFLATS = {}
NEWPFLATS["F098M"]="uc72113ni_pfl.fits"
NEWPFLATS["F105W"]="uc72113oi_pfl.fits"
NEWPFLATS["F110W"]="uc72113pi_pfl.fits"
NEWPFLATS["F125W"]="uc72113qi_pfl.fits"
NEWPFLATS["F126N"]="uc72113ri_pfl.fits"
NEWPFLATS["F127M"]="uc72113si_pfl.fits"
NEWPFLATS["F130N"]="uc721140i_pfl.fits"
NEWPFLATS["F132N"]="uc721141i_pfl.fits"
NEWPFLATS["F139M"]="uc721142i_pfl.fits"
NEWPFLATS["F140W"]="uc721143i_pfl.fits"
NEWPFLATS["F153M"]="uc721144i_pfl.fits"
NEWPFLATS["F160W"]="uc721145i_pfl.fits"
NEWPFLATS["F164N"]="uc721146i_pfl.fits"
NEWPFLATS["G141"]="u4m1335mi_pfl.fits"



def log_sky(filename,sky,mode,ngoods,start,end):
    """Record the dataset name, ra, dec and sky, exptime, and sky value in a SQL table"""
    import ephem
    import MySQLdb
    host = "asterix.stsci.edu"
    database = "wfc3"
    db=MySQLdb.connect(host=host,db=database,user="npirzkal",passwd="medusa")

    hdr = pyfits.open(filename)[0].header
    dataset = hdr["ROOTNAME"].strip()
    exptime = hdr["EXPTIME"]
    filter = hdr["FILTER"].strip()
    ra = hdr["RA_TARG"]
    dec = hdr["DEC_TARG"]
    #print '%s' % (ra),'%s' % (dec)
    p = ephem.Equatorial('%s' % (ra),'%s' % (dec))
    q = ephem.Ecliptic(p)
    lon = q.lon/(2*pi)*360.
    lat = q.lat/(2*pi)*360.
    proposid = hdr["PROPOSID"]
    sun_alt = hdr["SUN_ALT"]
    start = hdr["EXPSTART"]
    end = hdr["EXPEND"]
    
#SUBTYPE
#DETECTOR
    s = "delete from sky where dataset like \"%s\" and ext=\"FLT\"" % (dataset)
    c=db.cursor()
    c.execute(s)
    #print s
    lines = c.fetchall()
    c.close()
    print "sky",sky,mode
    
    now = DateTime.now().mjd
    id = dataset + "_"

    if sky==None or mode==None:
        sky = "NULL"
        mode = "NULL"
        print "sky",id,dataset,filter,ra,dec,long,lat,exptime,proposid,ngoods,start,end,sun_alt,now
        s = "insert into %s (id,dataset,filter,ra0,dec0,long0,lat0,exptime,sky,mode,proposid,ngoods,start,end,sun_alt,updated) values(\"%s\",\"%s\",\"%s\",%12.10f,%12.10f,%12.10f,%12.10f,%f,NULL,NULL,%d,%d,%f,%f,%f,%10.5f)" % ("sky",id,dataset,filter,ra,dec,lon,lat,exptime,proposid,ngoods,start,end,sun_alt,now)

    else:    
        s = "insert into %s (id,dataset,filter,ra0,dec0,long0,lat0,exptime,sky,mode,proposid,ngoods,start,end,sun_alt,updated) values(\"%s\",\"%s\",\"%s\",%12.10f,%12.10f,%12.10f,%12.10f,%f,%f,%f,%d,%d,%f,%f,%f,%10.5f)" % ("sky",id,dataset,filter,ra,dec,lon,lat,exptime,sky,mode,proposid,ngoods,start,end,sun_alt,now)
        
    print s
    c=db.cursor()
    c.execute(s)
    print s
    #lines = c.fetchall()
    c.close()
    db.close()

    
def imcopy(ifile,ext,ofile):
    """Copy ext of ifile into a new file called ofile, copying header too"""
    hdr = pyfits.getheader(ifile,ext)
    data = pyfits.open(ifile)[ext].data

    if os.path.isfile(ofile):
        os.unlink(ofile)

    pyfits.writeto(ofile,data,hdr)
  
def mode_conf(data0,conf=0.682,n=10,skip=0,plotit=0):
    from pygsl import statistics as gslstats
    
    data = ravel(data0[skip:])
    mean = Statistics.average(ravel(data[(data<20)&(data>0)]))
    stddev = Statistics.standardDeviation(ravel(data[(data<20)&(data>0)]))   
    mean = Statistics.average(ravel(data))
    stddev = Statistics.standardDeviation(ravel(data))   
    #print mean, stddev
    #data = data[(data<mean+3*stddev)&(data>mean-3*stddev)]
    if plotit==1:
        pylab.hist(ravel(data),n)
        pylab.show()

    print "data:",data
    num = len(data)
    
    a = sort(data)
    p25 = gslstats.quantile_from_sorted_data(a,0.25)
    p75 = gslstats.quantile_from_sorted_data(a,0.75)
    IQR = p75-p25
    h = 2*IQR/len(a)**(1./3)
    if min(a)==max(a) or h==0.:
        n = 1
    else:
        n = int((max(a)-min(a))/h)+1
    
    
    histo = histogram(data,n)
    vals = histo[0]
    lbins = histo[1][:-1] # left edges of bins # (histo[1][:-1]+histo[1][1:])/2.
    rbins = histo[1][1:] # right edges of bins # (histo[1][:-1]+histo[1][1:])/2.
    mode =  (lbins[vals==max(vals)][0] + rbins[vals==max(vals)][0])/2.

    for nn in arange(max(vals),0,-1):
        #print nn
        data1 =  vals[vals>=nn]
        rbins1 =  rbins[vals>=nn]
        lbins1 =  lbins[vals>=nn]

        s = sum(data1)
        #print nn,s,num*conf
        if s>(num*conf):
            break

    #print "95%:",num,s,min(lbins1),max(rbins1)

    #print len(vals),len(bins)
    a = (lbins>=min(lbins1))
    b = (rbins<=max(lbins1))
    #print "check:",sum(vals),num,sum(vals[a&b]),1.0*sum(vals[a&b])/num
    #if plotit==1:
    return mode,min(lbins1),max(rbins1)

def mode_conf2(data0,conf=0.682,n=10,skip=0,plotit=0):
    import pymc
    
    mode = pymc.utils.hpd(data0,0.995)    
    print "mode:",mode,sum(mode)/2.,pymc.utils.hpd(data0,1.-conf)
    conf = pymc.utils.hpd(data0,1.-conf)
    
    return sum(mode)/2.,conf[0],conf[1]

def check_if_done(dataset,updated=0):
    import MySQLdb
    host = "asterix.stsci.edu"
    database = "wfc3"
    db=MySQLdb.connect(host=host,db=database,user="npirzkal",passwd="medusa")
    c=db.cursor()

    s = "select filename,exptime,sky from cache as c, sky as s where exptime is not null and c.dataset=s.dataset and c.dataset=\"%s\" and s.updated>%f" % (dataset,updated)
    print s
    c.execute(s)
    l = c.fetchall()
    if len(l)==0:
        filename = None
        exptime = None
        sky = None
    else:
        filename = l[0][0]
        exptime = float(l[0][1])
	if l[0][2]!=None:
             sky = float(l[0][2])
        else:
             filename = None
             exptime = None
             sky = None 
    return filename,exptime,sky

def update_cache_info(orgfile,filename,ngoods):        
    import MySQLdb
    import time
    
    host = "asterix.stsci.edu"
    database = "wfc3"
    db=MySQLdb.connect(host=host,db=database,user="npirzkal",passwd="medusa")
    c=db.cursor()

    dataset = os.path.split(filename)[-1].split("_flt.fits")[0]
    print dataset
    s = "delete from cache where dataset=\"%s\"" % (dataset)
    print s
    c.execute(s)
    
    now = time.time()
    s = "insert into cache (orgfile,dataset, filename, time,ngoods) values (\"%s\",\"%s\",\"%s\",%f,%d)" % (orgfile,dataset, filename, now, ngoods)
    print s
    c.execute(s)

    
def run_sextractor(filename,badval = 1e8,logit=1,start=0,end=0,force=0,persist_mask=None,makeitflat=None):
    """Run SExtractor and returns a data array containing data where objects are set to badval"""
    import sextractor2,os
    from pyraf import iraf
    import tempfile,sys
    global ALLSUM
   
    
        
 # ttfile   checkfile check2file check3file check4file 
    ttfile = tempfile.mktemp(".fits") # this one is the one that will be used and returned. No FF applied
    ttfile2 = tempfile.mktemp(".fits") # this one is the one where latest FF is applied (old one de-applied)
    
    #if os.path.isfile("tt.fits"):
    #    os.unlink("tt.fits")
    #iraf.imcopy(filename+"[1]","tt.fits")
    print filename, ttfile
    imcopy(filename,1,ttfile)
    imcopy(filename,1,ttfile2)
    
    print ttfile,ttfile2
    
    ys,xs = shape(pyfits.open(ttfile)[0].data)
    if ys!=1014 or xs!=1014:
        print "Not a full array!"
        return None,None
#        sys.exit(1)
    # Generate persistance mask
    #startime = pyfits.open(filename)[0].header["STARTIME"].split()[0][0:3]

    
    
    calver = pyfits.open(filename)[0].header["CAL_VER"].split()[0][0:3]
    calver = float(calver)
    print calver
    
    if calver>= 1.8:
        print "ALL OK!"
    else:
        print "Correct for GAIN problem"
        fin = pyfits.open(ttfile,mode="update")
        fin[0].data[507:1015,0:508] = fin[0].data[507:1015,0:508] * 1.004
        fin[0].data[507:1015,507:1015] = fin[0].data[507:1015,507:1015] * 0.992
        fin[0].data[0:508,507:1015] = fin[0].data[0:508,507:1015] * 1.017
        fin[0].data[0:508,0:508] = fin[0].data[0:508,0:508] * 0.987
        fin.close()
        fin = pyfits.open(ttfile2,mode="update")
        fin[0].data[507:1015,0:508] = fin[0].data[507:1015,0:508] * 1.004
        fin[0].data[507:1015,507:1015] = fin[0].data[507:1015,507:1015] * 0.992
        fin[0].data[0:508,507:1015] = fin[0].data[0:508,507:1015] * 1.017
        fin[0].data[0:508,0:508] = fin[0].data[0:508,0:508] * 0.987
        fin.close()
        
    
    
    # De-aplying PFLAT and applying current one 
    
    
    FILTER = pyfits.open(filename)[0].header["FILTER"]
    if not NEWPFLATS.has_key(FILTER):
        print "Filter ",FILTER," is unkown. Moving on!"
        return None,None
        
    NEWPFILE = NEWPFLATS[FILTER]
    NEWPFILE = os.path.join("/grp/hst/cdbs/iref/",NEWPFILE)
    NEWPFLAT = pyfits.open(NEWPFILE)[1].data[5:1014+5,5:1014+5]
    
    OLDPFILE = pyfits.open(filename)[0].header["PFLTFILE"].split("iref$")[-1]
    OLDPFILE = os.path.join("/grp/hst/cdbs/iref/",OLDPFILE)
    OLDPFLAT = pyfits.open(OLDPFILE)[1].data[5:1014+5,5:1014+5]
                
    #print OLDPFILE,FILTER,NEWPFILE,OLDPFLAT,NEWPFLAT
    #print shape(OLDPFLAT),shape(NEWPFLAT)
    #print min(ravel(OLDPFLAT)),max(ravel(OLDPFLAT)),min(ravel(NEWPFLAT)),max(ravel(NEWPFLAT))
            
    # De-applying flat field and applying latest flat to ttfile2
    fin = pyfits.open(ttfile2,mode="update")
    #print shape(fin[0].data),shape(OLDPFLAT),shape(NEWPFLAT)
    fin[0].data = fin[0].data * OLDPFLAT / NEWPFLAT
    fin.close()
    # De-applying flat field to ttfile
    fin = pyfits.open(ttfile,mode="update")
    fin[0].data = fin[0].data * OLDPFLAT ## / NEWPFLAT # making ttfile FLAT which we do not really want (testing)
    if makeitflat!=None:
        print "MAKING IT FLAT FLAT FLAT LIKE A FLAT"
        fin[0].data = fin[0].data / NEWPFLAT
    fin.close()

    print ttfile,ttfile2
#sys.exit(1)

    # zeroing bad pixels in input data
    dq0 = pyfits.open(filename)["DQ",1].data
    bit_mask = (4+16+32+128)
    dq = bitwise_and(dq0,zeros(shape(dq0),'Int16')+ bit_mask)
    print dq
    fin = pyfits.open(ttfile2,mode="update")
    fin[0].data[dq!=0] = 0.
    fin.close()
    

    sex = sextractor2.SExtractor()
    
    sex.config["DETECT_TYPE"] = "CCD"
    sex.config["DETECT_MINAREA"] = 20.
    sex.config["DETECT_THRESH"] = .5
    sex.config["ANALYSIS_THRESH"] = 1

    sex.config["FILTER"] = "Y"

    sex.config["DEBLEND_NTHRESH"] = 64
    sex.config["DEBLEND_MINCONT"] = 0.001
 
    sex.config["CLEAN"] = "Y"
    sex.config["CLEAN_PARAM"] = 1.0

    sex.config["MASK_TYPE"] = "CORRECT"

    sex.config["PHOT_APERTURES"] = 5
    sex.config["PHOT_AUTOPARAMS"] = [2.5,3.5]

    sex.config["SATUR_LEVEL"] = 500000.

    sex.config["MAG_ZEROPOINT"] = 26.27
    sex.config["MAG_GAMMA"] = 4.0
    
    sex.config["GAIN"] = 2.3
    sex.config["PIXEL_SCALE"] = 0.

    sex.config["SEEING_FWHM"] = 0.27
 
    sex.config["BACK_TYPE"] = "AUTO"
    sex.config["BACK_SIZE"] = 64
    sex.config["BACK_FILTERSIZE"] = 3
    sex.config["BACKPHOTO_TYPE"] = "GLOBAL"
    sex.config["BACKPHOTO_THICK"] = 100
   
    sex.config["CHECKIMAGE_TYPE"] = "SEGMENTATION"
    sex.config["CHECKIMAGE_NAME"] = tempfile.mktemp(".fits",dir="/tmp/",prefix="py-check")
    
    sex.config["MEMORY_OBJSTACK"] = 20000
    sex.config["MEMORY_PIXSTACK"] = 2000000
    sex.config["MEMORY_BUFSIZE"] = 1024
    sex.config["VERBOSE_TYPE"] = "NORMAL"
    sex.update_config()

    # Running Sextractor on flatter version of the FLT file
    sex.run(ttfile2)
    #catalog = sex.catalog()
    #print catalog
    
    if os.path.isfile(sex.config["CATALOG_NAME"]):
        os.unlink(sex.config["CATALOG_NAME"])
    print "Opening ",sex.config["CHECKIMAGE_NAME"]
    mask = pyfits.open(sex.config["CHECKIMAGE_NAME"])[0].data


    # Grow the mask
    mask[mask!=0] = 1.
    
    check2file = tempfile.mktemp(".fits")

    write_fits(check2file,mask*1.)

    check3file = tempfile.mktemp(".fits")
        
    iraf.imfilter()
    iraf.gauss.input = check2file
    iraf.gauss.output = check3file
    iraf.gauss.sigma = 10.
    iraf.gauss.ratio = 1.0
    iraf.gauss.nsigma = 1.0
    iraf.gauss.run(mode='h')
    
    mask2 = pyfits.open(check3file)[0].data
    ##write_fits("mask2.fits",mask2*1.)
    mask2[mask2<0.03] = 0.
    mask2[mask2>0] = 1.
    ##write_fits("mask2b.fits",mask2*1.)
    ##write_fits("dq.fits",dq0)
    
    # Add the orginal mask to the smooth one to keep small makss too
    mask = mask + mask2
    mask[mask>0] = 1.

    check4file = tempfile.mktemp(".fits")

    
    ##write_fits(check4file,mask)
    
    # data contains the un flatfielded data stored in ttfile   
    data = pyfits.open(ttfile)[0].data
    ##write_fits("d.fits",data)

    
    
    print "Removing objects..",
    data[mask!=0] = badval
    ngoods = len(ravel(data[data!=badval]))
    print ngoods, "pixels left."
    
    print "Removing DQ pixels...",
    data[dq!=0] = badval
    ngoods = len(ravel(data[data!=badval]))
    print ngoods, "pixels left."

    print "Removing persistance affected pixels...",
    data[persist_mask>0] = badval
    ngoods = len(ravel(data[data!=badval]))
    print ngoods, "pixels left."
    data11 = data*1.0
    data1 = ravel(data)
    ngoods = len(data1[data1!=badval])
    print "ngoods: ",ngoods
       
    badpixmask = (data==badval)
    
    if ngoods>100:
        data2 = data1[data1!=badval]
        print data2
        sky = median(data2)
        print "Median Sky:",sky

        try:
            mode,m1,m2 = mode_conf2(data2)
            #mode,m1,m2 = mode_conf(data2)
            #mode,m1,m2 = mode_conf2(data2)

        except ValueError:
            mode = None
            m1 = None
            m2 = None
        print "Mode Sky:",mode
        data[data!=badval] = data[data!=badval]/sky
        
        
        
    else:
        sky = None
        mode = None
        data = None
        
    if logit!=None:
        log_sky(filename,sky,mode,ngoods,start,end)
        
    ##write_fits("mask2.fits",mask*1.)
    # Keeping track of total number of e-/s
    if ALLSUM==None:
        ALLSUM = data11 * 0.0
    fin = pyfits.open(filename)
    exptime = fin[0].header["EXPTIME"]
    thisdata = fin[1].data
    vg = data!=badval
    print "shapes:",shape(vg),shape(ALLSUM),filename
    ALLSUM[vg] = ALLSUM[vg] + thisdata[vg]*exptime
    #write_fits("ALLSUM.fits",ALLSUM)


    #print "continue?"
    #sys.stdin.readlines()
    

    #write_fits("dm.fits",data)
    #sys.stdin.readlines()
    #sys.exit(1)
    
    print check2file
    print check3file
    print check4file
    #sys.exit(1)
    sex.clean()
    if os.path.isfile(ttfile):
        os.unlink(ttfile)
    if os.path.isfile(check2file):
        os.unlink(check2file)
    if os.path.isfile(check3file):
        os.unlink(check3file)
    if os.path.isfile(check4file):
        os.unlink(check4file)
        
    # data contains a un-flattened/mask image, we remove the flatfielding to get an image of the background with the FF effect
    if data!=None:
        ###data = data * OLDPFLAT
        data = data #* NEWPFLAT
        # Remask the bad pixels in data
        data[badpixmask] = badval
    
    if os.path.isfile(ttfile2):
        os.unlink(ttfile2)
    if os.path.isfile(sex.config["CHECKIMAGE_NAME"]):
        os.unlink(sex.config["CHECKIMAGE_NAME"])
        
    return data,ALLSUM #[300:600,300:600]
    
def make_fake(n,xs,ys):
    """create n fake data sets of size (xs,ys)"""
    res = []
    for i in range(n):
        data = zeros([ys,xs],float32)+i
        res.append(data)
    return res
    
def combine_row(datalist, badpix=None):
    """Compute the mean of the data in datalist, ignoring any values equal to badpix
    Row by row"""
    # concatenate all data
    ys,xs = shape(pyfits.open(datalist[0])[0].data)
    n = len(datalist)    
    data = zeros([n,ys,xs],float32)

    newdata = zeros([ys,xs],float32)
    newdata2 = zeros([ys,xs],float32)
    newdata3 = zeros([ys,xs])
    for i in range(xs):
        print "Row ",i," of ",xs
        rows = []
        for ii in range(n):
            #print "Reading ",datalist[ii]
            row = pyfits.open(datalist[ii])[0].data[i]
            rows.append(row)
        rows = asarray(rows)
        print shape(rows)
        rowst = transpose(rows)
        for j in range(ys):
            col = rowst[j]
            col1 = col[col!=badpix]
            avg1,std1 = Stats2.avg(col1,nsig=3,stdev=1.)

            if len(col1)>1:
                histo = histogram(col1,100)
                vals = histo[0]
                lbins = histo[1][:-1] # left edges of bins # (histo[1][:-1]+histo[1][1:])/2.
                rbins = histo[1][1:]
                mode1 =  (lbins[vals==max(vals)][0] + rbins[vals==max(vals)][0])/2.


            if len(col1)<1:
                newdata[i][j] = 1.0
            else:
              #  newdata[i][j] = avg1
                newdata[i][j] = mode1

                newdata2[i][j] = std1
                newdata3[i][j] = len(col1)
    return newdata,newdata2,newdata3
 
class combine_pixel_func(Thread):
    def __init__(self,files,row,badpix=32):
        Thread.__init__(self)
        self.files = files
        self.row = row
        self.n = len(self.files)
        self.modes = []
        self.badpix = badpix
    def run(self):
        try:
            ys,xs = shape(self.files[0][0].data)
        except:
            print self.files
            print self.files[0]
            print self.files[0][0]
            print self.files[0][0].data

        for j in range(ys):
            pixels = []
            for ii in range(self.n):
                try:
                    pixel = self.files[ii][0].data[self.row][j]
                    pixels.append(pixel)
                except:
                    print ii
                    print self.files[ii]
                    print self.files[ii][0]
                    print self.files[ii][0].data
                    print self.files[ii][0].data[self.row]
                    print self.files[ii][0].data[self.row][j]

            pixels = asarray(pixels)
            col1 = pixels[pixels!=self.badpix]

            if len(col1)>2:
                col1 = pixels[pixels!=self.badpix]
                #histo = histogram(col1,10)
                #vals = histo[0]
                #lbins = histo[1][:-1] # left edges of bins # (histo[1][:-1]+histo[1][1:])/2.
                #rbins = histo[1][1:]
                #mode1 =  (lbins[vals==max(vals)][0] + rbins[vals==max(vals)][0])/2.
                mode1 = 3*median(col1) - 2*average(col1)
            else:
                mode1 = 1.
            self.modes.append(mode1)

def combine_pixel_threaded(datalist, badpix=None):
    """Compute the mean of the data in datalist, ignoring any values equal to badpix
    Row by row"""
    import sys
    
    # concatenate all data
    ys,xs = shape(pyfits.open(datalist[0])[0].data)
    n = len(datalist)    
    data = zeros([n,ys,xs],float32)

    newdata = zeros([ys,xs],float32)
    newdata2 = zeros([ys,xs],float32)
    newdata3 = zeros([ys,xs])
    
    files = []
    for ii in range(n):
        print "Reading ",datalist[ii]
        fits = pyfits.open(datalist[ii])
        files.append(fits)    
        
    sys.stdin.readlines()
    thrs = []
    for i in range(10):
        print "Row ",i," of ",xs
        thr = combine_pixel_func(files,i)
        thr.start()
        thrs.append(thr)
     
    running = oldrunning = 2
    while running>1:
        running = threading.activeCount()
        if running!=oldrunning:
            print running," threads running"
        oldrunning = running
        
    for thr in thrs:
        newdata[thr.row] = asarray(thr.modes)
        
    return newdata,newdata2,newdata3


def combine_pixel(datalist, badpix=None,makeitflat=0,cubename="cube.dat"):
    """Compute the median of the data in datalist, ignoring any values equal to badpix
    Row by row"""
    import copy
    import numpy.ma as ma


    # concatenate all data
    ys,xs = shape(pyfits.open(datalist[0][0])[0].data)
    n = len(datalist)    
    data = zeros([n,ys,xs],float32)

    ##modedata = zeros([ys,xs],float32)
    avgdata = zeros([ys,xs],float32)
    ndata = zeros([ys,xs])
    stdevdata = zeros([ys,xs],float32)
    meddata = zeros([ys,xs],float32)
    
    alldata = []
    files = []
    print "Opening ",n,"files"
    for ii in range(n):
        print "Reading ",ii,"of",n,datalist[ii]
        #fits = pyfits.open(datalist[ii])
        #files.append(fits)    
        d = copy.copy(pyfits.open(datalist[ii][0])[0].data)
        d2 = copy.copy(pyfits.open(datalist[ii][0])[0].data)
        print "adding ",datalist[ii][0],datalist[ii][1],datalist[ii][2],datalist[ii][1]*datalist[ii][2]
        print "Loading new ff:",datalist[ii][3]
        
        FILTER = pyfits.open(datalist[ii][3])[0].header["FILTER"]
        if not NEWPFLATS.has_key(FILTER):
            print "Filter ",FILTER," is unkown. Moving on!!"
            return None,None
        
        NEWPFILE = NEWPFLATS[FILTER]
        NEWPFILE = os.path.join("/grp/hst/cdbs/iref/",NEWPFILE)
        NEWPFLAT = pyfits.open(NEWPFILE)[1].data[5:1014+5,5:1014+5]

        d2[d2==badpix] = 0.
        if ii==0:
            allsum = d2*datalist[ii][1]*datalist[ii][2]
        else:
            allsum = allsum + d2*datalist[ii][1]*datalist[ii][2]
        print average(allsum)
        
###        # The input is now a normalized, un flatfielded sky image.. we only need to combine these to get an estimate of the flat-field/
        if makeitflat!=0:
            print "Applying new FF",
            d[d!=badpix] = d[d!=badpix] / NEWPFLAT[d!=badpix]
            print ".Done"
            #sys.exit(1)
        alldata.append(d)
          
    # Save the data cube...
    import pickle
    print "Saving Data Cube in ",cubename
    pickle.dump(alldata,open(cubename,"w"),pickle.HIGHEST_PROTOCOL)
    print "Done"

    for i in range(xs):
    #for i in range(300,310):
        #if i!=185: continue

        print "Row ",i," of ",xs
        rows = []
        print "Reading rows..",
        for ii in range(n):
            #rows.append(copy.copy(pyfits.open(datalist[ii])[0].data[i]))
            rows.append(alldata[ii][i])
            
        rows = asarray(rows)
        #rows[rows==badpix] = NaN 
        mx = ma.masked_values(rows,badpix)
        print shape(rows),shape(mx)
        print "Computing row median:",
        med0 = ma.median(mx,axis=0)
        std0 = ma.std(mx,axis=0)
        avg0 = ma.mean(mx,axis=0)
        print "===>",med0,shape(med0),shape(std0),shape(mx)
        print "Done."

        show = 0
        if min(med0)<0.1 and show==1:
            print min(med0)
            idx = nonzero(med0==min(med0))[0][0]
            print idx
            print shape(transpose(mx)[idx])
            xxxx  = transpose(mx)[idx]
            import pylab

            pylab.plot(range(len(xxxx)),xxxx)
            pylab.show()
            import np as norp
            pylab.hist(xxxx,norp.get_n(xxxx))
            pylab.show()
            raw_input("med<0")
            xxxx  = transpose(mx)[idx-2]
            import pylab

            pylab.plot(range(len(xxxx)),xxxx)
            pylab.show()
            import np as norp
            pylab.hist(xxxx,norp.get_n(xxxx))
            pylab.show()


        meddata[i] = med0
        stdevdata[i] = std0
        avgdata[i] = avg0
        print "Done."
        
        print "Computing n:",
        nd = sum(mx/mx,axis=0)
        print "Done"
        ndata[i] = nd
        
        s = []
        print shape(rows)
        r2= transpose(rows)
        if i==400:
            for j in range(len(r2[i])):
                print rows[j]
                ss = "%d %f\n" % (j,float(r2[i][j]))
                s.append(ss)
            open("col.txt","w").writelines(s)
        continue
        
        
        for j in range(ys):
            pixels = []
            for ii in range(n):
            #print "Reading ",datalist[ii]
                pixel = rows[ii][j]
#                pixel = files[ii][0].data[i][j]
                pixels.append(pixel)
            pixels = asarray(pixels)


            col1 = pixels[pixels!=badpix]
            #print "len:",len(pixels),len(col1),pixels
            avg1 = 0.
            std1 = 0.
            avg1,std1 = Stats2.avg(col1,nsig=3,stdev=1.)
            print "======>",avg1,std1

            if avg<0:
                raw_input("...")
            if len(col1)<=1:
                mode1 = 1.

            if len(col1)>1:
                #histo = histogram(col1,10)
                #vals = histo[0]
                #lbins = histo[1][:-1] # left edges of bins # (histo[1][:-1]+histo[1][1:])/2.
                #rbins = histo[1][1:]
                ##mode1a =  (lbins[vals==max(vals)][0] + rbins[vals==max(vals)][0])/2.
                med1 = median(col1)
                ##mode1 = 3*med1 - 2*sum(col1)/len(col1)
                #print "data:",avg1,med1,mode1a,mode1

            if len(col1)<=1:
                ##modedata[i][j] = 1.0
                avgdata[i][j] = 1.0
                ndata[i][j] = 0
                stdevdata[i][j] = 0.0
                meddata[i][j] = 0.0
            else:
                ##modedata[i][j] = mode1
                avgdata[i][j] = avg1
                ndata[i][j] = len(col1)
                stdevdata[i][j] = std1
                meddata[i][j] = med1
            #if abs(i-300)<5 and abs(j-300)<5:
            #    save_pixel("%d_%d.pixel" % (i,j),col1)
    #for f in files:
    #    f.close()
    ##return modedata,avgdata,ndata,stdevdata,meddata
    return ndata,avgdata,stdevdata,meddata,allsum

def save_pixel(name,data):
    s = []
    for i in range(len(data)):
        ss = "%d %f\n" % (i,data[i])
        s.append(ss)
    open(name,"w").writelines(s)
    
    
def combine(datalist,badpix=None):
    """Compute the mean of the data in datalist, ignoring any values equal to badpix"""
    # concatenate all data
    ys,xs = shape(datalist[0])
    n = len(datalist)    
    data = zeros([n,ys,xs],float32)

    for i in range(n):
        data[i] = datalist[i]
    datat = transpose(data)
    
    newdata = zeros([ys,xs],float32)
    newdata2 = zeros([ys,xs],float32)
    newdata3 = zeros([ys,xs])

    for i in range(xs):
        print "Row ",i," of ",xs
        for j in range(ys):
            #print i
            col = datat[i][j]
            col1 = col[col!=badpix]
            #avg1 = sum(col1)/len(col1)
            avg1,std1 = Stats2.avg(col1,nsig=3,stdev=1)
            if len(col1)<1:
                newdata[j][i] = 1.0
            else:
                newdata[j][i] = avg1
                newdata2[j][i] = std1
                newdata3[j][i] = len(col1)
    #newdata[newdata == nan] = 1.0
    
    return newdata,newdata2,newdata3

def get_data(filelist):
    """Load a set of FITS file into a list of array. Loads extension sci,1"""
    res = []
    for f in filelist:
        data = pyfits.open(f)[1].data
        res.append(data)
    return res
        
def write_fits(filename,data):
    import os
    phdu = pyfits.PrimaryHDU(data=asarray(data,float32))
    if os.path.isfile(filename):
        os.unlink(filename)
    phdu.writeto(filename)


if __name__=="__main__":
    import sys,glob
    
    files = glob.glob(sys.argv[1])
    print files
    res = []
    for f in files:
        data = run_sextractor(f,badval = 32)
        print shape(data)
        sys.exit(1)
        res.append(data)
#    res = make_fake(10,100,100)

    print "Combining...",
    newdata,newdata2,newdata3 = combine(res,badpix=32)
    print "Done."
    
    write_fits("test.fits",newdata)
    write_fits("test2.fits",newdata2)
    write_fits("test3.fits",newdata3)
    
