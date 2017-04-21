# 
from pylab import *
from matplotlib import *
import numpy as np
#************************************#
def GetCramer(nbin=100):
    T10,W10=getW(dirname='T=10');
    P10,w10=histogram(W10,normed=1,bins=nbin);
    T100,W100=getW(dirname='T=100');
    P100,w100=histogram(W100,normed=1,bins=nbin);
    T400,W400=getW(dirname='T=400');
    P400,w400=histogram(W400,normed=1,bins=nbin);
    C10=-log(P10)/10.;
    C100=-log(P100)/100.;
    C400=-log(P400)/400.;
    clf()
    plot(w10[1:],C10-min(C10),'sr')
    plot(w100[1:],C100-min(C100),'^b')
    plot(w400[1:],C400-min(C400),'ok')
    grid()
    return C10,C100,C400,w10,w100,w400
#************************************#
def Wall(nbin=100):
    dirname=['T=1','T=10','T=20','T=30','T=100','T=400'];
    s=shape(dirname);
    ndir=s[0];
    T1,W1=getW(dirname=dirname[0]);
    P1,w1=histogram(W1,normed=1,bins=nbin);
    T10,W10=getW(dirname=dirname[1]);
    P10,w10=histogram(W10,normed=1,bins=nbin);
    T20,W20=getW(dirname=dirname[2]);
    P20,w20=histogram(W20,normed=1,bins=nbin);
    T30,W30=getW(dirname=dirname[3]);
    P30,w30=histogram(W30,normed=1,bins=nbin);
    T100,W100=getW(dirname=dirname[4]);
    P100,w100=histogram(W100,normed=1,bins=nbin);
    T400,W400=getW(dirname=dirname[5]);
    P400,w400=histogram(W400,normed=1,bins=nbin);
    #plot(w1[1:],P1,'r.')
    return P1,P10,P20,P30,P100,P400,w1,w10,w20,w30,w100,w400
#************************************#
def getW(fname='tseries.out',dirname='T=10'):
    ts=np.loadtxt(dirname+'/'+fname);
    s=shape(ts);
    nt=s[0]-1;
    Np=s[1]-1;
    W=ts[2:,1:];
    T=ts[2,0]-ts[1,0];
    W=W/T;
    return T,W
#************************************#
def pW(fname='tseries.out',nbin=100):
    T,W=getW(fname);
    P,w=histogram(W,normed=1,bins=nbin);
    clf()
    plot(w[1:],log10(P),'.-b')
    grid()
#************************************#
def rts(fname='tseries.out',ip=1):
    ts=np.loadtxt(fname);
    time=ts[:,0];
    x=ts[:,2*ip-1];
    v=ts[:,2*ip];
    return time,x,v
#************************************#
def trev(ip=1,lplot=0):
    t,x1,v1=rts(ip=ip);
    nt=t.size;
    ntby2=nt/2;
    xfor=x1[0:ntby2];
    vfor=v1[0:ntby2];
    xback=zeros(ntby2);
    vback=zeros(ntby2)
    for it in range(0,ntby2):
        xback[it]=x1[nt-it-1];
        vback[it]=-v1[nt-it-1];
    if (lplot):
        clf()
        plot(xfor,vfor,'r.-')
        plot(xfor[0],vfor[0],'ro')
        plot(xfor[ntby2],vfor[ntby2],'rs')
        plot(xback,vback,'.-k')
        plot(xback[0],vback[0],'ok')
        plot(xback[-1],vback[-1],'sk')
    return xfor,vfor,xback,vback
#************************************#
def plot_traj(np=4):
    clf()
    for ip in range(1,np):
        xfor,vfor,xback,vback=trev(ip)
        plot(xfor,vfor,'r.-')
        plot(xfor[0],vfor[0],'ro')
        plot(xfor[-1],vfor[-1],'rs')
        plot(xback,vback,'.-k')
        plot(xback[0],vback[0],'ok')
        plot(xback[-1],vback[-1],'sk')
#************************************#
if __name__ == '__main__':
  main();
else:
  print 'Analysis routine for PiF'
