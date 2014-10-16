# -*- coding: utf-8 -*-

import h5py
import triangule as t
from myutils import ProgressMonitor
from pylab import *
#srcdir='/work/Projects/Lensing/data/'
srcdir='/data/raid3/jxhan/Lensing/data/'
f=h5py.File(srcdir+'shearmap_ZEBRA.mat','r')
data=list(range(3))
data[0]=f['/e1']
data[1]=f['/e2']
data[2]=f['/e3']
x0=data[0][0:2,:]
e0=data[0][2:4,:]
e0[0,:]*=-1
z0=data[0][5,:]
x1=data[1][0:2,:]
e1=data[1][2:4,:]
e1[0,:]*=-1
z1=data[1][5,:]
x2=data[2][0:2,:]
e2=data[2][2:4,:]
e2[0,:]*=-1
z2=data[2][5,:]

chi0=sqrt(e0[0,]**2+e0[1,]**2)
theta0=arctan2(e0[1,],e0[0,])/2
theta0=theta0*180/pi

raminmax=array([[129,141],[174,186],[211.5,223.5]])
decminmax=array([[-1,3],[-2,2],[-2,2]])

rad=logspace(-1,1,20)
average_shear=zeros(len(rad))
for i in xrange(len(rad)):
    a=0.
    b=0.
    print "%d/%d"%(i,len(rad)),rad[i],": ", 
    progress=ProgressMonitor(x0.shape[1])
    for galid in xrange(x0.shape[1]):
        progress.monitor_progress(galid)
        arcs=t.get_inner_arcs(x0[:,galid],rad[i],raminmax[0],decminmax[0])     
        shear_sum,angle_sum=t.integrate_shear_cart_in_arcs(e0[:,galid],arcs)
        a+=shear_sum
        b+=angle_sum
    if b==0:
        average_shear[i]=nan
    else:
        average_shear[i]=a/b*180/pi
    print average_shear[i]
    
plot(rad,average_shear)
print rad
print average_shear  
xscale('log')


#datadir='/work/Projects/Lensing/outputv4/data/'
datadir='/data/raid3/jxhan/Lensing/output/'
f=h5py.File(datadir+'randprof_Gal10_1234.hdf5','r')

plot(f['/shear/seperation'][:]*180/pi,f['/shear/profile'][:]*2,'o')
savefig('sys_shear_py.eps')
