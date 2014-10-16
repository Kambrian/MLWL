# -*- coding: utf-8 -*-

import h5py
import triangule 
import myutils
from pylab import *
import pp
srcdir='/work/Projects/Lensing/data/'
#srcdir='/data/raid3/jxhan/Lensing/data/'
f=h5py.File(srcdir+'shearmap_ZEBRA.mat','r')
#f=h5py.File(srcdir+'cfhtsrc.hdf5','r')
data=list(range(3))
data[0]=f['/e1']
#~ data[1]=f['/e2']
#~ data[2]=f['/e3']
ngal=data[0].shape[1]
#ngal=100
x0=data[0][0:2,:ngal]
e0=data[0][2:4,:ngal]
e0[0,:]*=-1
z0=data[0][5,:ngal]
#redshift filter:
zsrc=0.1
x0=x0[:,z0>zsrc+0.1]
e0=e0[:,z0>zsrc+0.1]

chi0=sqrt(e0[0,]**2+e0[1,]**2)
theta0=arctan2(e0[1,],e0[0,])/2
theta0=theta0*180/pi

raminmax=array([[129,141],[174,186],[211.5,223.5]])
decminmax=array([[-1,3],[-2,2],[-2,2]])
#raminmax=array([[132,137],[174.,186],[211.5,223.5]])
#decminmax=array([[-1,-0.94],[-2,2],[-2,2]])

ncpus=8
ppservers=()
job_server = pp.Server(ncpus,ppservers=ppservers)
print "Starting pp with", job_server.get_ncpus(), "workers"

def average_shear_at_r(rad,e,x,raminmax,decminmax):
	shear_sum=0.
	angle_sum=0.
#	progress=myutils.ProgressMonitor(x.shape[1])
	for galid in xrange(x.shape[1]):
#		progress.monitor_progress(galid)
		arcs=triangule.get_inner_arcs(x[:,galid],rad,raminmax,decminmax)     
		a,b=triangule.integrate_shear_cart_in_arcs(e[:,galid],arcs)
		shear_sum+=a
		angle_sum+=b	
	return (shear_sum/2,angle_sum*triangule.np.pi/180)
	
rad=logspace(-1,1,60)
#rad=arange(3.0,12.0,1.0)
average_shear=zeros(len(rad))
outputs=list(range(len(rad)))
for i in xrange(len(rad)):
    print "%d/%d"%(i,len(rad)),rad[i],": ", 
    outputs[i]=job_server.submit(average_shear_at_r,args=(rad[i],e0,x0,raminmax[0],decminmax[0]),modules=('triangule',)) 

for i in xrange(len(rad)):
    shear_sum,angle_sum=outputs[i]()
    if angle_sum==0:
        average_shear[i]=nan
    else:
        average_shear[i]=shear_sum/angle_sum
    print average_shear[i]

job_server.print_stats()
    
plot(rad,average_shear)
print rad
print average_shear  
xscale('log')

datadir='/work/Projects/Lensing/outputv4/data/'
#datadir='/data/raid3/jxhan/Lensing/output/'
f=h5py.File(datadir+'randprof_SKY0_1234.hdf5','r')
#f=h5py.File(datadir+'randprof_CFHT.CFHTSKY_1234.hdf5','r')
#savetxt(datadir+'sys_shear_CFHT.CFHTSKY.dat',(rad,average_shear))

plot(f['/shear/seperation'][:]*180/pi,f['/shear/profile'][:],'o')
xlabel('R/deg')
ylabel('Systematic Shear')
title('First GAMA Region')
legend(('Analytical','Random-Lens'))
savefig('sys_shear_py_SKY0.eps')
