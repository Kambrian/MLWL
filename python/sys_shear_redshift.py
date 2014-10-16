# -*- coding: utf-8 -*-

import h5py
import triangule 
import myutils
from pylab import *
import pp
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

chi0=sqrt(e0[0,]**2+e0[1,]**2)
theta0=arctan2(e0[1,],e0[0,])/2
theta0=theta0*180/pi

raminmax=array([[129,141],[174,186],[211.5,223.5]])
decminmax=array([[-1,3],[-2,2],[-2,2]])


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
	
rad=10
zbin=arange(0,1,0.02)
append(zbin,10)
nbin=len(zbin)-1
shear_sum=zeros(nbin)
average_shear=zeros(nbin)

ncpus=min(nbin,24)
ppservers=()
job_server = pp.Server(ncpus,ppservers=ppservers)
print "Starting pp with", job_server.get_ncpus(), "workers"
outputs=list(range(nbin))

for i in xrange(nbin):
    filt= logical_and(z0>zbin[i],z0<zbin[i+1])
    outputs[i]=job_server.submit(average_shear_at_r,args=(rad,e0[:,filt],x0[:,filt],raminmax[0],decminmax[0]),modules=('triangule',)) 

for i in xrange(nbin):
    shear_sum[i],angle_sum=outputs[i]()
    if angle_sum==0:
        average_shear[i]=nan
    else:
        average_shear[i]=shear_sum[i]/angle_sum
    print "[%.1f,%.1f]:"%(zbin[i],zbin[i+1]),shear_sum[i],average_shear[i]

job_server.print_stats()

figure()    
plot(zbin[:-1],shear_sum,'o-')
xlabel('Reshift')
ylabel('Systematic Shear Contribution')
title('First GAMA Region')
savefig('sys_shear_redshift_contrib.eps')
figure()
plot(zbin[:-1],average_shear,'o-')
xlabel('Reshift')
ylabel('Systematic Shear')
title('First GAMA Region')
savefig('sys_shear_redshift_average.eps')
