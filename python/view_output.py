import h5py
from pylab import *
datadir='/work/Projects/Lensing/outputv4/data/'
#datadir='/mnt/Bright/'
outdir='/work/Projects/Lensing/outputv4/'
##
ids=range(5)
f=list(ids)
for i in ids:
    f[i]=h5py.File(datadir+'WL_Z%d.hdf5'%(i+1),'r')
##    
figure()
for i in ids:
    plot(f[i]['/shear/seperation'][:],f[i]['/rand/profile'][:],'-o')
    hold('on')
legend(['Z%d'%(i+1) for i in ids])   
xscale('log')
yscale('symlog',linthres=1)
##    
figure()
for i in ids:
    errorbar(f[i]['/shear/seperation'][:],f[i]['/rand/profile'][:]/f[i]['/predict0/profile'][:],f[i]['/rand/profile_err'][:]/f[i]
['/predict0/profile'][:],fmt='o--',label='Z%d'%(i+1))

legend(numpoints=1)    
xlabel('Comoving Seperation[Mpc/h]')
ylabel('Relative background shear')
xscale('log')
yscale('symlog',linthres=1)

#savefig(outdir+'SysShear_RedShift.eps')
