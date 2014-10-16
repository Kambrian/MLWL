import h5py
import matplotlib
#matplotlib.use('agg')
from pylab import *
ion()
rc('text', usetex=True)

indir='/mnt/charon/Lensing/output/'
outdir='/work/Projects/Lensing/outputv4/paper/'
files=['WL_L1.DZ3.hdf5','WL_L2.DZ3.hdf5','WL_L3.DZ3.hdf5','WL_L3.DZ0.hdf5']

fmts=['ro','gs','bd','kd-']
for i in range(len(files)):
    f=h5py.File(indir+files[i],'r')
    n=f['/shear/numpair'][:]
    r=f['/shear/seperation'][:]
    nr=f['/rand/numpair'][:]
    nerr=f['/rand/numpair_err'][:]
    errorbar(r,n/nr-1,nerr/nr*2,fmt=fmts[i])
plot([1e-2,10],[0,0])    
legend(('L1','L2','L3',r'L3,$\Delta z$=0.01'))
xscale('log')
xlabel('comoving seperation [Mpc/h]')
ylabel(r'$n/n_{rand}-1$')
savefig(outdir+'Contamination.eps')