''' to load WL outputs and apply corrections '''
import h5py
import math
import numpy as np
import matplotlib.pyplot as plt

def plot_shearfile(infile, readyplot=1):
	''' to apply systematic shear and boost corrections '''
	f = h5py.File(infile, 'r')
	r=f['/shear/seperation'][:]
	p=f['/shear/responsivity'][:]
	s=f['/shear/profile'][:]
	
	if readyplot==1:
		cov=f['/shear/covariance'][:]
		es=np.sqrt(cov.diagonal())
		plt.errorbar(r,s,es,fmt='bx')
		print es
	else:
		es=f['/shear/profile_err'][:]
		s1=f['/predict0/profile'][:]
		s2=f['/predict1/profile'][:]
		print es/p
		plt.errorbar(r,s/p,es/p,fmt='ko-')
		plt.hold('on')
		plt.plot(r,s1,'r')
		plt.plot(r,s2,'g')
	
	plt.xscale('log')
	plt.yscale('symlog')
	
	f.close()


def weight_average_correction(w,ew,n):
	'''correct for averaging-bias due 
	to the exclusion of zero count realizations
	in outputing the average weight'''
	ew=(ew**2+w**2)*n/max(n)  # <w^2>
	w=w*n/max(n)
	ew=np.sqrt(ew-w**2)
	return (w,ew,max(n))

def reduce_shear(file):
	''' to apply systematic shear and boost corrections '''
	f = h5py.File(file, 'r')
	r=f['/shear/seperation'][:]
	p=f['/shear/responsivity'][:]
	prand=f['/rand/responsivity'][:]
	w=f['/shear/weight'][:]
	cov=f['/shear/covariance'][:]
	
	(wrand,ewrand,nrand)=weight_average_correction(f['/rand/weight'][:], \
    f['/rand/weight_err'][:],f['/rand/numrand'][:])
	
	#systematic-shear correction, as well as responsivity correction
	s=f['/shear/profile'][:]/p-f['/rand/profile'][:]/prand
	es=np.sqrt(cov.diagonal())/p
	#add photo-z err
	es2=es**2+(f['/sphoto/profile_err'][:]/f['/sphoto/responsivity'][:])**2
	#add background(systematic-shear) error
	es2=es2+(f['/rand/profile_err'][:]/prand/math.sqrt(nrand))**2
	
	#boost correction
	es2=es2/s**2+(ewrand/wrand)**2  #relative error
	s=s*w/wrand
	es=np.sqrt(es2)*s
	
	coveff=np.zeros(cov.shape)
	for i in range(0,cov.shape[0]):
		for j in range(0,cov.shape[1]):
			coveff[i,j]=cov[i,j]/math.sqrt(cov[i,i])/math.sqrt(cov[j,j])
	

	f.close()

	return (r,s,es,coveff)
	
def plot_shear(runname):
	''' plot shear and predictions from file '''
	#file='/mnt/Bright/WL_'+runname+'.hdf5'
	file='/work/Projects/Lensing/outputv4/data/WL_'+runname+'.hdf5'
	print "loading "+file+"..."
	(r,s,es,c)=reduce_shear(file)
	#h=plt.figure()
	l1=plt.errorbar(r,s,es,fmt='o')
	l1[0].set_label('WL')
	plt.hold('on')
	plt.xscale('log')
	#plt.yscale('symlog',linthreshy=10)
	#plt.ylim(0,1e5)
	plt.yscale('log',nonposy='clip') #this would not work for old matplotlib
	plt.ylim(1,1e5)
	
	f=h5py.File(file,'r')
	l2,=plt.plot(r,f['/predict0/profile'],'r',label='Dyn')
	l3,=plt.plot(r,f['/predict1/profile'],'g',label='Lum')
	l4,=plt.plot(r,f['/predict3/profile'],'r--',label='DynNew')
	l5,=plt.plot(r,f['/predict4/profile'],'g--',label='LumNew')
	#h.legend((l1[0],l2,l3),('WL','Dyn','Lum'),'upper right')
	plt.title(runname)
	#plt.show()

	#return (l1,l2,l3) 
	return (l1,l2,l3,l4,l5) 

def plot_runs(runs):
	''' plot many runs together 
	runs should be in the form of [(a,b,c),(d,'',f),...]
	where the each inner tuple will be plotted in a row.'''
	nrows=len(runs)
	ncols=len(runs[0])
	fig,ax=plt.subplots(ncols=ncols,nrows=nrows,sharey=True)
	ax=np.array(ax,ndmin=2)
	#fig.subplots_adjust(hspace=0)
	lines=[]
	for i in range(nrows):
		lrow=[]
		for j in range(ncols):
			if runs[i][j]!='':
				plt.axes(ax[i][j])
				l=plot_shear(runs[i][j])
			else:
				l=()
			lrow.append(l)
		lines.append(lrow)
	
	lshow=lines[0][-1]#legend on top right
	#lshow=(lshow[0][0],lshow[1],lshow[2])
	#ax[0][-1].legend(lshow,('WL','Dyn','Lum'),'lower left') 
	lshow=(lshow[0][0],lshow[1],lshow[2],lshow[3],lshow[4])
	ax[0][-1].legend(lshow,('WL','Dyn','Lum','DynNew','LumNew'),'lower left') 

	plt.setp([ax[i][j].get_yticklabels() for i in range(nrows) for j in range(1,ncols-1)], visible=False)
	
	if ncols>1:
		for i in range(nrows):
			ax[i][-1].yaxis.set_ticks_position('right')
			ax[i][-1].yaxis.set_ticks_position('both')
	
	#~if nrows>1:
		#~for j in range(ncols):
			#~ax[0][j].xaxis.set_ticks_position('top')
			#~ax[0][j].xaxis.set_ticks_position('both')	
	return (fig,ax,lines)

if __name__=="__main__":
	outdir='/work/Projects/Lensing/outputv4/PredictProf/'
	plot_runs([('D1','D2','D3'),('L1','L2','L3')])
	plt.savefig(outdir+'PredictProf_MassProxy.eps')
	plot_runs([('F1','F2','F3'),('F1N2','F2N2','F3N2')])
	plt.savefig(outdir+'PredictProf_Flux.eps')
	plot_runs([('Z1','Z2','Z3','Z4','Z5'),('N1','N2','N3','N4','N5')])
	plt.savefig(outdir+'PredictProf_Redshift.eps')
	plot_runs([('L1.CFHT','L2.CFHT','L3.CFHT'),('L1','L2','L3')])
	plt.savefig(outdir+'../CFHT/PredictProf_CFHT.eps')
	
	plt.show()


#rootdir='/mnt/charon/Lensing/'
#file1='randprof_Test.test_1024.hdf5'
#file1=rootdir+'code/v5.0/C/sheartest.hdf'
#file1=rootdir+'output/WL_Test.test.hdf5' #default
#file2=rootdir+'output/WL_Test.test2.hdf5' #rvircut=5
#file3=rootdir+'output/WL_Test.test3.hdf5' #thetaobscure=0
#file4=rootdir+'output/WL_Test.test4.hdf5' #offsetcut=0
#plot_shearfile(file1)
#plt.ion()
#plt.figure();
#plot_shearfile(file1,0)
#plt.figure();
#plot_shearfile(file2,0)
#plt.figure();
#plot_shearfile(file3,0)
#plt.figure();
#plot_shearfile(file4,0)