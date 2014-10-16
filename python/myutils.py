""" utility functions """


import sys
import numpy as np

class ProgressMonitor:
	"""monitor progress of your loops"""
	def __init__(self,total_steps, total_show=100, init_show=0):
		"""init_show: initial value progress percentage, set to 0 if no reason
		total_steps: maximum iteration steps to monitor
		total_show: number of revealing times
		"""
		self.current_show=int(init_show)
		self.total_steps=total_steps
		self.total_show=total_show
		print " %02d%%"%(self.current_show*100/self.total_show),
		sys.stdout.flush()
		
	def monitor_progress(self,current_step):
		"""put this inside loop to monitor progress, to print the percent of
		job finished."""
		#print when current_step first exceeds current show percent
		if current_step>=self.total_steps*self.current_show/self.total_show:
			print "\b\b\b\b\b %02d%%"%(self.current_show*100/self.total_show),
			sys.stdout.flush()
			self.current_show+=1
			

def skeleton(x,y,nbin=10,alpha=0.683):
	"""function [xmed,ymed,ylim,xm,ym,ysig,count]=skeleton(x,y,nbin,alpha)
	% to divide x into bins and give estimation of center and variance of y
	% inside each bin
	%input:
	% x,y: column vectors to extract skeleton from
	% nbin: number of bins or bin edges for x
	% alpha: confidence level for boundary estimation
	%"""

	x=np.array(x)
	y=np.array(y)
	
	[count,xbin]=np.histogram(x,nbin)
	nbin=len(xbin)-1
	bin=np.digitize(x,xbin)-1
	
	xm=np.empty(nbin)
	#ym=xm[:]     #this is wrong! even though id(ym)!=id(xm), and id(ym[0])!=id(xm[0])
	ym=np.empty_like(xm)
	ysig=np.empty_like(xm)
	xmed=np.empty_like(xm)
	ymed=np.empty_like(xm)
	ylim=np.empty([2,nbin]);
	alpha=(1-alpha)/2;
	
	for i in xrange(nbin):
		xm[i]=np.mean(x[bin==i])
		xmed[i]=np.median(x[bin==i])
		ym[i]=np.mean(y[bin==i])
		ymed[i]=np.median(y[bin==i])
		ysig[i]=np.std(y[bin==i])
		tmp=np.sort(y[bin==i])
		if count[i]:
			ylim[:,i]=[tmp[np.ceil(alpha*count[i])],tmp[np.ceil((1-alpha)*count[i])]]
		else:
			ylim[:,i]=[NaN,NaN]
			
	return {'x':{'median':xmed,'mean':xm,'bin':xbin,'hist':count},'y':{'median':ymed,'mean':ym,'std':ysig,'CI':ylim}}
