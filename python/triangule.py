""" triangular functions including sind, cosd, asind etc., 
and ccdistd()
"""
import numpy as np
#from math import *
#========================== auxilliary funcs====================
pi_h=np.pi/2
DEG2RAD=np.pi/180.

def sind(x):
   return np.sin(x*DEG2RAD)
def cosd(x):
   return np.cos(x*DEG2RAD)
def asind(x):
   return np.arcsin(x)/DEG2RAD
def acosd(x):
   return np.arccos(x)/DEG2RAD
            
def ccdistd(x,y):
	""" cosine of celestial distance(in degrees) between two points at x and y
	% 
	% calculate cos(d) where d is the angular seperation of two celestial points
	% with coordinates x=[ra1,dec1],y=[ra2,dec2]
	% x,y can should be [2xn] or [2x1] numpy.arrays
	%  *********** all angles in units of degrees ****************
	% to be improved: use Vincenty formula for numerical stability
	"""
	csd=sind(y[1,])*sind(x[1,])+cosd(y[1,])*cosd(x[1,])*cosd(y[0,]-x[0,]);
	np.where(csd>1,1,csd)
	np.where(csd<-1,-1,csd)
	return csd

def ccdist(x,y):
	""" cosine of celestial distance(in radian) between two points at x and y
	% 
	% calculate cos(d) where d is the angular seperation of two celestial points
	% with coordinates x=[ra1,dec1],y=[ra2,dec2]
	% x,y can should be [2xn] or [2x1] numpy.arrays
	%  *********** all angles in units of radians ****************
	% to be improved: use Vincenty formula for numerical stability
	"""
	csd=np.sin(y[1,])*np.sin(x[1,])+np.cos(y[1,])*np.cos(x[1,])*np.cos(y[0,]-x[0,]);
	np.where(csd>1,1,csd)
	np.where(csd<-1,-1,csd)
	return csd
		
def circle_intersect_ra(x, rad, ra, unit='deg'):
	""" intersecting points of a circle	with a RA line 
	x: (ra,dec) of the center of the circle, can be 2*n array
	rad: radius of the circle,scalar
	ra: ra of the RA line
	unit: 'radian' or 'deg'(default)
	return: (dec1,l), with dec1 the mid point dec of the chord,
			l=half_width of the chord, 
			so that the intersecting points are given
			as (ra, dec1+l) and (ra, dec1-l)
			l=0 if not intersecting, but dec1 still valid
	*****************LIMITATION:***************************
	not applicable if the distance between x and ra is >pi/2	
	"""
	#s={'radian':pi/2,'deg':90}
	u={'radian':1,'deg':DEG2RAD}
	x=np.array(x)
	
	a=pi_h-x[1,]*u[unit]
	beta=abs(ra-x[0,])*u[unit]
	b=np.arcsin(np.sin(beta)*np.sin(a)) #distance from x to RA line
	r=rad*u[unit]
	c=2*np.arctan(np.tan((a-b)/2)*np.sin((pi_h+beta)/2)/np.sin((pi_h-beta)/2)) #90-dec1,
	dec1=pi_h-c
	
	b=np.array(b,ndmin=1) #make sure b is transformed into an indexable array
	bb=b[b<r]
	gamma=np.arcsin(np.sin(bb)/np.sin(r)) #angle between rad and ra
	l=np.zeros(b.shape)
	l[b<r]=2*np.arctan(np.tan((r-bb)/2)*np.sin((pi_h+gamma)/2)/np.sin((pi_h-gamma)/2))
	
	if l.size==1: #beautify
		l=l[0]
		
	return (dec1/u[unit],l/u[unit])	
	
def circle_intersect_dec(x, rad, dec, unit='deg'):
	""" intersecting points of a circle	with a DEC line 
	x: (ra,dec) of the center of the circle, can be 2*n array
	rad: radius of the circle,scalar
	dec: dec of the DEC line
	unit: 'radian' or 'deg'(default)
	return: dRA, half_width of the chord, 
			so that the intersecting points are given
			as (ra0+dRA, dec) and (ra0-dRA, dec)
			dRA=0 if not intersecting
	*****************LIMITATION:***************************
	not applicable if the distance between x and ra is >pi/2	
	"""
	u={'radian':1,'deg':DEG2RAD}
	r=rad*u[unit]
	x=np.array(x)
	dec0=np.array(x[1,]*u[unit],ndmin=1)
	dec1=dec*u[unit]
	d=abs(dec0-dec1)
	dRA=np.zeros(d.shape)
	
	dRA[d<r]=np.arccos((np.cos(r)-np.sin(dec0[d<r])*np.sin(dec1))/(np.cos(dec0[d<r])*np.cos(dec1)))
	
	if dRA.size==1:
		dRA=dRA[0]
		
	return dRA/u[unit]	

def circle_inside_ra(x,rad,ra, side=+1, unit='deg'):
	""" position angle range for the inner part 
	(the part defined by side) of a circle 
	intersecting with a RA line
	
	x: (ra,dec) of the center of the circle, can be 2*n array
	rad: radius of the circle,scalar
	ra: ra of the RA line
	side: +1/-1, the side with smaller('-1') or greater('+1') RA, 
		  to be defined as the inner part; default: '+'
	unit: 'radian' or 'deg'(default)
	return: (flag,(tmin,tmax)), 
			flag: whether the circle intersects with the line, and if not,
				  whether the circle is completely inside or outside.
				  0: intersect.
				  1: inside
				  -1: outside
		    tmin,tmax: the position angle range, counting from x-axis 
				  (axis perpendicular to constant-RA line, 
				  positive with increasing RA, i.e., pointing eastward), 
				  rotating counter-clock-wise around the outward-pointing
				  z-axis.
	*****************LIMITATION:***************************
	not applicable if the distance between x and ra is >pi/2	
	or if pi/2-x[1]<rad
	"""
	u={'radian':1,'deg':DEG2RAD}
	x=np.array(x)
	
	a=np.array(pi_h-x[1,]*u[unit],ndmin=1)
	beta=np.array(abs(ra-x[0,])*u[unit],ndmin=1)
	r=rad*u[unit]
	#b=np.arcsin(np.sin(beta)*np.sin(a)) #distance from x to RA line
	sinabeta=np.sin(beta)*np.sin(a)
	sinr=np.sin(r)
	flag=sinabeta<sinr #sin(distance)<sin(r)
	alpha2=np.arcsin(sinabeta[flag]/sinr)
	alpha1=np.pi-alpha2
	dt1=np.zeros(a.shape)
	dt2=np.zeros(a.shape)
	tmp=np.sin((r+a[flag])/2)/np.sin((r-a[flag])/2)	
	dt1[flag]=np.pi-2*np.arctan(np.tan((beta[flag]-alpha1)/2)*tmp)
	dt2[flag]=np.pi-2*np.arctan(np.tan((beta[flag]-alpha2)/2)*tmp)
	sign=np.where(ra<x[0,],-1,1)
	tmin=pi_h-dt2*sign
	tmax=pi_h-dt1*sign
	if side<0:
		(tmin,tmax)=(tmax,tmin)
	
	if flag.size==1: #beautify
		flag=flag[0]
		tmin=tmin[0]
		tmax=tmax[0]
	
	flag=(not flag)*side*(-sign)
		
	return (flag,(tmin/u[unit],tmax/u[unit]))
	
def circle_inside_dec(x, rad, dec, side=+1, unit='deg'):
	""" position angle range for the inner part 
	(the part defined by side) of a circle 
	intersecting with a DEC line
	
	x: (ra,dec) of the center of the circle, can be 2*n array
	rad: radius of the circle,scalar
	dec: dec of the DEC line
	side: +1/-1, the side with smaller('-1') or greater('+1') DEC, 
		  to be defined as the inner part; default: '+'
	unit: 'radian' or 'deg'(default)
	return: (flag,(tmin,tmax)),
			flag: whether the circle intersects with the line, and if not,
				  whether the circle is completely inside or outside.
				  0: intersect.
				  1: inside
				  -1: outside
		    tmin,tmax: the position angle range, counting from x-axis 
				  (axis perpendicular to constant-RA line, 
				  positive with increasing RA, i.e., pointing eastward), 
				  rotating counter-clock-wise around the outward-pointing
				  z-axis.
	"""
	u={'radian':1,'deg':DEG2RAD}
	r=rad*u[unit]
	x=np.array(x)
	dec0=np.array(x[1,]*u[unit],ndmin=1)
	dec1=dec*u[unit]
	d=abs(dec0-dec1)
	
	flag=d<r
	dt=np.zeros(d.shape)	
	dt[flag]=np.arccos((np.sin(dec1)-np.cos(r)*np.sin(dec0[flag]))/(np.sin(r)*np.cos(dec0[flag])))
	
	tmin=pi_h-dt
	tmax=pi_h+dt
	if side<0:
		(tmin,tmax)=(tmax,tmin)
	
#	tmin=np.zeros(d.shape)
#	tmax=np.zeros(d.shape)
#	tmin[dt<=pi_h]=pi_h-dt[dt<=pi_h]
#	tmax[dt<=pi_h]=pi_h+dt[dt<=pi_h]	
#	tmin[dt>pi_h]=pi_h+dt[dt>pi_h]
#	tmax[dt>pi_h]=pi_h+2*np.pi-dt[dt>pi_h]	
	
	if dt.size==1:
		dt=dt[0]
		flag=flag[0]
		tmin=tmin[0]
		tmax=tmax[0]
	
	flag=(not flag)*side*np.sign(x[1,]-dec)

	return (flag,(tmin/u[unit],tmax/u[unit]))


Period={'radian':2*np.pi,'deg':360}

def angular_bin_and(bin0,bin1,unit='deg'):
	""" shared part of two angular bins 
	return a list of the AND-ed bins
	"""
	T=Period[unit]
	#unpack
	a,b=bin0[0],bin0[1]
	x,y=bin1[0],bin1[1]
	
	#a=a-np.floor(a/T)*T
	b=b-np.floor((b-a)/T)*T
	x=x-np.floor((x-a)/T)*T
	y=y-np.floor((y-a)/T)*T
	
	if y>x: #aligned bins
		u=min(b,y)
		if u>x:
			return [(x,u)]
		
	if y<x:
		if y>b:
			return [(a,b)]
		else:
			out=[(a,y)]
			if x<b:
				out.append((x,b))
			return out
	
	return []

def normalize_bin(bin,unit='deg'):	
	#shift the bin to start in a~[0,T] and make it aligned, i.e, b~[a,a+T]
	T=Period[unit]

	a=bin[0]-np.floor(bin[0]/T)*T #shift bin[0] to a~[0,T]
	b=bin[1]-np.floor((bin[1]-a)/T)*T #shift bin[1] to b~[a,a+T]
	return (a,b)
	
def normalize_binlist(binlist,unit='deg'):
	newlist=[]
	for i in xrange(len(binlist)):
		newlist.append(normalize_bin(binlist[i],unit))
	return newlist
	
def angular_binlist_and(binlist,bin1,unit='deg'):
	""" differ from angular_bin_and in that binlist is a list of bins (OR-ed) 
	binlist: [(a0,b0),(a1,b1),...]
	bin1: (x,y) or [x,y]
	"""
	
	out=[]
	for i in xrange(len(binlist)):
		out.extend(angular_bin_and(binlist[i],bin1,unit))
	
	return out

def angular_binlist_inner_and(binlist,unit='deg'):
	""" return the AND result of all the bins in binlist together
	"""
	
	n=len(binlist)
	if n<1:
		return []
		
	out=[binlist[0]]
	for i in xrange(1,n):
		out=angular_binlist_and(out,binlist[i],unit)
	
	return out

def get_inner_arcs(x,rad,raminmax,decminmax):
    binflag={}
    binbdry={}
    binflag['l'],binbdry['l']=circle_inside_ra(x,rad,raminmax[0],+1)
    binflag['r'],binbdry['r']=circle_inside_ra(x,rad,raminmax[1],-1)
    binflag['b'],binbdry['b']=circle_inside_dec(x,rad,decminmax[0],+1)
    binflag['t'],binbdry['t']=circle_inside_dec(x,rad,decminmax[1],-1)
    binlist=[]  
    for side in binflag.keys():
        if binflag[side]<0: #completely outside
            return []
        if binflag[side]==0: #intersect
            binlist.append((binbdry[side][0],binbdry[side][1]))
    if binlist==[]: #fully inside the region
        arcs=[(0,360)]
    else:    
        arcs=normalize_binlist(angular_binlist_inner_and(binlist))
    return arcs

# <codecell>

def integrate_shear_polar_in_arcs(shear0,theta0,arcs):
    """ angle in degrees
    Note the average shear evaluates to shear_sum/(angle_sum*pi/180) """
    if ( len(arcs)==1 and arcs[0][1]-arcs[0][0]==360 ):
        return (0.,360.)
    shear_sum=0.
    angle_sum=0.
    for i in range(len(arcs)):
        shear_sum+=sind(2*(theta0-arcs[i][1]))-sind(2*(theta0-arcs[i][0]))
        angle_sum+=arcs[i][1]-arcs[i][0]
    shear_sum*=shear0/2.0
    return (shear_sum,angle_sum)

def integrate_shear_cart_in_arcs(e,arcs):
    """ angle in degrees
    Note the average shear evaluates to shear_sum/(angle_sum*pi/180) """
    if ( len(arcs)==1 and arcs[0][1]-arcs[0][0]==360 ):
        return (0.,360)
    shear_sum=0.
    angle_sum=0.
    for i in range(len(arcs)):
        shear_sum-=sind(arcs[i][1]-arcs[i][0])*(e[1]*sind(arcs[i][1]+arcs[i][0])+e[0]*cosd(arcs[i][1]+arcs[i][0]))
        angle_sum+=arcs[i][1]-arcs[i][0]
    return (shear_sum,angle_sum)

def tangential_shear_cart_along_fai(e,fai):
        return -(e[1]*sind(2*fai)+e[0]*cosd(2*fai))

	
from matplotlib import pyplot as plt

def plot_circle(cen,r,linestyle='-'):
    t=np.arange(0,2*np.pi+0.2,0.1)
    x=cen[0]+r*np.cos(t)
    y=cen[1]+r*np.sin(t)
    plt.plot(x,y,linestyle)

def plot_brd(raminmax,decminmax,linestyle='k-'):
	x=[raminmax[0],raminmax[0],raminmax[1],raminmax[1],raminmax[0]]
	y=[decminmax[0],decminmax[1],decminmax[1],decminmax[0],decminmax[0]]
	plt.plot(x,y,linestyle)
    #~ plt.plot([raminmax[0],raminmax[0]],decminmax,linestyle)
    #~ plt.hold('on')
    #~ plt.plot([raminmax[1],raminmax[1]],decminmax,linestyle)
    #~ plt.plot(raminmax,[decminmax[0],decminmax[0]],linestyle)
    #~ plt.plot(raminmax,[decminmax[1],decminmax[1]],linestyle)
