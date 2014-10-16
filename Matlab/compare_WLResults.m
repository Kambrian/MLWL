cd /mnt/Bright/
file='WL_D3.hdf5';
r=hdf5read(file,'/shear/seperation');
s=hdf5read(file,'/shear/profile');
s0=hdf5read(file,'/predict0/profile');
s1=hdf5read(file,'/predict1/profile');
es=hdf5read(file,'/shear/profile_err');
z=hdf5read(file,'/predict0/z');
M0=hdf5read(file,'/predict0/Mmean');
M1=hdf5read(file,'/predict1/Mmean');
%%
figure;
ploterr(r,s,[],es,'logxy')
hold on;
plot(r,s0,'r',r,s1,'g')
%%
[cf0,gof0,out0]=createFit(r,s0,es.^-2,z);
[cf1,gof1,out1]=createFit(r,s1,es.^-2,z);
[cf,gof,out]=createFit(r,s,es.^-2,z);
%%
par=coeffvalues(cf);
plot(r,lensing_rfunc(r,par(1),0,z,par(3)),'--');
par0=coeffvalues(cf0);
plot(r,lensing_rfunc(r,par0(1),0,z,par0(3)),'r--');
par1=coeffvalues(cf1);
plot(r,lensing_rfunc(r,par1(1),0,z,par1(3)),'g--');