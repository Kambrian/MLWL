LumPred=h5read('data/WL_L1.hdf5','/predict1/profile');
DynPred=h5read('data/WL_L1.hdf5','/predict0/profile');
rPred=h5read('data/WL_L1.hdf5','/shear/seperation');
Signal=h5read('data/WL_L1.hdf5','/shear/profile');
Err=h5read('data/WL_L1.hdf5','/shear/profile_err');
%%
figure;
loglog(rPred,LumPred,'b-')
hold on;
loglog(rPred,DynPred,'r-')
plot(r1,nfw_DeltSig(r1./(1+Zbin(1)),MLbin(1)/1e10,Zbin(1),2)/(1+Zbin(1))^2,'b--');
plot(r1,nfw_DeltSig(r1./(1+Zbin(1)),MDbin(1)/1e10,Zbin(1),2)/(1+Zbin(1))^2,'r--');
% plot(r1,nfw_DeltSig(r1./(1+Zbin(1)),par1(1),Zbin(1),2,par1(3))/(1+Zbin(1))^2,'g-');
% plot(r1,nfw_DeltSig(r1./(1+Zbin(1)),par1(1),Zbin(1),2,par1(3)),'g-');
plot(r1,(lensing_rfunc(r1,par1(1),par1(2),Zbin(1),par1(3))),'k-');
errorbar(r1,s1,es1,'ks')
xlabel('R/Mpc');ylabel('$\Delta\Sigma$');
legend('Lum Stack','Dyn Stack','Lum Average','Dyn Average','Fit');
print('-depsc','/work/Projects/Lensing/outputv4/Average_vs_stack_L1.eps')
%%
LumPred=h5read('data/WL_L3.hdf5','/predict1/profile');
DynPred=h5read('data/WL_L3.hdf5','/predict0/profile');
rPred=h5read('data/WL_L3.hdf5','/shear/seperation');
Signal=h5read('data/WL_L3.hdf5','/shear/profile');
Err=h5read('data/WL_L3.hdf5','/shear/profile_err');
%%
figure();
loglog(rPred,LumPred,'b-')
hold on;
loglog(rPred,DynPred,'r-')
plot(r3,nfw_DeltSig(r3./(1+Zbin(3)),MLbin(3)/1e10,Zbin(3),2)/(1+Zbin(3))^2,'b--');
plot(r3,nfw_DeltSig(r3./(1+Zbin(3)),MDbin(3)/1e10,Zbin(3),2)/(1+Zbin(3))^2,'r--');
% plot(r1,nfw_DeltSig(r1./(1+Zbin(1)),par1(1),Zbin(1),2,par1(3))/(1+Zbin(1))^2,'g-');
% plot(r1,nfw_DeltSig(r1./(1+Zbin(1)),par1(1),Zbin(1),2,par1(3)),'g-');
plot(r3,(lensing_rfunc(r3,par3(1),par3(2),Zbin(3),par3(3))),'k-');
errorbar(r3,s3,es3,'ks')
xlabel('R/Mpc');ylabel('$\Delta\Sigma$');
legend('Lum Stack','Dyn Stack','Lum Average','Dyn Average','Fit');
print('-depsc','/work/Projects/Lensing/outputv4/Average_vs_stack_L3.eps')
%%
figure;
x=10:0.1:14
nl=hist(log10(lummass(1:n1)),x)
nd=hist(log10(dynmass(1:n1)),x)
stairs(x,nl,'r-')
hold on
stairs(x,nd,'g-')
% loglog(lummass(1:n1),dynmass(1:n1),'.')
legend('LumMass','DynMass')
xlabel('log(M/Msun)')
ylabel('Count')
title('Mass Distribution for lowest LumMass bin')
print('-depsc','/work/Projects/Lensing/outputv4/MassDitribtuionL1.eps')