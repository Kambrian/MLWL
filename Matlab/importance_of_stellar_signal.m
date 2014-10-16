MvInMvb=0.893;
MvcInMvb=0.733;
halo2star=@(M,A,M0,alpha,beta,gamma) A./((M/M0).^alpha+(M/M0).^beta).^gamma.*M;
halo2starP=@(M,P) halo2star(M,P(1),P(2),P(3),P(4),P(5)); %Msun/h^2 for Mstar, Msun/h for Mhalo
Pguo=[0.129*0.73*MvcInMvb,10^11.4*0.73/MvcInMvb,-0.926,0.261,2.440]; %close to Yang at >M0 ; M200c

mh=10.^(10:14);
ms=halo2starP(mh,Pguo);
r=logspace(log10(0.02),1);
c='rgbck';
myfigure;
for i=1:numel(mh)
%     loglog(r,ms(i)/1e10/4/pi./r.^2,'r');
%     hold on;
%     loglog(r,nfw_DeltSig(r,mh(i)/1e10,0,2),'g');
loglog(r,ms(i)/1e10/4/pi./r.^2./nfw_DeltSig(r,mh(i)/1e10,0,2),c(i),'displayname',printexp10(mh(i)));
hold on;
end
legend('show')
xlabel('R[Mpc/h]');ylabel('$\Delta\Sigma_{\star}/\Delta\Sigma_{halo}$');
print('-depsc','/work/Projects/Lensing/outputv4/stellarmass_importance.eps');