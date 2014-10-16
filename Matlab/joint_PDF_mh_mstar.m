function p=joint_PDF_mh_mstar(logmh,logmstar,meanh,sigh,sigstar)

MvInMvb=0.893;
MvcInMvb=0.733;
halo2star=@(M,A,M0,alpha,beta,gamma) A./((M/M0).^alpha+(M/M0).^beta).^gamma.*M;
halo2starP=@(M,P) halo2star(M,P(1),P(2),P(3),P(4),P(5)); %Msun/h^2 for Mstar, Msun/h for Mhalo
Pwang4=[2*10^10.23/3.43e11*0.73^2*MvcInMvb,3.43e11/MvcInMvb,1-2.56,1-0.34,1];%DR7, 2013
Pguo=[0.129*0.73*MvcInMvb,10^11.4*0.73/MvcInMvb,-0.926,0.261,2.440]; %close to Yang at >M0 ; M200c
Pmoster=[2*0.0282*0.72*MvInMvb,10^11.884*0.72/MvInMvb,-1.057,0.556,1]; %close to ling at low M; tophat Mvir
Pguofit=[exp(-2.594)*0.73*MvcInMvb,exp(26.32)/MvcInMvb,-1.743,0.592,1]; % refit guo with 4-par model
Pling=[2*10^-1.73,10^11.70,-1.16,0.71,1];%M200b
z=0.;
parmosterZ=[0.0351-0.0247*z/(1+z),11.590+1.195*z/(1+z),1.376-0.826*z/(1+z),0.608+0.329*z/(1+z)];
PmosterZ0=[2*parmosterZ(1)*0.72*MvcInMvb,10^parmosterZ(2)*0.72/MvcInMvb,-parmosterZ(3),parmosterZ(4),1];%Moster13, M200c

Gaussian=@(x,m,s) exp(-(x-m).^2./s.^2/2)./s/sqrt(2*pi);
p=Gaussian(logmh,meanh,sigh).*Gaussian(logmstar,log10(halo2starP(10.^logmh,Pwang4)),sigstar);