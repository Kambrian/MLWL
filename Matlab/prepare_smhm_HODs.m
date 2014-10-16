data=importdata('/work/Projects/Lensing/code/v6.0/Matlab/cosmology/dndlogM_ST_3778.dat',' ',1);
hmf=data.data;
sp=spline(log10(hmf(:,1)),log(hmf(:,3)));
halomfunc=@(logm) exp(ppval(sp,logm)); %dn/dlog10(M)/dV, for halo mass func
clear data
newMF=load('/work/Projects/Lensing/data/StellarMassFunc/massfun.dr72bbright0.modelmass');
%
MvInMvb=0.893;
MvcInMvb=0.733;
halo2star=@(M,A,M0,alpha,beta,gamma) A./((M/M0).^alpha+(M/M0).^beta).^gamma.*M;
halo2starP=@(M,P) halo2star(M,P(1),P(2),P(3),P(4),P(5)); %Msun/h^2 for Mstar, Msun/h for Mhalo
dlnMhdlnMs=@(Mh,P) 1./(1-P(5).*(P(3)*(Mh/P(2)).^P(3)+P(4)*(Mh/P(2)).^P(4))./((Mh/P(2)).^P(3)+(Mh/P(2)).^P(4)));
Pwang4=[2*10^10.23/3.43e11*0.73^2*MvcInMvb,3.43e11/MvcInMvb,1-2.56,1-0.34,1];%DR7, 2013
Pguo=[0.129*0.73*MvcInMvb,10^11.4*0.73/MvcInMvb,-0.926,0.261,2.440]; %close to Yang at >M0 ; M200c
Pling=[2*10^-1.73,10^11.70,-1.16,0.71,1];%M200b

z=0.;
parmosterZ=[0.0351-0.0247*z/(1+z),11.590+1.195*z/(1+z),1.376-0.826*z/(1+z),0.608+0.329*z/(1+z)];
PmosterZ0=[2*parmosterZ(1)*0.72*MvcInMvb,10^parmosterZ(2)*0.72/MvcInMvb,-parmosterZ(3),parmosterZ(4),1];%Moster13, M200c

z=0.3;
parHudsonZ=[0.036+(z-0.5)*0.018, 12.40+(z-0.5)*0.25, 0.68, 0.8];
PHudson=[2*parHudsonZ(1)*0.7*MvcInMvb,10^parHudsonZ(2)*0.7/MvcInMvb,-parHudsonZ(3),parHudsonZ(4),1];%Hudson13, M200c
%
powerlaw_err=@(x,lnAerr,Berr,corrAB) sqrt(lnAerr^2+(log(x)*Berr).^2+2*log(x)*corrAB*lnAerr*Berr);
powerlawmodel_err=@(m,perr,c) powerlaw_err(m/1e12,perr(1)/log10(exp(1)),perr(2),c); %error in ln(y) for y=10^p(1)*(m/1e12).^p(2)
powerlawmodel=@(m,p) 10^p(1)*(m/1e12).^p(2);
%
logMmin=7;
logMmax=17;
sigma0=0.0;%additional scatter due to our SM estimation method
dndlogMs=@(logms,sigma,par) quadgk(@(logm) halomfunc(logm).*normpdf(logms-log10(halo2starP(10.^logm,par)),0,sqrt(sigma^2+sigma0^2)),logMmin,logMmax);
dndlogMsdlogMh=@(logms,logmh,sigma,par) halomfunc(logmh).*normpdf(logms-log10(halo2starP(10.^logmh,par)),0,sqrt(sigma^2+sigma0^2));
medianMh=@(logMs,sigma,par) quadgk(@(logmh) dndlogMsdlogMh(logMs,logmh,sigma,par)./dndlogMs(logMs,sigma,par).*logmh, logMmin, logMmax);
halo2smfunc=@(mh,par) halomfunc(log10(mh)).*abs(dlnMhdlnMs(mh,par));
%--
star2halo=@(M,M1,M0,beta,delta,gamma) 10.^(log10(M1)+beta*log10(M/M0)+(M/M0).^delta./(1+(M/M0).^-gamma)-1/2); %Berhoozi10
parCosmos=[10^12.520*0.72, 10^10.916*0.72^2, 0.457, 0.566, 1.53]; %Leauthaud12
star2haloP=@(M,P) star2halo(M,P(1),P(2),P(3),P(4),P(5));
x=logspace(0,12);y=star2haloP(x,parCosmos);
halo2starB=@(M) pchip(y,x,M);
dndlogMsB=@(logms,sigma) quadgk(@(logm) halomfunc(logm).*normpdf(logms-log10(halo2starB(10.^logm)),0,sqrt(sigma^2+sigma0^2)),logMmin,logMmax);
dndlogMsdlogMhB=@(logms,logmh,sigma) halomfunc(logmh).*normpdf(logms-log10(halo2starB(10.^logmh)),0,sqrt(sigma^2+sigma0^2));
medianMhB=@(logMs,sigma) quadgk(@(logmh) dndlogMsdlogMhB(logMs,logmh,sigma)./dndlogMsB(logMs,sigma).*logmh, logMmin, logMmax);
%--

ms=logspace(8,12,20);
ywang4=ms;sigma_wang4=0.17; %Wang13
yling=ms;sigma_ling=0.22;
ymosterZ0=ms;sigma_mosterZ0=0.1;
yguo=ms;sigma_guo=0.0;
ycosmos=ms; sigma_cosmos=0.206;
yhudson=ms; sigma_hudson=0.15;
for i=1:numel(ms)
    ywang4(i)=medianMh(log10(ms(i)),sigma_wang4,Pwang4);
    yling(i)=medianMh(log10(ms(i)),sigma_ling,Pling);
    ymosterZ0(i)=medianMh(log10(ms(i)),sigma_mosterZ0,PmosterZ0);
    yhudson(i)=medianMh(log10(ms(i)),sigma_hudson,PHudson);
    yguo(i)=medianMh(log10(ms(i)),sigma_guo,Pguo); 
    ycosmos(i)=medianMhB(log10(ms(i)),sigma_cosmos);
end