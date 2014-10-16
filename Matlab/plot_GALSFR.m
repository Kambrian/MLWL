f1=gal.logSFR>-10&gal.logSFR<0.38;
f2=gal.logSFR>0.38&gal.logSFR<1.0;
f3=gal.logSFR>1.0&gal.logSFR<10;
Zbin=[mean(gal.Zspec(f1)),mean(gal.Zspec(f2)),mean(gal.Zspec(f3))];

f=~isnan(gal.SM);
ff=gal.LumMass>0;
SFRbin=[median(10.^gal.logSFR(f1)),median(10.^gal.logSFR(f2)),median(10.^gal.logSFR(f3))];
SMbin=[median(10.^gal.SM(f1&f)),median(10.^gal.SM(f2&f)),median(10.^gal.SM(f3&f))];
LMbin=[median(gal.LumMass(f1&ff)),median(gal.LumMass(f2&ff)),median(gal.LumMass(f3&ff))];

efun=@(x) prctile(x,[(100-68.3)/2;100-(100-68.3)/2]);
eSFRbin=[efun(10.^gal.logSFR(f1)),efun(10.^gal.logSFR(f2)),efun(10.^gal.logSFR(f3))];
eSMbin=[efun(10.^gal.SM(f1&f)),efun(10.^gal.SM(f2&f)),efun(10.^gal.SM(f3&f))];
eLMbin=[efun(gal.LumMass(f1&ff)),efun(gal.LumMass(f2&ff)),efun(gal.LumMass(f3&ff))];
%%
figure;
f=~isnan(gal.SM)&gal.logSFR>-10&gal.logSFR<10;
plot(gal.SM(f),gal.logSFR(f),'.')
%%
OmegaM=0.25;H0=100;G=43.0071;
rhob=OmegaM*3*H0^2/8/pi/G;
binname='GALSFR';

myfigure;
[h,r1,s1,es1,cov1]=plotgama_prof('GSFR1.ZGAMA.V4','r.');
hold on;
[h,r2,s2,es2,cov2]=plotgama_prof('GSFR2.ZGAMA.V4','go');
[h,r3,s3,es3,cov3]=plotgama_prof('GSFR3.ZGAMA.V4','bs');
%%
[cf1,gof1,out1]=createFit(r1,s1,es1.^-2,Zbin(1));
ci1=confint(cf1,0.683);
par1=coeffvalues(cf1);
[cf2,gof2,out2]=createFit(r2,s2,es2.^-2,Zbin(2));
ci2=confint(cf2,0.683);
par2=coeffvalues(cf2);
[cf3,gof3,out3]=createFit(r3,s3,es3.^-2,Zbin(3));
ci3=confint(cf3,0.683);
par3=coeffvalues(cf3);

mfit=[par1(1);par2(1);par3(1)]*1e10;
bfit=[par1(2);par2(2);par3(2)];
errm=[ci1(:,1),ci2(:,1),ci3(:,1)]*1e10;
errb=[ci1(:,2),ci2(:,2),ci3(:,2)];
%%
J=out1.Jacobian;
ok= isfinite(r1) & isfinite(s1) & isfinite(es1);  % to account for number of observations in fitting
err1=J\cov1(ok,ok)/J';
err=inv(J'*J);
J=out2.Jacobian;
ok= isfinite(r2) & isfinite(s2) & isfinite(es2);
err2=J\cov2(ok,ok)/J';
J=out3.Jacobian;
ok= isfinite(r3) & isfinite(s3) & isfinite(es3);
err3=J\cov3(ok,ok)/J';
em=sqrt([err1(1),err2(1),err3(1)])';
eb=sqrt([err1(4),err2(4),err3(4)])';
cmb=[err1(2),err2(2),err3(2)]';  %correlation between m and b
cmb=cmb./em./eb;
em=em*1e10;
% plot(r1,log10(nfw_DeltSig(r1./(1+Zbin(1)),Mbin(1)/1e10,Zbin(1),2)/(1+Zbin(1))^2),'r-');
% plot(r2,log10(nfw_DeltSig(r2./(1+Zbin(2)),Mbin(2)/1e10,Zbin(2),2)/(1+Zbin(2))^2),'g-');
% plot(r3,log10(nfw_DeltSig(r3./(1+Zbin(3)),Mbin(3)/1e10,Zbin(3),2)/(1+Zbin(3))^2),'b-');
% plot(r1,log10(2*rhob*ppval(DSig_interp_spline,r1).*growth_factor(0.3,Zbin(1)).^2),'r--');
% plot(r2,log10(2.5*rhob*ppval(DSig_interp_spline,r2).*growth_factor(0.3,Zbin(2)).^2),'g--');
% plot(r3,log10(4*rhob*ppval(DSig_interp_spline,r3).*growth_factor(0.3,Zbin(3)).^2),'b--');

plot(r1,signlog(lensing_rfunc(r1,par1(1),par1(2),Zbin(1))),'r-');
plot(r2,signlog(lensing_rfunc(r2,par2(1),par2(2),Zbin(2))),'g-');
plot(r3,signlog(lensing_rfunc(r3,par3(1),par3(2),Zbin(3))),'b-');
legend('SFR<2.4','2.4<SFR<10','SFR>10');
print('-depsc','wlSig_GALSFR.ZGAMA.eps');
%%
myfigure;
subplot(3,1,1);
ploterr(mfit,SFRbin,em,{eSFRbin(1,:),eSFRbin(2,:)},'.','logxy');
subplot(3,1,2);
ploterr(mfit,SMbin,em,{eSMbin(1,:),eSMbin(2,:)},'.','logxy');
subplot(3,1,3);
ploterr(mfit,LMbin,em,{eLMbin(1,:),eLMbin(2,:)},'.','logxy');
