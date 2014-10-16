clear global
macros.PROXY='halomass';  %mass proxy for stacking
macros.LENSCATID=0; % catalogue used as lens: 0 means real GAMA, 1 means mockGAMA, <0 means random catalogue
macros.STACK_COMOV=1;  %comoving stacking
macros.RADIAL_LOGBIN=1;  %radial bin in logscale
macros.BCG_CENTER=1;     %use BCG as center
macros.MULT_MIN=3;      %minimum multiplicity of groups
macros.PHOTOZ_DENS_ERR2_MAX=0;  %maximum photoz uncertainty in SigmaCrit, set to 0 means do not calc nor include this error
macros.PHOTOZ_DENS_ERR_INC=0; %include photoz error in calc; PHOTOZ_DENS_ERR2_MAX must be >0 for this to take effect
macros.RMIN=0.02;  %rmin for radial bin, in Mpc/h
macros.RMAX=-3;   %rmax for radial bin, in Mpc/h; if <0, set rmax=-RMAX*rvir;
% macros.ZMAX;  %zmax of groups
macros.ZGAMA=1;  %use gama's spec-z to update source cat's photoz
if macros.PHOTOZ_DENS_ERR2_MAX==0
    macros.PHOTOZ_DENS_ERR_INC=0;  %disable zerr
end

ltyp={'ro-','ro--';'g>-','g>--';'bs-','bs--';'kd-','kd--';};
ldtyp={'ro:';'g>:';'bs:';'kd:';};
lltyp={'ro-';'g>-';'bs-';'kd-';};
mtyp={'ro';'g>';'bs';'kd';};
color={'r';'g';'b';'k'};
nbinvar=[3,6,12,15];
%
macro=macros;
% macro.MULT_MIN=2; 
data=gama_rebin(macro,nbinvar);
macro2=macro;
% macro2.LENSCATID=-1;
% macro2.BCG_CENTER=1;
% macro2.MULT_MIN=3; 
macro2.ZGAMA=2;
data2=gama_rebin(macro2,nbinvar);
% the inclusion of photo-z error overestimates the fluctuation
% however, the shear-error along tend to over-estimate the error in the
% inner part for the largest mass bin
myfigure;
for l=1:4
h1=plot(data.rvar{l},signlog(data.svar{l}),lltyp{l},'markerfacecolor',color{l});hold on;
% plot(data.rvar{l},signlog(data.esvar{l}),mtyp{l},'markerfacecolor',color{l});hold on;
end

for l=1:4
h2=plot(data2.rvar{l},signlog(data2.svar{l}),ldtyp{l},'markersize',8);hold on;
% plot(data2.rvar{l},signlog(data2.esvar{l}),mtyp{l},'markersize',8);hold on;
end
plot([1e-2,1e1],[-4.5,-1.5],'c:',[1e-2,1e1],[4.5,1.5],'c:');
set(gca,'xscale','log');
legend([h1,h2],'ZEBRA','CC2');
print('-depsc','wlSig_cmp_ZEBRA_CC2.eps');

myfigure;
for l=1:4
% plot(data.rvar{l},signlog(data.svar{l}),lltyp{l},'markerfacecolor',color{l});hold on;
plot(data.rvar{l},data2.esvar{l}./data.esvar{l},mtyp{l},'markerfacecolor',color{l});hold on;
end
xlabel('r');ylabel('$\sigma_{new}/\sigma_{old}$');
set(gca,'xscale','log');
print('-depsc','wlErr_cmp_ZEBRA_CC2.eps');