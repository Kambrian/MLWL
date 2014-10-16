function macro=default_WLparam()

macro.PROXY='halomass';  %mass proxy for stacking
macro.LENSCATID=0; % catalogue used as lens: 0 means real GAMA, 1 means mockGAMA, <0 means random catalogue, >1 means random permutation catalogue
macro.STACK_COMOV=1;  %comoving stacking
macro.RADIAL_LOGBIN=1;  %radial bin in logscale
macro.BCG_CENTER=1;     %0: use Cen; 1 :use BCG as center; 2: use IterCen
macro.MULT_MIN=3;      %minimum multiplicity of groups
macro.PHOTOZ_DENS_ERR2_MAX=0;  %maximum photoz uncertainty in SigmaCrit, set to 0 means do not calc nor include this error
macro.PHOTOZ_DENS_ERR_INC=0; %include photoz error in calc; PHOTOZ_DENS_ERR2_MAX must be >0 for this to take effect
macro.RMIN=0.02;  %rmin for radial bin, in Mpc/h
macro.RMAX=10;   %rmax for radial bin, in Mpc/h; if <0, set rmax=-RMAX*rvir;
macro.ZGAMA=0;  % >0: use gama's spec-z to update source cat's photoz; 
                                % abs()=2: use zcc2 to update photoz's first;
                                % abs()=3: use zd1
                                %            4:       sdss_tmpl
                                %            1:  GAMAspec-Z only
                                %            -1 or 0: nothing, keep ZEBRA
% macro.ZMAX;  %zmax of groups
% macro.BCG_Z=1; %use BCG redshift as center redshift; ===== deprecated because it has NO discernable difference from MedianZ
% macro.ZDIFF=0.01; %minimum redshift difference for source from lens
if macro.PHOTOZ_DENS_ERR2_MAX==0
    macro.PHOTOZ_DENS_ERR_INC=0;  %disable zerr
end