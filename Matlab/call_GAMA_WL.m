macro=default_WLparam();
macro.LENSCATID=0; % catalogue used as lens: 0 means real GAMA, 1 means mockGAMA, <0 means random catalogue, >1 means random permutation catalogue
% macro.RMIN=
macro.RMAX=10;   %rmax for radial bin, in Mpc/h; if <0, set rmax=-RMAX*rvir;
macro.PROXY='lummass';
macro.ZGAMA=0; %cc2
%macro.RTRIM=0.5;
disp(macro)
GAMA_WL(macro);
% macro.PROXY='dynmass';
% disp(macro)
% GAMA_WL(macro);
% macro.PROXY='bcgmag';
% disp(macro)
% GAMA_WL(macro);
% macro.PROXY='comvsize';
% disp(macro)
% GAMA_WL(macro);
