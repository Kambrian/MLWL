cd /work/Projects/Lensing/data
%%
galv6=fits_load_bintable('/work/Projects/Lensing/data/201308/TilingCatv40_kcorr_z00v03.fits');
hubble=1.0; %assume H=100, so all the magnitudes and luminosities should retain h in the units.
dm=@(z) 5*log10(lum_dist(0.3,z)*1e6/hubble/10); %should use asinh() DM??
kcorr=@(p1,p2,p3,p4,p5,z) p1.*z.^4+p2.*z.^3+p3.*z.^2+p4.*z+p5;
zgal=galv6.Z_TONRY;
DM=dm(zgal);DM(zgal<=0)=0;
% Umodel_abs=galv6.U_MODEL-DM-kcorr(galv6.PCOEFF_U_1,galv6.PCOEFF_U_2,galv6.PCOEFF_U_3,galv6.PCOEFF_U_4,galv6.PCOEFF_U_5,zgal);
Gmodel_abs=galv6.G_MODEL-DM-kcorr(galv6.PCOEFF_G_1,galv6.PCOEFF_G_2,galv6.PCOEFF_G_3,galv6.PCOEFF_G_4,galv6.PCOEFF_G_5,zgal);
Rmodel_abs=galv6.R_MODEL-DM-kcorr(galv6.PCOEFF_R_1,galv6.PCOEFF_R_2,galv6.PCOEFF_R_3,galv6.PCOEFF_R_4,galv6.PCOEFF_R_5,zgal);
% Imodel_abs=galv6.I_MODEL-DM-kcorr(galv6.PCOEFF_I_1,galv6.PCOEFF_I_2,galv6.PCOEFF_I_3,galv6.PCOEFF_I_4,galv6.PCOEFF_I_5,zgal);
% Zmodel_abs=galv6.Z_MODEL-DM-kcorr(galv6.PCOEFF_Z_1,galv6.PCOEFF_Z_2,galv6.PCOEFF_Z_3,galv6.PCOEFF_Z_4,galv6.PCOEFF_Z_5,zgal);
%%
load G3Cv4up8/G3Cv4up8

%%
file='GAMA-mag.hdf5';
dsize=size(gal.Rpetro');
h5create(file, '/GAMAI/rPetro', dsize, 'datatype','single');
h5write(file,'/GAMAI/rPetro', single(gal.Rpetro'));
h5writeatt(file, '/GAMAI/rPetro', 'Info', 'r band apparent petrosian magnitude. extinction corrected. GAMA selection magnitude. use rPetro<19.4 for a flux-limited sample.');

h5create(file, '/GAMAI/Redshift', dsize, 'datatype','single');
h5write(file, '/GAMAI/Redshift', single(gal.Zspec'));
h5writeatt(file, '/GAMAI/Redshift', 'Info', 'Redshift');
h5create(file, '/GAMAI/RModel', dsize, 'datatype','single');
h5write(file, '/GAMAI/RModel', single(gal.Rmodel_abs'));
h5writeatt(file, '/GAMAI/RModel', 'Info', 'r band absolute model magnitude. k-corrected.');
h5create(file, '/GAMAI/GModel', dsize, 'datatype','single');
h5write(file, '/GAMAI/GModel', single(gal.Gmodel_abs'));
h5writeatt(file, '/GAMAI/GModel', 'Info', 'g band absolute model magnitude. k-corrected.');
f=gal.RA>160&gal.DEC<200;
h5create(file, '/GAMAI/DeepFieldFlag', dsize, 'datatype','int32');
h5write(file, '/GAMAI/DeepFieldFlag', int32(f'));
h5writeatt(file, '/GAMAI/DeepFieldFlag', 'Info', 'Flag to indicating whether the galaxy is in the r<19.8 field. use rPetro<19.8 and DeepFieldFlag==1 to select a deeper flux-limited sample than the default one.')
h5writeatt(file,'/GAMAI', 'Info', 'GAMA-I galaxy magnitudes. Galaxies from the G3Cv4 catalogue, while absolute magnitudes from TilingCatv40_kcorr_z00v03.fits. Missing values are filled with NaN.');


f=galv6.SURVEY_CLASS>=3&galv6.NQ>=3;
dsize=[1,sum(f)];
h5create(file, '/GAMAII/rPetro', dsize, 'datatype','single');
h5write(file, '/GAMAII/rPetro', single(galv6.R_PETRO(f))');
h5writeatt(file, '/GAMAII/rPetro', 'Info', 'r band apparent petrosian magnitude. probably NOT extinction corrected.');
h5create(file, '/GAMAII/RModel', dsize, 'datatype','single');
h5write(file, '/GAMAII/RModel', single(Rmodel_abs(f))');
h5writeatt(file, '/GAMAII/RModel', 'Info', 'R band absolute model magnitude. k-corrected.');
h5create(file, '/GAMAII/GModel', dsize, 'datatype','single');
h5write(file, '/GAMAII/GModel', single(Gmodel_abs(f))');
h5writeatt(file, '/GAMAII/GModel', 'Info', 'g band absolute model magnitude. k-corrected.');
h5create(file, '/GAMAII/Redshift', dsize, 'datatype','single');
h5write(file, '/GAMAII/Redshift', single(galv6.Z_TONRY(f))');
h5writeatt(file, '/GAMAII/Redshift', 'Info', 'Redshift. (Z_TONRY)');
h5writeatt(file,'/GAMAII','Info','GAMA-II galaxy magnitudes extracted from TilingCatv40_kcorr_z00v03.fits, with SURVEY_CLASS >= 3 and nQ>=3 selection, in the three equatorial regions.');
