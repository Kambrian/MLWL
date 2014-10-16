b*lin_corr_DeltSig(r)*rhob*growth_factor(OmegaM,z)^2;
%% Output DSig, growth factor, and bias, for calculation (interpolation) of the two-halo term
%%
sig8=0.8;
z=0:0.1:0.5;
c='rgbkcmy';
logM=linspace(0,16,10);
b=zeros(numel(logM),numel(z));
figure;
for i=1:numel(z)
b(:,i)=linear_bias(10.^(logM-10),z(i),sig8);
plot(logM,b(:,i),c(i));
hold on;
end
yscale('log')
%% at most 20% difference to z=0.2 case. good enough
for i=1:numel(z)
plot(logM,b(:,i)./b(:,3),c(i));
hold on;
end
% accurate to 5% percent
% x=log10(M), y=log10(b)
% biasfit = 
% 
%      Linear model Poly5:
%      biasfit(x) = p1*x^5 + p2*x^4 + p3*x^3 + p4*x^2 + p5*x + p6
%      Coefficients (with 95% confidence bounds):
%        p1 =   4.378e-06  (3.673e-06, 5.083e-06)
%        p2 =  -3.705e-05  (-6.717e-05, -6.94e-06)
%        p3 =  -0.0003215  (-0.000785, 0.000142)
%        p4 =    0.003687  (0.0005968, 0.006777)
%        p5 =    -0.02313  (-0.0315, -0.01476)
%        p6 =    -0.07584  (-0.08276, -0.06892)
%% compute and save bias data 
z=0.2;
sig8=0.8;
logM=linspace(0,17);
b=linear_bias(10.^(logM-10),z,sig8);
% figure;semilogy(logM,b);
x=[logM;b]';
save('/work/Projects/Lensing/data/correlation/bias_logM.dat','x','-ascii');
%% compute and save correlation data
OmegaM=0.3;H0=100;G=43.0071;
rhob=OmegaM*3*H0^2/8/pi/G; %comoving matter density, 10^10Msun/h/(Mpc/h^3)

r=logspace(-3,2); %comoving radius, Mpc/h
s=lin_corr_DeltSig(r)*rhob; %comoving density, 10^10Msun/h/(Mpc/h)^2
%%
loglog(r,s);
% well fitted by polynomial, within 5 percent (0.02dex).
% y=log10(s), x=log10(r):
% corrfit = 
% 
%      Linear model Poly5:
%      corrfit(x) = p1*x^5 + p2*x^4 + p3*x^3 + p4*x^2 + p5*x + p6
%      Coefficients (with 95% confidence bounds):
%        p1 =    -0.00325  (-0.003798, -0.002701)
%        p2 =    -0.02437  (-0.02591, -0.02283)
%        p3 =    -0.05605  (-0.05912, -0.05298)
%        p4 =      -0.156  (-0.1622, -0.1498)
%        p5 =      0.3879  (0.3823, 0.3935)
%        p6 =       1.518  (1.513, 1.522)
%%
x=[log10(r);s]';
save('/work/Projects/Lensing/data/correlation/linDSig_logR_Cmv.dat','x','-ascii');
%% save growth rate
OmegaM=0.3;
z=linspace(0,1);
d=growth_factor(OmegaM,z);
x=[z;d]';
save('/work/Projects/Lensing/data/correlation/GrowthFactor_z.dat','x','-ascii');
% x=z, y=D.
% well fitted to within 0.1 percent (extremely accurate)
% growthfit = 
% 
%      Linear model Poly2:
%      growthfit(x) = p1*x^2 + p2*x + p3
%      Coefficients (with 95% confidence bounds):
%        p1 =      0.1318  (0.1313, 0.1324)
%        p2 =     -0.5195  (-0.52, -0.5189)
%        p3 =           1  (1, 1)