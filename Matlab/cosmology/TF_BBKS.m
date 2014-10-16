function tf=TF_BBKS(OmegaM,h,fb,k,~)
%BBKS (Bardeen et al. 1986) transfer function for adiabatic CDM model
% OmegaM: matter param
% h:  dimensionless hubble param
% fb: baryon fraction, OmegaB/OmegaM
% k: wavenumber, (Mpc/h)^-1
% output TF at k

Theta27=2.725/2.7; % CMB temperature parameter, Tcmb/2.7K
OmegaB=fb*OmegaM;
Gamma=OmegaM*h*exp(-OmegaB*(1+sqrt(2*h)/OmegaM)); % Sugiyama,N. 1995, ApJS, 100, 281; 
                                                                                                          % the exponetial factor corrects for baryon effect on top of 
                                                                                                          % the original BBKS fit
q=k/Gamma*Theta27^2;
tf=log(1+2.34*q)./(2.34*q).*(1+3.89*q+(16.1*q).^2+(5.46*q).^3+(6.71*q).^4).^(-0.25); % BBKS, G3

if nargin>4 %output baryon TF, BBKS, G4
    Rj=1.6/sqrt(OmegaM)/1000;  %Mpc/h
    tf=tf./(1+(k.*Rj).^2/2);
end