function [Mout,Cout]=convert_NFWmc(Min,Cin,vir_in,vir_out,z)
% convert the mass concentration from one definition to another
% Min, Cin: input Mass(10^10Msun/h) and concentration
% vir_in: input mass virial definition. 0: spherical collapse; 1: 200*crit; 2: 200*mean
% vir_out: mass definition to convert to
% z: redshift
% Mout, Cout: output mass and concentration

G=43.0071;
HUBBLE0=100;  %km/s/(Mpc/h)
if nargin<5
    z=0.;
end
Omega0=0.3;OmegaLambda=0.7;scaleF=1./(1+z);
Hz=HUBBLE0 * sqrt(Omega0 ./scaleF.^3+ (1 -Omega0 -OmegaLambda) ./ scaleF.^2 +OmegaLambda);

Hratio=Hz/HUBBLE0;
OmegaZ=Omega0./scaleF.^3./Hratio.^2;

% virFin=vir_in;virFout=vir_out;
virFin=vir_factor(vir_in);
virFout=vir_factor(vir_out);
    function virialF=vir_factor(virtype)
        switch virtype   %virial factor with respect to critical density
            case 0
                virialF=18.0*pi^2+82.0*(OmegaZ-1)-39.0*(OmegaZ-1).^2;
            case 1
                virialF=200.;
            case 2
                virialF=200*OmegaZ;
            otherwise
                error('virialtype must be 0/1/2')
        end
    end

rhoc=(3*Hz.^2)/(8*pi*G);
rhos=virFin/3.*Cin.^3./(log(1+Cin)-Cin./(1+Cin)).*rhoc;
Rvin=(Min./(4*pi/3*virFin.*rhoc)).^(1/3);
rs=Rvin./Cin;

Cout=Cin;
for i=1:numel(Cin)
    Cout(i)=fzero(@(x) virFout/3.*x.^3./(log(1+x)-x./(1+x))-rhos(i)./rhoc,Cin(i));
end
Rvout=Cout.*rs;
Mout=4*pi/3*rhoc*virFout.*Rvout.^3;
end