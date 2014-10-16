function r=comoving_function(chi,OmgTotal)
% chi: comoving radial distance, in units Mpc/h
% K: curvature param
% r:  comoving radial coordinate, in units Mpc/h

c=3e5; %light speed, km/s
% G=43007.1*1000; %gadget unit modified from kpc/h to Mpc/h
Hubble=100; %km/s/(Mpc/h)
K=(OmgTotal-1)*Hubble^2/c^2;

if abs(OmgTotal-1)<eps
    r=chi;
    return
else if K>0
    r=sin(chi*sqrt(K))/sqrt(K);
    return
    else
        r=sinh(chi*sqrt(-K))/sqrt(-K);
    end
end

