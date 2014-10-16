function A=dynmass_calibration(cat)
%the prefactor A for dynamical mass calibration
%Robotham et al 2011, Table 3, r_AB<19.4.

mult=[2,5,10,20,1000];[~,ibin]=histc(cat.Mult,mult);
z=[0,0.1,0.2,0.3,0.5];[~,jbin]=histc(cat.Zfof,z);

C=[19.0,10.8,12.0,12.6;
    19.5,10.5,11.1,10.4;
    21.5,10.3,8.6,8.3;
    17.4,6.1,5.4,5.6;
    ];

A=ibin;
for i=1:numel(ibin)
A(i)=C(ibin(i),jbin(i));
end