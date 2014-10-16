function a=numprof_sup_factor()

sigmz=0.1;
zc=0.15;
OmegaM=0.3;H=100;c=3e5;
h=@(z) H*sqrt(OmegaM .*(1+z).^3+(1-OmegaM));
g=@(z) exp(-z.^2./2/sigmz^2)/sqrt(2*pi)/sigmz;
d=@(z) exp(-8.8*z);
p=@(z) g(z-zc).*d(z);

a=quad(p,zc+0.1,1.5)/quad(@(z) hc(z).*d(z),zc+0.1,1.5);

function q=hc(z)
q=z;
for i=1:numel(z)
    q(i)=quad(@(zp) c*g(z(i)-zp)/h(z(i)),0,2);
end
end

end