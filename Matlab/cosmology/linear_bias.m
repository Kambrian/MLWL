function bias=linear_bias(M,z,sig8)
% linear bias b=\ksi_hm/\ksi_hg
% M: in units 10^10M_sun/h

OmegaM=0.3;
h=0.7;
fb=0.168;

a=0.707;b=0.35;c=0.8;
% a=0.707;b=0.5;c=0.6;
sa=sqrt(a);

if numel(M)==1
    M=repmat(M,size(z));
end
if numel(z)==1
    z=repmat(z,size(M));
end

if numel(M)~=numel(z)
    error('M and z must have equal size or be scalar');
end

w=collapse_barrier(1./(1+z),OmegaM);

bias=M;
for i=1:numel(M)
v=w(i)/sqrt(mass_variance(M(i),OmegaM,h,fb,sig8));
bias(i)=1+1/sa/1.686*(sa*a*v^2+b*sa*(a*v^2)^(1-c)-(a*v^2)^c/((a*v^2)^c+b*(1-c)*(1-c/2)));
end

