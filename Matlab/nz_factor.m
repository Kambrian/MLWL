function f=nz_factor(z)
f=z;
for i=1:numel(z)
f(i)=(AD_dist_flat(0.3,0,z(i))*(1+z(i)))^-2*quad(@source_zdistr,z(i)+0.1,20);
end
end