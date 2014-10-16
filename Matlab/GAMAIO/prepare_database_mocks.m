%%
bower=importdata('additional/smhm_bower.millimil.dat');
bower.Rmag=bower.data(:,5);
bower.LumHalpha=bower.data(:,4)*1e33; %J/s/h2
bower.sfr=bower.LumHalpha/1.27e34*1e9; %Msun/h2/Gyr
bower.sm=bower.data(:,3)*0.7; %Msun/h2
bower.mh=bower.data(:,1)*8.6e8; %Msun/h
bower.md=bower.data(:,2); %Msun/h, the dhalo mass used in the galform model, required to be non-decreasing so is no smaller than mh.
%%
bower=importdata('additional/smhm_guo13.dat');
bower.Rmag=bower.data(:,5);
bower.sfr=bower.data(:,4)*0.7^2*1e9; %Msun/h2/Gyr
bower.sm=bower.data(:,3)*0.7*1e10; %Msun/h2
bower.mh=bower.data(:,1)*1e10; %Msun/h, 200mean
bower.md=bower.data(:,2)*8.6e8; %Msun/h, the subhalo bound mass
%% assign redshift; not necessary, produces little difference
dist2z=@(d) fzero(@(z) comoving_dist(0.3,0.7,z)-d,0.2);
h=0.7;
dm=@(z) 5*log10(lum_dist(0.3,z)*1e6/h/10);

zmax=0.5;
rmax=comoving_dist(0.3,0.7,zmax);
Vmax=rmax.^3;
bower.dist_comov=(rand([numel(bower.Rmag),1])*Vmax).^(1/3); %random comoving dist
bower.z=bower.dist_comov;
for i=1:numel(bower.Rmag)
    bower.z(i)=dist2z(bower.dist_comov(i));
end
bower.rmag=bower.Rmag+dm(bower.z);
%%
guo=bower;
save('additional/smhm_guo13.mat','guo');
%%
font=bower;
save('additional/smhm_font.milli.mat','font');
%%
save('additional/smhm_bower.milli.mat','bower'); %little difference in our observable between millimill and mill.