data=importdata('LumFunc/lf_gama.cat.2.0.01_z_0.51_swml_z.0.05.0.10.best.txt',' ',2);
lumfunc_data=data.data;
lumfunc=@(m,a,mc,fai) 0.4*log(10)*fai*10.^(0.4*(mc-m)*(1+a)).*exp(-10.^(0.4*(mc-m)));
dm=@(z) 5*log10(lum_dist(0.3,z)*1e6/10); %DM+5log10(h)
pkcorr=[2.3843,3.5902,0.5237,1.0226,0.2085]; %k correction polynomial
kecorr=@(z) polyval(pkcorr,z-0.2)+1.75*z;
Mlim=@(z) 19.8-dm(z)-kecorr(z);
figure;
semilogy(lumfunc_data(:,1),lumfunc_data(:,[5:end,3])./repmat(lumfunc_data(:,3),[1,6]));
hold on;
plot(repmat(Mlim(0.0:0.1:0.4),[2,1]),repmat([1e-1;1e1],[1,5]),'--');
% plot(lumfunc_data(:,1),lumfunc(lumfunc_data(:,1),-1.26,-20.73,0.9/1e2),'ko');
%%
mask=lumfunc_data(:,5)>0;
fun = fit(lumfunc_data(mask,2),(lumfunc_data(mask,5)),'linearinterp');
figure;
plot(lumfunc_data(:,1),(lumfunc_data(:,5)),'o');hold on;
plot(fun);
yscale('log')
%%
colors='rgbcmk';
figure;
cols=[5:9];
zbin=0.001:0.1:0.5;
x=lumfunc_data(:,1);
x(x==0)=lumfunc_data(x==0,2);
x0=-26;
% fun0=fit(x, lumfunc_data(:,3), @(fai,x) 0.4*log(10)*fai*10.^(0.4*(-20.73-x)*(1-1.26)).*exp(-10.^(0.4*(-20.73-x))));
% y0=integrate(fun0,x,x0);
fun0=fit(x, lumfunc_data(:,3), 'linearinterp');y0=integrate(fun1,x,x0);
% plot(fun0,'k--','integral');hold on;
for i=1:numel(cols)
    mask=x<Mlim(zbin(i));
    y=lumfunc_data(:,cols(i));
    y(~mask)=lumfunc_data(~mask,3);
fun = fit(x,y,'linearinterp');
y=integrate(fun,x,x0)./y0;
plot(x,y,colors(i));hold on;
% plot(fun,colors(i),'integral');hold on;
end
yscale('log');
plot(repmat(Mlim(0.001:0.1:0.5),[2,1]),repmat([1e-1;1e1],[1,5]),'linestyle','--');
%%
z=0.01:0.02:0.5;
zbin=0:0.1:0.5;
x=lumfunc_data(:,1);
x(x==0)=lumfunc_data(x==0,2);
[~,ibins]=histc(z,zbin);
nav=z;
for i=1:numel(z)
    y=lumfunc_data(:,cols(ibins(i)));
    fun = fit(x,y,'linearinterp');
    nav(i)=integrate(fun,Mlim(z(i)),x0);
end
figure;
loglog(z,nav,'.');
    
    