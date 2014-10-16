clear
load GAMA_SRC_matchZ.mat
load shearmap
%%
zp=[];
zs=[];
for i=1:3
    zp=[zp;e{i}(source_idmatch{i},6)];
    zs=[zs;gama_zmatch{i}(:,1)];
end
%%
figure;
plot(zs,zp,'.')
hold on;
% plot([0,0.7],[0,0.7],'r')
% xlim([0,.7]);
% ylim([0,.7]);
%%
nbin=20;
[xmed,ymed,ylim,xm,ym,ysig,count]=skeleton(zs,zp,nbin,0.683);
plot(xmed,ymed,'r--');
plot(xm,ym,'k-');
plot(xm,ym+ysig,'k-',xm,ym-ysig,'k-');
plot(xm,ym+sigz(xm),'g',xm,ym-sigz(xm),'g')
%%
figure;
plot(xm,ysig)
hold on;
plot(xm,sigz(xm),'r')
