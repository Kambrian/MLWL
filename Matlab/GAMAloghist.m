function ax=GAMAloghist(grpmock,grp,nbin,masstype)
%nbin=15;
[grp09,grp12,grp15]=split_mockcat(grp);

ax=ones(2,1);
ax(1)=subplot(2,1,1);
[xm,ym,dym]=loghist(grpmock.(masstype),nbin,'stairs','k--');
hold on;
[x1,y1,dy1]=loghist(grp09.(masstype),nbin,'stairs','r');
[x2,y2,dy2]=loghist(grp12.(masstype),nbin,'stairs','g');
[x3,y3,dy3]=loghist(grp15.(masstype),nbin,'stairs','b');
[x,y,dy]=loghist(grp.(masstype),nbin,'stairs','k');
ylabel('dN');
legend('Mock','G09','G12','G15','GAMA','location','northwest');
% axis tight;
% yl=get(gca,'ylim');
% set(gca,'ylim',[yl(1)/2,yl(2)*2]);

ax(2)=subplot(2,1,2);
dym=dym/sum(ym);dy1=dy1./sum(y1);dy2=dy2./sum(y2);dy3=dy3/sum(y3);dy=dy/sum(y);
loglog(xm,dym,'k.--');
hold on;
loglog(x1,dy1,'r.-');
loglog(x2,dy2,'g.-');
loglog(x3,dy3,'b.-');
loglog(x,dy,'k.-');
xlabel(masstype);
ylabel('$dN/dx/N_{tot}$','interpreter','latex');
% legend('Mock','G09','G12','G15','GAMA','location','southwest');
% axis tight;
% yl=get(gca,'ylim');
% set(gca,'ylim',[yl(1)/2,yl(2)*2]);
TightenSubplots(ax);

