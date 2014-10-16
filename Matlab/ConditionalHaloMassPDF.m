logms=linspace(8,12);
[xm,ym,dyx]=linhist(log10(grp.IterCenSM(grp.Mult>2)),logms);
plot(xm,dyx./sum(ym),'o');
hold on;
plot(xm,normpdf(xm,10.72,0.5224/sqrt(2)),'-')
%%
meanh=13.2;
sigh=0.705/sqrt(2);
sigstar=0.45/sqrt(2);
pdf_star=@(logms) quadgk(@(logmh) joint_PDF_mh_mstar(logmh,logms,meanh,sigh,sigstar), -10, 20);

p=logms;
for i=1:numel(logms)
    p(i)=pdf_star(logms(i));
end
plot(logms,p,'r-');

