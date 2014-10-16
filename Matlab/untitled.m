myfigure;
for i=1:4
    semilogx(a.rvar{i},a.y{i}./b.y{i},lltyp{i});
    hold on;
end
set(gca,'xscale','log');
% legend([h1,h2],'mock','gama');
xlabel('r/(Mpc/h)');ylabel('$n_{new}/n_{old}$');
cd /home/kam/Projects/Lensing/output
% print('-depsc','n_prof_zcmp.eps');
%%
catall=structcat([cat{1},cat{2},cat{3}]);
f=catall.Mult>2;
% nbin=[30,50,50,10];

myfigure;
z=zeros(50,4);n=zeros(50,4);
for i=1:4
zmean(i)=mean(catall.MedianZ(f&catall.mid==i));
[z(:,i),n(:,i)]=linhist(catall.MedianZ(f&catall.mid==i),0:0.01:0.5,'stairs',lltyp{i});
hold on;
end
data=[(0:0.01:0.49)',z,n];
% save gamaWL_zdist.txt -ascii data
% z(isnan(z))=-1;

