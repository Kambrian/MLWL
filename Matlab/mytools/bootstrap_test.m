%compare the distribution of  bootstrap sampled gamma-order moment with that from nboot real realizations
% ==conclusion: bootstrap does a good job in general; even though the
% resampled data may have the same bias as the original data, the scatter
% is correctly captured
% in the extrem case when you have a infinitely large sample size, your
% bootstrap would be identical to the ideal sample, and hence identical to
% the underlying distribution
nsamp=500;
nboot=50*nsamp;
gamma=4;  
xt=normrnd(0,1,nboot,nsamp);
x=xt(1,:);

xx=zeros(nboot,numel(x));
for i=1:nboot
    xx(i,:)=x(randi(nsamp,1,nsamp));
end

% figure;hist(x);
% figure;hist(xx(:))

% xm=mean(xx,2);
% xtm=mean(xt,2);
xm=mean(xx.^gamma,2);
xtm=mean(xt.^gamma,2);
[std(xm),std(xtm)]

figure;
[xi,yi,dyx]=linhist(xm,20);
plot(xi,dyx/sum(yi),'.');
hold on;
% mu=mean(x);
% sig=1/sqrt(nsamp);
% xp=mu+(-3*sig:0.1*sig:3*sig);
% plot(xp,normpdf(xp,mu,sig),'r');

[xi,yi,dyx]=linhist(xtm,20);
plot(xi,dyx/sum(yi),'o');