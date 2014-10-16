nbin=100;
wght=tdata(:,3);%.^-2/sum(tdata(:,3).^-2);
signal=tdata(:,1);%.*tdata(:,4);
wmin=min(wght);
wmax=max(wght);
wbin=logspace(log10(wmin),log10(wmax),nbin+1);
[n,bin]=histc(wght,sbin);
emean=zeros(nbin,1);evar=emean;weff=emean;wmean=emean;
ebrd=zeros(nbin,2);
alpha=(1-0.683)/2;
for i=1:nbin
    if n(i)
    tmp=sort(signal(bin==i));
    ebrd(i,:)=[tmp(ceil(alpha*n(i))),tmp(ceil((1-alpha)*n(i)))];
    emean(i)=mean(signal(bin==i));
    evar(i)=std(signal(bin==i));
    weff(i)=sum(signal(bin==i).*wght(bin==i))/sum(signal(bin==i));
    wmean(i)=mean(wght(bin==i));
    end
end
figure;semilogx(abs(wght),signal,'.');
hold on;
plot(weff,emean,'r-');
plot(wmean,emean+evar,'g-');
plot(wmean,emean-evar,'g-');
plot(wmean,ebrd(:,1),'k--',wmean,ebrd(:,2),'k--');
figure;plot(wmean,evar,'.');
figure;plot(weff.*n(1:end-1),emean,'.')


    
    