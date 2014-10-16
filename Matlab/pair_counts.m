function [totalcounts,totalrsum]=pair_counts(data1,data2,rbin)
%count pairs
%data: struct{RA,DEC,z}
%rbin: radial bins
% return: rsum,counts
ZDIFF=0.005;
OmegaM=0.3;
FlagComoving=1;
% global OmegaM FlagComoving
nbin=numel(rbin)-1;
d1=cell(3,1);
d2=cell(3,1);
[d1{1},d1{2},d1{3}]=split_mockcat(data1);
[d2{1},d2{2},d2{3}]=split_mockcat(data2);
% global ll
ll=cell(3,1);
ngrids=[40,20];
for sky=1:3
    disp(['making linklist for sky',num2str(sky),'...']);
    ll{sky}.ngrid=ngrids;
    [ll{sky}.grids,ll{sky}.xrange,ll{sky}.yrange,ll{sky}.step]=linklist(d2{sky}.RA,d2{sky}.DEC,ll{sky}.ngrid);
end
%     eval(['save GAMAlinklist_',num2str(ngrids(1)),'_',num2str(ngrids(2)),'.mat ll'])

totalcounts=zeros(nbin,1);
totalrsum=zeros(nbin,1);
for sky=1:3
    switch FlagComoving  %this only affects the surface density estimator, but not the mass estimator
        case 0
            scale=ones(size(d1{sky}.RA));
        case 1
            scale=1./(1+d1{sky}.z);
        otherwise
            error('flag_comoving must be 0 or 1');
    end
    dl=AD_dist_flat(OmegaM,0,d1{sky}.z);  %angular diameter distance for lens
    dl=dl./scale; %proper dl to use with rbin
    prog=0;
    fprintf(1,'sky %02d:    ',sky);
    ncen=numel(d1{sky}.RA);
    for cenid=1:ncen
        if cenid>ncen/10*prog, fprintf(1,['\b\b\b',num2str(prog*10,'%02d'),'%%']);prog=prog+1; end
        xcen=[d1{sky}.RA(cenid),d1{sky}.DEC(cenid)];
        tbin=rbin./dl(cenid)*180/pi;
        thetamax=tbin(end);
        if thetamax>90, thetamax=90;end % we do not need to search half of the sky
        sources=select_grids(xcen,thetamax,[ll{sky}.xrange(1),ll{sky}.yrange(1)],ll{sky}.step,ll{sky}.ngrid,ll{sky}.grids);
        sources(abs(d2{sky}.z(sources)-d1{sky}.z(cenid))>ZDIFF)=[];
        cosls=ccdist(xcen,[d2{sky}.RA(sources),d2{sky}.DEC(sources)]);
        sinls=sqrt(1-cosls.^2);
        theta=asind(sinls);
        if isempty(theta) continue;end
        [counts,bin]=histc(theta,tbin);
        counts(end)=[];
        if numel(theta)<=1, counts=reshape(counts,nbin,1);end
        if nargout>1
            rsum=zeros(nbin,1);
            for i=1:nbin
                if counts(i)
                    rsum(i)=sum(theta(bin==i));
                end
            end
            rsum=rsum/180*pi*dl(cenid);
            totalrsum=totalrsum+rsum;
        end
        totalcounts=totalcounts+counts;
    end
end
