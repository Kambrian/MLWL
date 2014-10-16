function cen=GAMArand(ngrps,seed,ngrp0)
% generate random group centers for gama
% input: 
%     ngrps(3): number of groups for the three sky patch
%     seed: an integer seed to grow the random numbers, the same
%           seed gives you the same random catalogue; 
%           seed<0 means generate random centers
%           seed>1 returns a random permutation sequence to permute the centers
%     ngrp0: number of real groups, to draw random redshift from
% output:
%     cen{3}(ngrps,3): random [ra,dec] for three sky patches, and random
%     sample to draw redshift from ngrp0.

% G09:[129,141]x[-1,3]
% G12:[174,186]x[-2,2]
% G15:[211.5,223.5]x[-2,2]
if nargin<3
    ngrp0=[0,0,0];
end
cen=cell(3,1);
if seed<0
    seed=-seed;
    raminmax=[129,141;
        174,186;
        211.5,223.5];
%     decminmax=[-1,3;-2,2;-2,2];
    decminmax=[-0.7,3;-2,2;-2,2];
    
%     raminmax(:,1)=raminmax(:,1)+0.5;
%     raminmax(:,2)=raminmax(:,2)-0.5;
%     decminmax(:,1)=decminmax(:,1)+0.5;
%     decminmax(:,2)=decminmax(:,2)-0.5;
    
    s=RandStream('mt19937ar', 'Seed', seed);
    for sky=1:3
        cen{sky}=rand(s,ngrps(sky),3);
        cen{sky}(:,1)=cen{sky}(:,1)*(raminmax(sky,2)-raminmax(sky,1))+raminmax(sky,1);
        cen{sky}(:,2)=asind(cen{sky}(:,2)*(sind(decminmax(sky,2))-sind(decminmax(sky,1)))+sind(decminmax(sky,1)));
        if ngrp0(sky)>0
        cen{sky}(:,3)=randi(s,ngrp0(sky),ngrps(sky),1); %random sampling, to draw redshift.
        end
    end
elseif seed>1
    s=RandStream('mt19937ar', 'Seed', seed);
    for sky=1:3
        cen{sky}=[randperm(s,ngrps(sky))',randperm(s,ngrps(sky))'];
    end
end
