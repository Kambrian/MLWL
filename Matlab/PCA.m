f=grpmock.Mult>2&grpmock.MassProxy>1e6;%&~isnan(grpmock.IterCenSM)&grpmock.IterCenSM>1e8&grpmock.IterCenSM<1e12;
X=[grpmock.Mult(f),1+grpmock.Zfof(f),grpmock.VelDisp(f),grpmock.Rad50(f),grpmock.TotFluxProxy(f)];%,grpmock.IterCenSM(f)];
y=mockmass.MIter(f);
X=log(X);
y=log(y);
%%
f=grp.Mult>2&grp.MassProxy>1e6&~isnan(grp.IterCenSM)&grp.IterCenSM>1e8&grp.IterCenSM<1e12;
X=[grp.Mult(f),1+grp.Zfof(f),grp.VelDisp(f),grp.Rad50(f),grp.TotFluxProxy(f),grp.IterCenSM(f)];
X=log(X);
coeffit=[0.65,0,1.15,-0.4,0,1.1]';
y=X*coeffit+14.26*log(10); %best fit WL mass
%%
[pc,score,latent,tsquare]=princomp((X));
c=X*pc;
%%
yref=X(:,1:5)*refcoef;
figure;
plot(yref,c,'.')
%%
figure;
plot(y,c,'.');
%%
[LoadingsPM,specVarPM,TPM,stats,F] = factoran(X, 2,'rotate','promax');
%%
figure;
% plot(y,F,'.')
plot(y,X*LoadingsPM,'.');
%%
figure
plot(y,X,'.')
%%
XX=[X,ones(size(X,1),1)];
XX\y