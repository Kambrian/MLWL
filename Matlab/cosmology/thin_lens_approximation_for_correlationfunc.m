%% to investigate the thin lens approximation for the large scale halo-matter cross correlation
% Conclusion: even for the large-scale halo matter cross correlation, thin
% lens approximation is still an excellent approximation!!
rp=logspace(-1,2,10);
y=proj_corr_lin(rp);
y2=proj_corr_lin_LensSource(rp,0.1,0.7);
%%
figure;
plot(rp,y,'g-',rp,y2,'ro');
legend('Thin Lens','Lensing Kernel Weighted');
xlabel('Projected Distance(Mpc/h)');ylabel('Projected Correlation');