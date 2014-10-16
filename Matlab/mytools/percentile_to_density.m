function levels=percentile_to_density(counts, percents)
% levels=percentile_to_density(counts, percents)
% to convert a containment percent to density contour level

n=sort(counts(:),'descend');
frac=cumsum(n)/sum(n);
[m,ind]=unique(n,'last','legacy');
m(1)=[];
ind(1)=[];
levels=spline(frac(ind),m,percents);