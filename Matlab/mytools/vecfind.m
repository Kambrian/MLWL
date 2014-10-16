function i=vecfind(y,x)
%return the indexes of x in y

y=y(:);
yx=repmat(y,1,numel(x));
xy=repmat(x(:)',numel(y),1);
[i,j]=find(yx==xy);