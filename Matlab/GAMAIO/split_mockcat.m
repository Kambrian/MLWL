function [cat09,cat12,cat15]=split_mockcat(catmock)
% to split gama mock catalogue to three patches corresponding to three GAMA
% also applicable to grpcat and galcat

% n=numel(catmock);
if isfield(catmock,'CenRA')
    x=catmock.CenRA;
elseif isfield(catmock,'RA')
        x=catmock.RA;
else
    error('wrong catalogue type to split');
end

names=fieldnames(catmock(1));
m=numel(names);
f1=x>120&x<150;
f2=x>150&x<200;
f3=x>200&x<240;
for i=1:m
    cat09.(names{i})=catmock.(names{i})(f1);
    cat12.(names{i})=catmock.(names{i})(f2);
    cat15.(names{i})=catmock.(names{i})(f3);
end

