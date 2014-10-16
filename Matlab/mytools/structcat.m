function cat=structcat(s)

n=numel(s);
names=fieldnames(s(1));
m=numel(names);
for i=1:m
    tmp=[];
for j=1:n
    tmp=[tmp;s(j).(names{i})];
end
cat.(names{i})=tmp;
end
