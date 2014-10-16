files=dir('*.fig');
for i=1:numel(files)
    f=files(i).name;
    openfig(f);print('-depsc',[f(1:end-4),'.eps']);
end