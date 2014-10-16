function cat=gama_trim(cat,r)

names=fieldnames(cat{1});
raminmax=[129,141;
        174,186;
        211.5,223.5];
decminmax=[-1,3;
        -2,2;
        -2,2];
    for sky=1:3
    f=(cat{sky}.BCGRA<raminmax(sky,1)+r)|(cat{sky}.BCGRA>raminmax(sky,2)-r)|(cat{sky}.BCGDEC<decminmax(sky,1)+r)|(cat{sky}.BCGDEC>raminmax(sky,2)-r);
    for i=1:numel(names)
    cat{sky}.(names{i})(f)=[];
    end
    end