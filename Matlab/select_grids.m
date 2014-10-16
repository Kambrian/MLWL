function ilist=select_grids(xc,rmax,xmin,xstep,n,grids)
%select grids intervening with rmax centered at xc
%grids are specified with xmin,xstep and index
%return all the member in the selected grids as a single list
%**** all angles in units of degrees ****

rm=[asind(sind(rmax)/cosd(xc(2))),rmax];
if ~isreal(rm(1)), rm(1)=180; end
imin=(xc-rm-xmin)./xstep;
imin=ceil(imin);
imax=(xc+rm-xmin)./xstep;
imax=floor(imax)+1;
imin(imin<1)=1;
imax(imax>n)=n(imax>n);
ilist=grids(imin(1):imax(1),imin(2):imax(2));
ilist=cell2mat(ilist(:));

        