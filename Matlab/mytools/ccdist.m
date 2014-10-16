function csd=ccdist(x,y)
% cosine of celestial distance(in degrees) between two points at x and y
% 
% calculate cos(d) where d is the angular seperation of two celestial points
% with coordinates x=[ra1,dec1],y=[ra2,dec2]
% x,y can be matrixes of the same size or one matrix one two-element vector
%  *********** all angles in units of degrees ****************
if isempty(x)|isempty(y)
    csd=[];
    return;
else
csd=sind(y(:,2)).*sind(x(:,2))+cosd(y(:,2)).*cosd(x(:,2)).*cosd(y(:,1)-x(:,1));
csd(csd>1)=1;
csd(csd<-1)=-1;
end