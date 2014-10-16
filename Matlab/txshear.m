function [et,ex,r,f]=txshear(x,y,e,rmax)
% rotate ellipticity e to give tangential and cross shear for lens at x and
% source at y in (ra,dec) coordinates
% input ellipticity e in Mandelbaum convention: in (ra,dec) coordinates,
% e1+ means N-S (or x) direction (which is different from usual definition) and
% e2+ means NE-SW direction (this is the same as usual definition).
% by-product: angular seperation r
%             filter array f, used exclude r>rmax
%  *********** all angles in units of degrees ****************

cosls=ccdist(x,y);
sinls=sqrt(1-cosls.^2);
r=asind(sinls);

f=r>rmax;
y(f,:)=[];
e(f,:)=[];
r(f)=[];
cosls(f,:)=[];
sinls(f,:)=[];

cosf=cosd(x(2)).*sind(y(:,1)-x(1))./sinls;
sinf=(cosls.*sind(y(:,2))-sind(x(2)))./cosd(y(:,2))./sinls;

sin2f=2*sinf.*cosf;
cos2f=cosf.^2-sinf.^2;

e(:,1)=-e(:,1); %change to usual convention
et=-e(:,1).*cos2f-e(:,2).*sin2f;
ex=e(:,1).*sin2f-e(:,2).*cos2f;



