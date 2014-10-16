function [w,r]=pair_decomposition(x,y)
% lens at x and source at y in (ra,dec) coordinates
% et=e[0]*w[0]+e[1]*w[1], in Mandelbaum convention, where chi[0] differs by
% a sign from the usual convention.
% return Lens-Source angular distance r and w[2].
% y can be vector of size n, then w is of shape [n,2]
%  *********** all angles in units of degrees ****************

cosls=ccdist(x,y);
sinls=sqrt(1-cosls.^2);
r=asind(sinls);

cosf=cosd(x(2)).*sind(y(:,1)-x(1))./sinls;
sinf=(cosls.*sind(y(:,2))-sind(x(2)))./cosd(y(:,2))./sinls;

sin2f=2*sinf.*cosf;
cos2f=cosf.^2-sinf.^2;

w=[cos2f,-sin2f];




