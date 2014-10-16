function [theta,fai]=sphshift2(ra,dec,ra0,dec0)
% shift the coordinates (ra,dec)  to a non-rotated local reference frame
% centered at (ra0,dec0), with e_fai aligned to e_dec locallly 
% this function shifts the reference frame's origin to (ra0,dec0), so that
%    sphshift(ra0,dec0,ra0,dec0)=[0,0]
%i.e., it finds the new (local) coordinate of point (ra,dec) in a new
%reference frame centered at (ra0,dec0), and have e_fai aligned with e_dec
%locally
%
%the transformation is a rotation of the global axes around z with ra0, 
%followed by a rotation around y with -dec0. Note that axes rotation is
%reverse to object rotation
% units: degrees

Ry=[cosd(dec0),0,sind(dec0);
    0 1 0;
    sind(-dec0),0,cosd(dec0)];
Rz=[cosd(ra0),sind(ra0),0;
       sind(-ra0),cosd(ra0),0;
       0 0 1];

dim=size(ra);
ra=ra(:)*pi/180;
dec=dec(:)*pi/180;   
[x,y,z]=sph2cart(ra,dec,1);
M=[x,y,z]*Rz'*Ry';
[theta,fai]=cart2sph(M(:,1),M(:,2),M(:,3));
%theta(theta<0)=theta(theta<0)+2*pi;
theta=theta*180/pi;
fai=fai*180/pi;
theta=reshape(theta,dim);
fai=reshape(fai,dim);

