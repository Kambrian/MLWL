function h=myfigure(args)
% h=myfigure(args)
% same usage as figure(), but produces publication quality plots
%
%---------------------
% Jiaxin Han, 05/2008
%---------------------

if nargin<1
    h=figure;
else
h=figure(args);
end
set(gcf,...
    'DefaultLineLineWidth',2,'DefaultAxesLineWidth',.5,...
    'DefaultAxesFontName','Helvetica',...
    'DefaultAxesFontSize',20,...
    'DefaultAxesTickLength',[0.02,0.02],... 
    'DefaultAxesXMinorTick','on','DefaultAxesYMinorTick','on');
set(gcf,'DefaultLineMarkerSize',6);
set(gcf,'DefaultTextInterpreter','latex');
set(gcf,'paperunits','centimeters','paperposition',[0.6,6,20,17]);

% Font: Bitstream Vera Sans as in matplotlib
%     'DefaultTextFontName','Bitstream Vera Sans',...
