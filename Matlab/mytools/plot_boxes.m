function plot_boxes(x,y,linespec)
%plot_boxes(x,y,linespec)
%plot x and y data into boxes, with lower and upper boundaries set by first
%and second cols in x and y
% x: [n*2]
% y: [n*2]


n=size(x,1);
for i=1:n
    plot(x(i,[1,2,2,1,1]),y(i,[1,1,2,2,1]),linespec);
    hold on;
end