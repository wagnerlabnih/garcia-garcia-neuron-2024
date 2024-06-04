function varargout = errorbar_shadeSEM(x,yall,color,o)

arguments
    x
    yall
    color=[1 0 0];
    o.linewidth = 1;
    o.linespec = '-';
    o.markersize = 10;
    o.dy = [];
end

% function h = errorbar_patch(x,y,dy,color,linewidth)
% plots the line x vs y and a shaded background defined by dx
%       and returns the handles to the line & patch objects
% y + dy(:,1) is used for the upper boundary of the shaded region
% y - dy(:,2) is used for the upper boundary of the shaded region
% if dy is a vector y - dy(:,1) is used for the lower boundary 
% if color has two rows the 1st is for the line & the 2nd for the shading
% if color has one row a lighter version of the line color is used for the
% shading
% example: 
% figure; errorbar_patch(1:100,rand(100,1)+10,+.1,[0 .6 1],2); shg
% - MAS 9/6/2007

if size(x,1)==1; x=x'; end
% y = squeeze(mean(yall,2));
y = squeeze(nanmean(yall,2));
if isempty(o.dy)
o.dy = squeeze(nanstd_e(yall,2));
end
if size(o.dy,2)==1
o.dy = repmat(o.dy,[1 2]);
end

ylow = y-o.dy(:,1); yhigh=y+o.dy(:,2);
grr = find(isnan(ylow));
ylow(grr) = []; yhigh(grr) = []; x2 = x; x2(grr) = [];
coloreb = color;
coloreb(2,1:3)=1-[1-coloreb(1:3)]/4;
h2 = patch([x2;flipud(x2)],[ylow;flipud(yhigh)],coloreb(2,1:3)); hold on; 
h1 = plot(x,y,o.linespec,'linewidth',o.linewidth,Color=color);
xlim tight; figureformat_forsaving;
try
uistack(h1,"top");
end
set(h2,'edgecolor','none');
set(h1,'markerfacecolor',color(1,1:3),'markersize',o.markersize);
varargout{1} = h1;
varargout{2} = h2;