function [h,s]=boxplotcell(varargin)

c = varargin{1};
assert(iscell(c));
varargin{1} = {};

catv = [];
catg = [];
for k = 1:length(c)
    if size(c{k},2)~=1, c{k} = c{k}'; end
    catv = cat(1,catv,c{k});
    catg = cat(1,catg,k*ones(size(c{k})));
end
figure(pos=[1138          22         75*(length(c))         420])
% h=boxplot(catv,catg,notch=true);
h=boxplot(catv,catg,notch=true,whisk=inf);
box off;
% h=boxplot(catv,catg,'whisker',100);
if length(c)==2
[~,s]=ranksum_disp(c);
title(s);
elseif length(c)>2
[~,~,s] = anova_disp(c);
title(s);
end
ylabel(inputname(1));
set(gcf,"Renderer","painters");