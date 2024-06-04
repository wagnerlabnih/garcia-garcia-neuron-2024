function [xout_all,dotsize_norm] = spread_dots(all_x,all_y,o)
arguments
    all_x
    all_y
    o.dotsize = 0.1;
    o.dotcalc = "min";
    o.overlap = 1;
    o.szlims = [0 inf];
end
if ~iscell(all_x)
    assert(length(all_x)==length(all_y));
    all_x_o = all_x;
    all_x = cell(size(all_y));
    for i = 1:length(all_y)
        all_x{i} = all_x_o(i)*ones(size(all_y{i}));
    end
end
xout_all = cell(size(all_y));
ngtot = length(all_x);
dotsize_norm_all = zeros(length(all_y),1);
for k = 1:ngtot
all_y{k} = all_y{k}(:);
all_x{k} = all_x{k}(:);
y = all_y{k};
x = all_x{k};
y = y(:);
x = x(:);
if length(x)==1, x = repmat(x,size(y,1),1); end
ymin = min(y); ymax = max(y);
yl = ylim;
ymin = min(ymin,yl(1));
ymax = max(ymax,yl(2));
rn = range(y);
ax = gca();
xl = xlim;
rnx = range(xl);
set(ax,'Units','points');
ax_pos = tightPosition(ax);
xw_pts = ax_pos(3);
% o.dotsize_xax = o.dotsize*xw_pts/rnx;
o.dotsize_yax = rn*o.dotsize;

groupvar = ymin:o.dotsize_yax:ymax;
groupvar = linspace(ymin,ymax,length(groupvar));
binsz = mean(diff(groupvar));
ng = length(groupvar)-1;
[gcounts,~,gbins] = histcounts(y,groupvar);
dotsize_norm_all(k) = min(o.dotsize*xw_pts/rnx,...
    1/max(gcounts)*xw_pts/rnx)*o.overlap;
gbins_all{k} = gbins;
gcounts_all{k} = gcounts;
end
if o.dotcalc=="min"
dotsize_norm = min(dotsize_norm_all);
elseif o.dotcalc=="90"
    dotsize_norm = prctile(dotsize_norm_all,90);
end
dotsize_norm = max(1,dotsize_norm);
dotsize_norm = coerce(dotsize_norm,o.szlims);
for k = 1:ngtot
gbins = gbins_all{k};
gcounts = gcounts_all{k};
y = all_y{k};
x = all_x{k};
ng = length(gcounts);
xout = zeros(size(x));
for g= 1:ng
    cur_n=gcounts(g);
    if cur_n>0
        curix = find(gbins==g);
        jitrng = 0.4-dotsize_norm*rnx/xw_pts;
        jitmin = max(-jitrng,-o.dotsize*round(cur_n/2));
        jitmax = min(jitrng,o.dotsize*round(cur_n/2));
        jitters = (linspace(jitmin,jitmax,cur_n))';
        if cur_n==1, jitters=0; end
        % plot(x(curix)+jitters,y(curix),'.k',...
        %     markersize=dotsize_norm);
        xout(curix) = x(curix)+jitters;
    end
end
xout_all{k} = xout;
end