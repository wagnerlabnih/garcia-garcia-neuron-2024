function [out,h] = histogram_cell(q,edges,o)
arguments
    q
    edges = [];
    o.cdf = true;
    o.color = colororder;
    o.linestyle = "-";
    o.pos = [];
    o.fig = figure;
    o.kruskwall=false;
    o.title = "";
    o.lw = 1;
    o.medians=false;
end
if ~iscell(q), q = {q}; end
if size(o.color,1)==1, o.color = repmat(o.color,length(q),1); end
if length(o.linestyle)==1, o.linestyle = repmat(o.linestyle,length(q),1); end
if length(o.lw)==1, o.lw = repmat(o.lw,length(q),1); end
if ~isempty(o.pos), set(o.fig,"position",o.pos); end
q = cellfun(@(x)x(:),q,unif=false);
if nargin<2 || isempty(edges)
if size(q{1},1) >1
mn = prctile(cat(1,q{:}),0.5);
mx = prctile(cat(1,q{:}),99.5);
edges = linspace(mn,mx,1000);
else
[~,edges] = histcounts(cat(2,q{:}));
end
end
for g = 1:length(q)
    if o.cdf
    h(g)=histogram(q{g},edges,norm="cdf",displayst="stairs",...
        edgecol=o.color(g,:),linestyle=o.linestyle(g),linewidth=o.lw(g));
    ylim([0 1]);
    else
        h(g)=histogram(q{g},edges,norm="probability",facecol=o.color(g,:));
    end
    hold on;
end
if ~o.cdf, ylim tight; end
box off;
if length(q)==2
[~,s]=ranksum_disp(q,meds=o.medians);
elseif length(q)>2 && length(unique(cat(1,q{:})))<5
[~,~,s] = chisq_disp(q);
s=s+newline+mean_stde_disp(q,meds=o.medians);
elseif length(q)>2 && o.kruskwall==false
[~,~,s] = anova_disp(q);
s=s+newline+mean_stde_disp(q,meds=o.medians);
elseif length(q)>2 && o.kruskwall==true
    [~,s] = kruskwall_disp(q);
    s=s+newline+mean_stde_disp(q,meds=o.medians);
else, s=mean_stde_disp(q,meds=o.medians);
end
if o.title~="", s=o.title+newline+s; end
title(s);
xlabel(inputname(1))
out = s;
set(gcf,"Renderer","painters");