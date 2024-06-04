function out = barwitherr_stde(grps,o)
arguments 
    grps
    o.dots = [];
    o.fsize = [1138          22         max(75*(length(grps)),210)         420];
    o.fhandle = [];
    o.color = [0 0 1];
    o.xcoords = 1:length(grps);
    o.msz = 0.1;
    o.dotfc = [0.7 0.7 0.7];
    o.dotec = [0 0 0];
    o.dotcalc="90";
    o.dotgrp=[];
    o.kruskwall=false;
    o.signrank=false;
    o.ranksumgrps=[];
    o.anovagrp = [];
    o.altplt="no";
    o.overlap=1;
    o.szlims=[0 inf];
    o.medians = false;
    o.title="";
end
if ~iscell(grps), grps = num2cell(grps,1); end
grps = cellfun(@(x)x(:),grps,unif=false);
means = cellfun(@(x)mean(double(x)),grps);
ses = cellfun(@(x)std_e(double(x),[],true),grps);
if ~isempty(o.fhandle), figure(o.fhandle); hold on; 
else, figure(pos=o.fsize)
end
if size(o.color,1)~=length(means), o.color = repmat(o.color(1,:),length(means),1); end
out = []; herrs = [];
if o.altplt=="no"
for i = 1:length(means)
[out(i),herrs(i)] = barwitherr(ses(i),o.xcoords(i),means(i),edgecol="none");
hold on;
set(out(i),"FaceColor",o.color(i,:));
end
% disp(o.xcoords); pause;
xticks(sort(unique(o.xcoords)));
% for i = 1:length(means)
%     out.CData(i,:) = o.facecol(i,:);
% end
box off;
if isempty(o.fhandle)
if isempty(o.xcoords)
xlim([0.5 length(grps)+0.5])
else
    xlim([min(o.xcoords)-0.5 max(o.xcoords)+0.5])
end
end
elseif o.altplt=="box"
    boxplotcell(grps);
elseif o.altplt=="violin"
    violin(grps,'edgecolor',o.color,'facecolor','none','medc',[]);
end
if ~isempty(o.dots)
    min_n = min(cellfun(@(x)length(x),grps));
    if (islogical(o.dots) && o.dots) || (isstring(o.dots) && o.dots=="all")
        grps_use = grps;
    elseif (isstring(o.dots) && o.dots=="equal")
        grps_use = cellfun(@(x)shufflev(x,1,min_n),grps,unif=false);
    elseif isnumeric(o.dots)
        min_n = min(min_n,o.dots);
        grps_use = cellfun(@(x)shufflev(x,1,min_n),grps,unif=false);
    end
    [xspread_all,mszn] = spread_dots(o.xcoords,grps_use,dotsize=o.msz,...
        dotcalc=o.dotcalc,overlap=o.overlap,szlims=o.szlims);
    for i = 1:length(means)
        if isempty(o.dotgrp) || ismember(i,o.dotgrp)
    plot(xspread_all{i},grps_use{i},'o',markersize=mszn,...,
        markeredgecolor=o.dotec, markerfacecolor=o.dotfc);
        end
    end
end
ylim tight;
if length(grps)==2
[~,s]=ranksum_disp(grps);
elseif length(grps)>2 && length(unique(cat(1,grps{:})))<5
[~,~,s] = chisq_disp(grps);
s=s+newline+mean_stde_disp(grps,meds=o.medians);
elseif length(grps)>2 && o.kruskwall==false
    if isempty(o.anovagrp)
        [~,~,s] = anova_disp(grps);
    else
        [~,~,s] = anova_disp(grps(o.anovagrp));
    end
s=s+newline+mean_stde_disp(grps,meds=o.medians);
elseif length(grps)>2 && o.kruskwall==true
    [~,s] = kruskwall_disp(grps);
    s=s+newline+mean_stde_disp(grps,meds=o.medians);
else
    s=mean_stde_disp(grps,meds=o.medians);
end
if o.signrank
    s= s+newline;
    for i = 1:length(grps)
        s=s+sprintf("p%d=%1.1g, ",i,signrank(grps{i}));
    end
end
if ~isempty(o.ranksumgrps)
    s=s+newline;
for i = 1:length(o.ranksumgrps)
    cg = o.ranksumgrps{i};
    s=s+sprintf("p-"+strvec(cg)+"=%1.1g, ",ranksum(grps{cg(1)},grps{cg(2)}));
end
end
if o.title~="", s=o.title+newline+s; end
title(s);
ylabel(inputname(1),interp="none");
uistack(herrs,"top");
figureformat_forsaving;
if nargout==0, clear out; end