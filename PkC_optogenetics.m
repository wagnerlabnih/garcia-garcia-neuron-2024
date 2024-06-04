%%
clear; clc;
load("C:\Users\wagnermj\Desktop\data to share\PkC_optogenetics.mat")
%% Fig 2I, 2J, S2J, S2L, S2K, 2M
for kops = [1:4 8]
    if kops<5, inc=0; else, inc=2; end
    clear h;
switch kops
    case 1
        groupstoplt = {[1 3],[2 4]};
    case {2,3,4}
        groupstoplt = {[2 4]};
    case 8
        groupstoplt = {2};
end
nplts = length(groupstoplt);
[varsToPlt,xlms] = deal(cell(1,nplts));
varsToPlt(1:end) = {lick_alltrs_rate_filt{kops}};
leglbls = cell(1,nplts);
ttls = ["lick rate (hz)","lick rate (hz)","forelimb y"];
hlfmxa={};
for i = 1:nplts
h = []; leglabcur = [];
curvar = varsToPlt{i};
figure(pos=[564    82   380   420]);
h(end+1)=errorbar_shadeSEM(curd.tmpxx(xplt{kops}),curvar{groupstoplt{i}(1)}',[0 0 0]);
if kops~=8, legh = h(1);leglabcur = trgrplabls(groupstoplt{i}(1)); end
titstr = sprintf("%s\nn-"+trgrplabls{groupstoplt{i}(1)}+"=%d",ttls(i),size(curvar{groupstoplt{i}(1)},1));
hold on;
if length(groupstoplt{i})>1, 
    h(end+1)=errorbar_shadeSEM(curd.tmpxx(xplt{kops}),curvar{groupstoplt{i}(2)}',[0 0.678 0.933]); 
    titstr = titstr+sprintf("\nn-"+trgrplabls{groupstoplt{i}(2)}+"=%d",size(curvar{groupstoplt{i}(2)},1));
    legh = [legh h(2)];
    leglabcur = [leglabcur trgrplabls(groupstoplt{i}(2))];
end
if kops==8
    tg = groupstoplt{i}(1); tstr = ["D5-6","D3-4","D1-2"]; 
    tstr = [tstr "\nRecovery exp"];
    cc = repmat([1 1 1],3,1).*[.2; .4; .7]; cc = [cc; 0 .5 .5];
    kkuse=[7:-1:5 9]'; 
    lk=length(kkuse);
    for kuse = 1:lk
    kk = kkuse(kuse);
    h(end+1)=errorbar_shadeSEM(curd.tmpxx(xplt{kops}),lick_alltrs_rate_filt{kk}{tg}',cc(kuse,:));
    titstr = titstr + sprintf(", n-"+tstr(kuse)+"=%d",size(lick_alltrs_rate_filt{kk}{tg},1));
    end
    legh = h; 
    leglabcur = dsetlabls([8:-1:5 9]); 
end
ylabel(ttls(i))
if ~isempty(xlms{i}), xlim(xlms{i}); end
ylim tight;
yl2 = ylim; yl2(1) = 0; ylim(yl2);
if i==1, yl = ylim; 
elseif i<nplts, ylim(yl); 
end
hold on; plot([-1 -1]+inc,ylim,'k')
hold on; plot([0 0]+inc,ylim,'k')
hp = []; 
hold on; 
firston = arrayfun(@(x)find(laser_trs{kops}{6}(x,:)>0,1),1:size(laser_trs{kops}{6},1));
if ~isempty(firston)
laser_on_mean = curd.tmpxx(xplt{kops}(round(median(firston))));
laston = arrayfun(@(x)find(laser_trs{kops}{6}(x,:)>0,1,'last'),1:size(laser_trs{kops}{6},1));
laser_off_mean = curd.tmpxx(xplt{kops}(round(median(laston))+1));
if i==nplts, inc2=1; else, inc2=0; end
hp=patch([laser_on_mean laser_on_mean laser_off_mean laser_off_mean]+inc2,[yl yl(2) yl(1)],[1 0.9 0.9],edgecolor="none"); 
laser_mean_all(kops,:) = [laser_on_mean laser_off_mean ]
titstr = titstr+newline+"laser: "+strvec([laser_on_mean laser_off_mean])+"s";
end
uistack(h,"top");
uistack(hp,"bottom")
titstr = dsetlabls(kops)+": "+titstr;
legend(legh,leglabcur);
title(titstr);
figureformat_forsaving;
xlabel("Time rel to reward (s)")
% pause; 
end
end
%% Stats/n
for kops = 1:length(laser_trs)
    ntr_on = size(laser_trs{kops}{6},1);
    ntr_off = size(laser_trs{kops}{5},1);
    ntr_tot = ntr_on+ntr_off;
    optofrac = round(ntr_on/ntr_tot*100);
    nooptofrac = round(ntr_off/ntr_tot*100);
    nrew = size(laser_trs{kops}{1},1)+size(laser_trs{kops}{3},1);
    optorewfrac = round(size(laser_trs{kops}{3},1)/nrew*100);
    nomit = size(laser_trs{kops}{2},1)+size(laser_trs{kops}{4},1);
    optoomitfrac = round(size(laser_trs{kops}{4},1)/nomit*100);
    npulses = sum(diff(laser_trs{kops}{6},[],2)>0,2);
    firston = arrayfun(@(x)find(laser_trs{kops}{6}(x,:)>0,1),1:size(laser_trs{kops}{6},1));
    if ~isempty(firston)
    laser_on_mean = curd.tmpxx(xplt{kops}(round(firston)));
    laston = arrayfun(@(x)find(laser_trs{kops}{6}(x,:)>0,1,'last'),1:size(laser_trs{kops}{6},1));
    laser_off_mean = curd.tmpxx(xplt{kops}(round(laston)+1));
    end
    disp(dsetlabls(kops)+": optofrac:"+optofrac+", optorewfrac:"+optorewfrac+...
        ", optoomitfrac:"+optoomitfrac+", laseron:"+mean_stde_disp(laser_on_mean)+...
        ", laseroff:"+mean_stde_disp(laser_off_mean)+", npulses:"+mean_stde_disp(npulses))
end
%% Fig 2K Main opto quantification
ccplt = [0 0 0; 0 0.678 0.933];
kops = 1;
winearly = winfn(curd.tmpxx(xplt{kops}),[0.35 0.8]);
winlate = winfn(curd.tmpxx(xplt{kops}),[1.7 2.1]);
qchange = laseronminusoff;
rewlaserchange_early = mean(qchange{1}{1}(:,winearly),2);
rewlaserchange_late = mean(qchange{1}{1}(:,winlate),2);
omitlaserchange_early = mean(qchange{1}{2}(:,winearly),2);
omitlaserchange_late = mean(qchange{1}{2}(:,winlate),2);
omitlaserchange_early_short = mean(qchange{2}{2}(:,winearly),2);
omitlaserchange_late_short = mean(qchange{2}{2}(:,winlate),2);
omitcontrolchange_early = mean(qchange{4}{2}(:,winearly),2);
omitcontrolchange_late = mean(qchange{4}{2}(:,winlate),2);
omitLIXchange_early = mean(qchange{3}{2}(:,winearly),2);
omitLIXchange_late = mean(qchange{3}{2}(:,winlate),2);
%
barwitherr_stde({rewlaserchange_early,rewlaserchange_late,omitlaserchange_early,...
    omitlaserchange_late,omitlaserchange_early_short,omitlaserchange_late_short,...
    omitcontrolchange_early, omitcontrolchange_late,omitLIXchange_early,omitLIXchange_late},...
    color=[0 0 1; 0 0 1; 1 0 0; 1 0 0; 1 .4 .4; 1 .4 .4;...
    0 0 0; 0 0 0; 0 0 0; 0 0 0],...
    ranksumgrps={[1,2],[3,4],[5,6],[7,8],[9,10]},...
    signrank=true,dots=true,dotcalc="90");
ylabel("Laser ΔLick rate (Hz)")
%% Fig 2N, 2-s RELEARN QUANTIFICATION
wino_s = [2 2.75];
winomit = winfn(curd.tmpxx(xplt{9}),wino_s);
kuse = [5:9];
winbase_s = [0.3 0.9];
basewin = winfn(curd.tmpxx(xplt{9}),winbase_s);
omitlick = {};
omitquant = {};
nku = length(kuse);
for kk = 1:nku
    k=kuse(kk);
    oix = 2;
    if k<9, rix = 3; else, rix=1; end
omitlick{kk} = lick_alltrs_rate_filt{k}{oix};
rewlick{kk} = lick_alltrs_rate_filt{k}{rix};
basequant = mean(omitlick{kk}(:,basewin),2);
basequantrew = mean(rewlick{kk}(:,basewin),2);
scur = unique(s_ix{k}{oix});
for s = scur'
    scurix_o = find(s_ix{k}{oix}==s); 
    scurix_r = find(s_ix{k}{rix}==s); 
    basequant(scurix_o) = mean(basequantrew(scurix_r));
end
omitquant{1}{kk} = mean(omitlick{kk}(:,winomit),2);
omitquant{2}{kk} = mean(omitlick{kk}(:,winomit),2)-basequant;
end
ccplt = [0.75 0.75 0.75; ...
    0.6 0.6 0.6; 0.3 0.3 0.3; 0 0 0];
ccplt = [ccplt; 0 .5 .5];
ls = ["-" "-" "-" "-" "-"];
[~,h]=histogram_cell(omitquant{2},cdf=true,color=ccplt,linestyle=ls);
legend(h,dsetlabls(kuse))
xlabel("Omission lick rate "+strvec(wino_s)+"s minus "+strvec(winbase_s)+"s")
hold on; plot([0 0],ylim,'k')

ranksum_disp(omitquant{2}([4 5]))
ranksum_disp(omitquant{2}([3 5]))
ranksum_disp(omitquant{2}([2 5]))
ranksum_disp(omitquant{2}([1 5]))
[~,~,s]=anova_disp(omitquant{2}(1:4))
%% Fig S2M, single trials
for kops=1:4
switch kops
    case 1, ixsrt=27;
    case 2, ixsrt=7;
    case 3, ixsrt=20;
    case 4, ixsrt=21;
end
figure
    plot(curd.tmpxx(xplt{1}),lick_alltrs_raw{kops}{4}(ixsrt,:));
    hold on; 
    plot(curd.tmpxx(xplt{1}),laser_trs{kops}{4}(ixsrt,:)-1,'r');
    hold on; plot([0 0],ylim,'k'); plot([-1 -1],ylim,'k');
    figureformat_forsaving; yticks([]);
    xlim tight;
    title(ixsrt);
end
