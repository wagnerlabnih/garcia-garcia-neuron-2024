%% 
clear; clc;
load("C:\Users\wagnermj\Desktop\data to share\PkC_ephys.mat");
%% Initialize vars
%Pkj metrics list:
%   1    2     3       4           5                6             7
% Epoch; #; SpRate; timecorr; delayslopesign; prewrewactsign; prewrewmag;
%   8       9    10     11        12         13
% ConfPkj; ISI; Ampl; chansub; presence; presenceFrac
strs = ["","expert"];
pcells_pos = cell(ng,1);
pcells_neg = cell(ng,1);
pcells_p_ix = cell(ng,1);
pcells_n_ix = cell(ng,1);
pcells_all = cell(ng,1);
pcells_all_sprate = cell(ng,1);
pcells_all_omit = cell(ng,1);
pcells_all_ix = cell(ng,1);
pcells_all_cs_rew = cell(ng,1);
pcells_all_cs_omit = cell(ng,1);
pcells_all_cs_ix = cell(ng,1);
for g = 1:ng
    disp(strs(g))
    for j = 1:groupn(g)
        curd = dats{g}{j};
algnuse = "rewAlgn";
tr_use = find(curd.mvlen>6 & curd.goodmvdir);
tr_rew = intersect(tr_use,find(curd.rewarded));
tr_omit = intersect(tr_use,find(~curd.rewarded));
twin = winfn(curd.tmpxx,[-2 2]);
twin_s = curd.tmpxx(twin);
fprintf("#good trials: %d\n",length(tr_use));

twinb = winfn(curd.tmpxx,[-0.4 0.4]);

sp_rew = double(curd.(algnuse).sp_Npixl(tr_rew,:,twin));
sp_omit = double(curd.(algnuse).sp_Npixl(tr_omit,:,twin));

spzsc_rew = curd.(algnuse).spRate_zsc(tr_rew,:,twin);
spzsc_rew_travg = squeeze(mean(spzsc_rew,1));

spzsc_omit = curd.(algnuse).spRate_zsc(tr_omit,:,twin);
spzsc_omit_travg = squeeze(mean(spzsc_omit,1));

sprate_rew = curd.(algnuse).spRate_Npixl(tr_rew,:,twin);
sprate_omit = curd.(algnuse).spRate_Npixl(tr_omit,:,twin);

[~,mxix] = max(spzsc_rew_travg,[],2);
[~,srtix] = sort(mxix);

cellsuse = [];
for i = 1:curd.nCells
    cuse = i;
    per_tr_sprate = squeeze(mean(sprate_rew(:,cuse,:),3));
    per_tr_sprate_75 = prctile(per_tr_sprate(per_tr_sprate>5),75);
    pkjix = find(cuse==pairs{g}(:,3) & j==pairs{g}(:,2) & g==pairs{g}(:,1),1);
    if ~isempty(pkjix), ispkj = 1; else, ispkj=0; end
    if curd.matchedcf(i), ismatched=1; else, ismatched= 0; end
    if curd.inpkjlayer(i), curinpkjlayer = 1; else, curinpkjlayer = 0; end
    if curd.goodacorr(i), curacorr = 1; else, curacorr = 0; end
    goodcell = curd.new_isi(cuse)<2;
    if (goodcell && curinpkjlayer) || ispkj || ismatched
    if per_tr_sprate_75>40 || ispkj || ismatched
    cellpres_1st = tr_rew(find(per_tr_sprate>(per_tr_sprate_75/2),1));
    cellpres_last = tr_rew(find(per_tr_sprate>(per_tr_sprate_75/4),1,'last'));
    if ~isempty(cellpres_1st)
        tr_use_cur = cellpres_1st:cellpres_last;
        [~,tr_rew_use,~] = intersect(tr_use_cur,tr_rew);
        [~,tr_omit_use,~] = intersect(tr_use_cur,tr_omit);
    if numel(tr_rew_use)>5
    sprate_cur = mean(mean(curd.(algnuse).spRate_Npixl(tr_use_cur,cuse,:)));
    if sprate_cur>40 || ispkj || ismatched
        zsc_temp = zscore(squeeze(curd.(algnuse).spRate_Npixl(tr_use_cur,cuse,:)),[],'all');
        rate_temp = squeeze(curd.(algnuse).spRate_Npixl(tr_use_cur,cuse,:));
        spmat_temp = double(squeeze(curd.(algnuse).sp_Npixl(tr_use_cur,cuse,:)));
        bin_len = round(length(twin)/8);
        bins = arrayfun(@(x)mean(spmat_temp(:,bin_len*(x-1)+(1:bin_len)),2),1:8,UniformOutput=false);
        [~,taskmod,s] = anova_disp(bins,disp=0);
        spzsc_rew_cur = zsc_temp(tr_rew_use,twin);
        spzsc_omit_cur = zsc_temp(tr_omit_use,twin);
        spmat_rew_cur = spmat_temp(tr_rew_use,twin);
        spmat_omit_cur = spmat_temp(tr_omit_use,twin);
        sprate_rew_cur = rate_temp(tr_rew_use,twin);
        sprate_omit_cur = rate_temp(tr_omit_use,twin);
twinant = winfn(twin_s,[-1 0]);
typ = "Pearson";
sigscur = reshape(squeeze(cat(1,spzsc_rew_cur(:,twinant), spzsc_omit_cur(:,twinant)))',[],1);
twinantrep = repmat(twin_s(twinant),1,length([tr_rew_use; tr_omit_use]))';
[rcur,pcur] = corr(sigscur,twinantrep,type=typ);
slopesign = (rcur>0) - (rcur<0);
pre_rew_win = winfn([-0.1 -0.025],twin_s);
earlydelwin = winfn([-0.95 -0.8],twin_s);
prerewmag = mean(mean(spzsc_rew_cur(:,pre_rew_win),2));
prerewsign = (mean(mean(spzsc_rew_cur(:,pre_rew_win),2))>0) - ...
    (mean(mean(spzsc_rew_cur(:,pre_rew_win),2))<0);

rcur2= rcur;
isiviol = curd.isi_viol(cuse);
new_isi = curd.new_isi(cuse);
ampcut = curd.ampcutoff(cuse);
cellsuse = cat(1,cellsuse,[cuse rcur2]);
lenuse = diff(curd.midpt(tr_use_cur(tr_rew_use([1 end]))))*.005/60;
totlen = diff(curd.midpt([1 end]))*.005/60;
chansub = curd.channum(cuse);
ixlist = [g j sprate_cur rcur2 slopesign prerewsign prerewmag ispkj||ismatched ...
    new_isi ampcut chansub ...
    lenuse lenuse/totlen];
if slopesign>0 && prerewsign>0
    pcells_pos{g} = cat(1,pcells_pos{g},mean(spzsc_rew_cur,1));
    pcells_p_ix{g} = cat(1,pcells_p_ix{g},ixlist);
elseif slopesign<0 && prerewsign<0
    pcells_neg{g} = cat(1,pcells_neg{g},mean(spzsc_rew_cur,1));
    pcells_n_ix{g} = cat(1,pcells_n_ix{g},ixlist);
end
pcells_all{g} = cat(1,pcells_all{g},mean(spzsc_rew_cur,1));
pcells_all_sprate{g} = cat(1,pcells_all_sprate{g},mean(sprate_rew_cur,1));
twindel = winfn(twin_s,[-1 0]);
pcells_all_sprate{g} = cat(1,pcells_all_sprate{g},mean(sprate_rew_cur,1));
pcells_all_omit{g} = cat(1,pcells_all_omit{g},mean(spzsc_omit_cur,1));
if ispkj
    csix = pairs{g}(pkjix,4);
    sprate_cs = mean(mean(curd.(algnuse).spRate_Npixl(tr_use_cur,csix,:)));
    pcells_all_cs_rew{g} = cat(1,pcells_all_cs_rew{g},...
        squeeze(mean(curd.(algnuse).spRate_Npixl(tr_use_cur(tr_rew_use),csix,twin),1))');
    pcells_all_cs_ix{g} = cat(1,pcells_all_cs_ix{g},[sprate_cs curd.new_isi(csix) curd.ampcutoff(csix) curd.channum(csix)]);
    pcells_all_cs_omit{g} = cat(1,pcells_all_cs_omit{g},...
        squeeze(mean(curd.(algnuse).spRate_Npixl(tr_use_cur(tr_omit_use),csix,twin),1))');
end
pcells_all_ix{g} = cat(1,pcells_all_ix{g},ixlist);
    elseif g==1, disp(i);
    end
    elseif g==1, disp(i+", tr:"+numel(tr_rew_use)); 
    end
    elseif g==1, disp(i);
    end
    end
end
end
end
end
pcells_all_ix
%% Fig 7B, S7B
%Pkj metrics list:
%   1    2     3       4           5                6             7
% Epoch; #; SpRate; timecorr; delayslopesign; prewrewactsign; prewrewmag;
%   8       9    10     11        12         13
% ConfPkj; ISI; Ampl; chansub; presence; presenceFrac
figure(pos=[333         377        600         392]);
ax = [];
for g = 1
    if ~isempty(pcells_all{g})
    ixplt = {}; ccplt = [1 0 0; 0 0 1; 0 1 1; 1 0 1; 0 0 0; 0 1 0];
    ixplt{1} = find(pcells_all_ix{g}(:,5)==1 & pcells_all_ix{g}(:,6)==1);
    ixplt{2} = find(pcells_all_ix{g}(:,5)==-1 & pcells_all_ix{g}(:,6)==-1);
    h = [];
    for k = 1:length(ixplt)
    h(k) = errorbar_shadeSEM(twin_s,pcells_all{g}(ixplt{k},:)',ccplt(k,:)); 
    hold on;
    end
    ylim tight; yl = ylim;
    hold on; plot([-1 -1],yl,'k')
    hold on; plot([0 0],yl,'k')
    uistack(h,'top')
    szs = cellfun(@(x)length(x),ixplt);
    pcts = cellfun(@(x)length(x)/size(pcells_all{g},1)*100,ixplt);
    frs = cellfun(@(x)mean(pcells_all_ix{g}(x,3)),ixplt);
    frs = [frs mean(pcells_all_ix{g}(:,3))];
    title(sprintf("%s. npergrp: ",strs(g))+sprintf("%d ",szs)+...
        sprintf("\npos/neg%%: ")+sprintf("%.3g ",pcts)+", fr: "+sprintf("%.3g ",frs))
    xlabel("Time rel to reward (s)")
    ylabel("Spike rate (zsc)")
    box off
    end
end
figureformat_forsaving
%% Fig S7C, D
%Pkj metrics list:
%   1    2     3       4           5                6             7
% Epoch; #; SpRate; timecorr; delayslopesign; prewrewactsign; prewrewmag;
%   8       9    10     11        12         13
% ConfPkj; ISI; Ampl; chansub; presence; presenceFrac
f1=figure(pos=[333         377        400         392]);
ax = [];
for g = 1
    if ~isempty(pcells_all{g})
    figure(f1);
    ixplt = {}; ccplt = [1 0 0; 0 0 1; 0 1 1; 1 0 1; 0 0 0; 0 1 0];
    ixplt{1} = find(pcells_all_ix{g}(:,5)==1 & pcells_all_ix{g}(:,6)==1);
    ixplt{2} = find(pcells_all_ix{g}(:,5)==-1 & pcells_all_ix{g}(:,6)==-1);
    p1 = -pcells_all{g}(ixplt{1},:);
    p2 = pcells_all{g}(ixplt{2},:);
    p12 = [p1; p2];
    h = errorbar_shadeSEM(twin_s,p12',[0 0 1]); 
    hold on;
    p1o = -pcells_all_omit{g}(ixplt{1},:);
    p2o = pcells_all_omit{g}(ixplt{2},:);
    p12o = [p1o; p2o];
    ho = errorbar_shadeSEM(twin_s,p12o',[1 0 0]);
    ylim tight; yl = ylim;
    hold on; plot([-1 -1],yl,'k')
    hold on; plot([0 0],yl,'k')
    uistack(h,'top')
    szs = sum(cellfun(@(x)length(x),ixplt));
    pcts = cellfun(@(x)length(x)/size(pcells_all{g},1)*100,ixplt);
    frs = cellfun(@(x)mean(pcells_all_ix{g}(x,3)),ixplt);
    frs = [frs mean(pcells_all_ix{g}(:,3))];
    title(sprintf("%s. nposneg: ",strs(g))+sprintf("%d ",szs)+...
        sprintf("\n%%: ")+sprintf("%.3g ",pcts)+", fr: "+sprintf("%.3g ",frs))
    xlabel("Time rel to reward (s)")
    ylabel("Spike rate (zsc)")
    box off
    twin_comp = winfn(twin_s,[0.15 0.35]);
    histogram_cell({mean(p12(:,twin_comp),2),mean(p12o(:,twin_comp),2)},color=[0 0 1; 1 0 0])
    signrank_disp(mean(p12(:,twin_comp),2),mean(p12o(:,twin_comp),2))
    end
end
%% Fig S6H-J
%Pkj metrics list:
%   1    2     3       4           5                6             7
% Epoch; #; SpRate; timecorr; delayslopesign; prewrewactsign; prewrewmag;
%   8       9    10     11        12         13
% ConfPkj; ISI; Ampl; chansub; presence; presenceFrac
ixuse = 7;
g=1;

ixpkj = find(pcells_all_ix{g}(:,8));
ixrest = find(~pcells_all_ix{g}(:,8));

barwitherr_stde({pcells_all_ix{g}(ixpkj,3),pcells_all_ix{g}(ixrest,3)},dots=true); 
barwitherr_stde({pcells_all_ix{g}(ixpkj,9),pcells_all_ix{g}(ixrest,9)},dots=true); 
histogram_cell({pcells_all_ix{g}(ixpkj,10),pcells_all_ix{g}(ixrest,10)});
%% Fig S7A, S7E
%Pkj metrics list:
%   1    2     3       4           5                6             7
% Epoch; #; SpRate; timecorr; delayslopesign; prewrewactsign; prewrewmag;
%   8       9    10     11        12         13
% ConfPkj; ISI; Ampl; chansub; presence; presenceFrac
g=1;
for k_run = 1:2
ixuse = 7;
% 
if k_run==1
ixc_use = 1:size(pcells_all_ix{g},1);
ylabstr = "All PkCs";
else
ixc_use = find(pcells_all_ix{g}(:,8));
ylabstr = "Confirmed PkCs";
end
[~,ixsrt] = sort(pcells_all_ix{g}(ixc_use,ixuse));
ixsrt = ixc_use(ixsrt);
figure;
subaxis(1,2,1);
imagesc(twin_s,[],pcells_all{g}(ixsrt,:),[-0.4 0.6]);
ylim tight; yl = ylim;
    hold on; plot([-1 -1],yl,'k')
    hold on; plot([0 0],yl,'k')
xlabel("Time rel to reward (s)")
ylabel(ylabstr)
title("Rewarded")
subaxis(1,2,2);
imagesc(twin_s,[],pcells_all_omit{g}(ixsrt,:),[-0.4 0.6]);
ylim tight; yl = ylim;
    hold on; plot([-1 -1],yl,'k')
    hold on; plot([0 0],yl,'k')
title("Omitted")
colorbar_noresize;
end
%% Fig S6G
figure;
g=1;
tmp = pcells_all_cs_rew{g}';
h1=errorbar_shadeSEM(twin_s,tmp,[0 0 1]); 
xlim([-1.95 2]);
ylim tight; yl = ylim;
    hold on; plot([-1 -1],yl,'k')
    hold on; plot([0 0],yl,'k')
uistack([h1],"top")
figureformat_forsaving;
win_rew = winfn([0.09,0.12],twin_s);
ylabel("Complex spike rate (Hz)");
xlabel("Time relative to reward (s)")