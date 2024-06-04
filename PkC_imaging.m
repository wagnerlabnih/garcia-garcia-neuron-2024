%%
clear; clc
load("C:\Users\wagnermj\Desktop\data to share\PkC_imaging.mat");
nep = 2;
nd = cellfun(@(x)length(x),dats);
%% Initialize
travg_rew_neg_all = cell(1,nep);
travg_omit_neg_all = cell(1,nep);
travg_rew_all = cell(1,nep);
trlickavg_rew_all = cell(1,nep);
corrs_travg_all = cell(1,nep);
corrs_travg_all_1s = cell(1,nep);
corrs_travg_all_2s = cell(1,nep);
travg_rew_pos_all = cell(1,nep);
prerewmag_all = cell(1,nep);
for ep = 1:nep
for j = 1:nd(ep)
    curd = dats{ep}{j};

curcs = 1:curd.nIC_GrC;
twin = winfn(curd.tmpxCb,[-1 3]);
twin_s = curd.tmpxCb(twin);
twin_b = winfn(curd.tmpxx,[-1 3]);
twin_b_s = curd.tmpxx(twin_b);

tr_use = curd.mvlen>6 & curd.goodmvdir;
tr_rew = find(tr_use & curd.rewarded);

tr_omit = find(tr_use & ~curd.rewarded);% tr_omit(end)=[];
cursigs_o = curd.midAlgn.sigFilt_GrC;
[b,a] = butter(2,2/15);
cursigs = cursigs_o;
travg_rew = squeeze(mean(cursigs(tr_rew,curcs,twin),1));
travg_omit = squeeze(mean(cursigs(tr_omit,curcs,twin),1));
ix = 1:length(curcs);
if ep==2
twin_prerew = winfn(twin_s,[1.7 2]);
twinant = winfn(twin_s,[0 2]);
else
    twin_prerew = winfn(twin_s,[0.7 1]);
    twinant = winfn(twin_s,[0 1]);
end
twin_1s = winfn(twin_s,[0 1]);
twin_2s = winfn(twin_s,[0 2]);
corrs_travg = corr(travg_rew(:,twinant)',twin_s(twinant)');
corrs_travg_all{ep} = cat(1,corrs_travg_all{ep},corrs_travg);
negcs = curcs(corrs_travg<0);

corrs_travg_1s = corr(travg_rew(:,twin_1s)',twin_s(twin_1s)');
corrs_travg_all_1s{ep} = cat(1,corrs_travg_all_1s{ep},corrs_travg_1s);
corrs_travg_2s = corr(travg_rew(:,twin_2s)',twin_s(twin_2s)');
corrs_travg_all_2s{ep} = cat(1,corrs_travg_all_2s{ep},corrs_travg_2s);

prerewmag = mean(travg_rew(:,twin_prerew),2);
prerewmag_all{ep} = cat(1,prerewmag_all{ep},prerewmag);
negcs_ix = find(corrs_travg<0);
poscs_ix = find(corrs_travg>0);
signeg_rs = reshape(squeeze(mean(cursigs(tr_rew,negcs,twin(twinant)),2))',[],1);
signeg = squeeze(mean(cursigs(tr_rew,negcs,twin(twinant)),2))';
corrs_singtr_cellavg = corr(signeg,twin_s(twinant)');
[~,ixtrsrt] = sort(corrs_singtr_cellavg);
travg_rew_neg_all{ep} = cat(1,travg_rew_neg_all{ep},travg_rew(negcs_ix,:));
travg_rew_pos_all{ep} = cat(1,travg_rew_pos_all{ep},travg_rew(poscs_ix,:));
travg_omit_neg_all{ep} = cat(1,travg_omit_neg_all{ep},travg_omit(negcs_ix,:));
travg_rew_all{ep} = cat(1,travg_rew_all{ep},travg_rew);
lickrate = curd.midAlgn.lick;
lickrate = diff(lickrate,[],2)>0;
lickrate= filter(ones(20,1)/20,1,lickrate')'./curd.dtb;
trlickavg_rew_all{ep} = cat(1,trlickavg_rew_all{ep},lickrate(:,twin_b));
end
end
%% S7F, 7D, S7G, 7E, S7I, S7H
figure;
if 1
for ep = 1:nep
if ep==2
winprerew = winfn(twin_s,[1.9,2]);
else
winprerew = winfn(twin_s,[0.9,1]);
end
[junk,ixsrt] = sort(mean(travg_rew_all{ep}(:,winprerew),2));
posline = find(junk>0,1);
subaxis(1,2,ep);
imagesc(twin_s,[],travg_rew_all{ep}(ixsrt,:),[-0.3 0.6])
hold on;
yl = ylim;
plot([0 0],yl,'k')
plot([1 1],yl,'r--')
plot([2 2],yl,'r--')
plot(xlim,posline*[1 1],'k')
end
colorbar_noresize
figureformat_forsaving
end

ccplt = [0.2 0.2 0.2; 0 0.5 0.5]; %0.5 0.5 0; 
figure; 
for ep = 1:nep
errorbar_shadeSEM(twin_s,travg_rew_neg_all{ep}',ccplt(ep,:));
hold on;
axis tight
end
yl = ylim; 
plot([0 0],yl,'k')
plot([1 1],yl,'k--')
plot([2 2],yl,'k--')
nc_ep = cellfun(@(x)size(x,1),travg_rew_neg_all);
nc_ep_all = cellfun(@(x)size(x,1),travg_rew_all);
title(sprintf("%d %.3g%%, %d %.3g%%",nc_ep(1),nc_ep(1)/nc_ep_all(1)*100,nc_ep(2),nc_ep(2)/nc_ep_all(2)*100))
figureformat_forsaving

histogram_cell(corrs_travg_all)
axis tight;
hold on;
ylim([0 1])
plot([0  0],ylim,'k')
xlabel('Delay correlation with time (r)')

corrs_travg_all_2s_r2 = cellfun(@(x)x.^2,corrs_travg_all_2s,unif=0);
barwitherr_stde(corrs_travg_all_2s_r2,dots=true)
axis tight;
hold on;
xlabel('Delay correlation with time (R2)')

histogram_cell(corrs_travg_all_2s_r2)
axis tight;
hold on;
ylim([0 1])
plot([0  0],ylim,'k')
xlabel('Delay correlation with time (r)')

histogram_cell(prerewmag_all)
signrank_disp(prerewmag_all{1});
signrank_disp(prerewmag_all{2});
axis tight;
hold on;
ylim([0 1]);
plot([0  0],ylim,'k')
xlabel('Pre-reward fluorescence (zsc)')