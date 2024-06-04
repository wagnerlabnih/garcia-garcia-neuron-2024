%% plot example data (Fig 1)
algnuse = "rewAlgn";
muse = "R44";
muse_ix = find(mice.name==muse);
d_use = find(mice.days{muse}.date=="112720");
basedir = "Z:\GrC-CF data\Main";
wd = fullfile(basedir,muse,mice.days{muse}.date(d_use));
curdate = mice.days{muse}.date(d_use);
curd = mice.data{muse}{d_use};
disp(muse),disp(curdate)
tpltCb = winfn([-2 2],curd.tmpxCb);

winuse = winfn([-0.3 0.3],curd.tmpxCb);
alltrs = 1:length(curd.rewarded);
rewtrs = find(curd.rewarded);
[~,plotsort_grc] = max(mean(curd.(algnuse).sigFilt_GrC(rewtrs,:,tpltCb),1),[],3);
[~,plotsortix_grc] = sort(plotsort_grc);
[~,plotsort_cf] = max(mean(curd.(algnuse).sigFilt_CF(rewtrs,:,tpltCb),1),[],3);
[~,plotsortix_cf] = sort(plotsort_cf);

figure(pos=[1335          51         560         420])
    subaxis(1,2,1);
imagesc(tpltCb_s,1:curd.nIC_GrC,squeeze(mean(curd.(algnuse).sigFilt_GrC(rewtrs,plotsortix_grc,tpltCb),1)),[-0.5 0.7])
hold on;plot([0 0],ylim,'k')
subaxis(1,2,2);
imagesc(tpltCb_s,1:curd.nIC_CF,squeeze(mean(curd.(algnuse).sigFilt_CF ...
    (rewtrs,plotsortix_cf,tpltCb),1)),[-0.3 1.5]); 
hold on;plot([0 0],ylim,'k')
title(muse+" "+curdate)
drawnow;

avim_grc = curd.f0_grc;
avim_cf = curd.f0_cf;
clim_grc = [-10 225];
avim_grc8 = uint8im(avim_grc,clim_grc); 
avim_grc8c = repmat(avim_grc8,[1 1 3]); avim_grc8c(:,:,[1 3])=0;
clim_cf = [-10 125];
avim_cf8 = uint8im(avim_cf,clim_cf);
avim_cf8c = repmat(avim_cf8,[1 1 3]); avim_cf8c(:,:,[2])=0;

% R44 Fig 1 GrC subset: 
csubs_mv{1}=[62 16 40 113 81 85 37 50 69 119];
% R44 Fig 1 CF subset: 
csubs_mv{2} = [14 25 35 15 43 44 26 18 19 42];
% pause
ncuse = length(csubs_mv{2});
ccall = ccall(randperm(ncuse*2),:);
ccmult1 = ccall(1:ncuse,:);
ccmult2 = ccall((ncuse+1):end,:);
sig1 = curd.sigFilt_GrC(csubs_mv{1},:);
sig2 = curd.sigFilt_CF(csubs_mv{2},:);
sig4 = round(curd.frameCb(curd.rewtimes)/curd.timesharefac);
sig3 = round(curd.frameCb(curd.truestart)/curd.timesharefac);

stackedTracePlot([sig2; sig1],dtimCb,colors=[ccmult1; ccmult2],events={sig3,sig4},one_t_range=true,scale=true);

ccgrc = repmat([0.9 0.9 .8],curd.nIC_GrC,1);
ccgrc(csubs_mv{1},:) = ccmult1;
coloredMaskPlot(curd.ICmat_GrC,meanim=double(avim_grc8),colors=ccgrc); axis image; truesize; axis off; 
cccf = repmat([.3 .3 .25],curd.nIC_CF,1);
cccf(csubs_mv{2},:) = ccmult2;
tmp = double(curd.ICmat_CF);
for ccur = 1:curd.nIC_CF, 
    tmp2 = tmp(:,:,ccur); 
    tmp2 = (tmp2 - min(tmp2(:)))./sum(tmp2(:));
    tmp2(isinf(tmp2) | isnan(tmp2) | tmp2<0) = 0;
    tmp(:,:,ccur) = tmp2; 
end
tmp = clearMaskBackground(double(tmp),dispresults=false);
coloredMaskPlot(tmp,meanim=double(avim_cf8),colors=cccf); axis image; truesize; axis off;

clear muse d_use plotallcells_trialavg plotmasks_and_traces muse muse_ix wd curdate curd ccgrc cccf sig1 sig2 sig3 sig4 ccall ccmult
clear avim_grc avim_cf clim_grc avim_grc8 avim_grc8c avim_grc8c clim_cf avim_cf8 avim_cf8c avim_cf8c csubs_mv csubs_rew
%% Example body trajectory S3J and mean
for m = 1:nmice
    for d = 1:mice.ndays(m)
        curd = mice.data{m}{d};
        if mice.days{m}.epoch(d)=="expert" && isfield(curd,'goodDLC')
            disp(mice.name(m)+" "+mice.days{m}.date(d)+" "+nnz(isnan(curd.midAlgn.dlcnew(curd.goodDLC,:,:))));
        end
    end
end
%
m = 4;
d = 7;
dlccoorduse = [1:2 5:8];
curd = mice.data{m}{d};
curd.rewarded = logical(curd.rewarded);
tmp = curd.midAlgn.dlcnew;
truse = find(curd.goodmvdir & curd.goodDLC);
truse_rew = truse(find(curd.rewarded(truse)));
truse_omit = truse(find(~curd.rewarded(truse)));
cc = rand(size(tmp,3),3)*0.8;
tmpxx = curd.tmpxDLCnew;
ixuse = winfn([-3 3],tmpxx);
figure(pos=[1063         150         425         710]);
for k = 1:length(dlccoorduse)
    kk = dlccoorduse(k);
    tmp3 = zscore(squeeze(tmp(truse,:,kk)),0,'all');
    errorbar_shadeSEM(tmpxx(ixuse),(squeeze(tmp3(curd.rewarded(truse),ixuse-2)))'+(k-1),cc(k,:));
    hold on;
end
yticks(0:(length(dlccoorduse)-1))
yticklabels(dlclabls(dlccoorduse))
ylim tight;
hold on; plot([0 0],ylim,'k');
hold on; plot([1 1]*meanmvsttime_all,ylim,'k--');
title(sprintf("%s %s m%d d%d, %d trials",mice.name(m),mice.days{m}.date{d},m,d,size(truse_rew,1)))
%% Cell stats Fig S1
nGrC_all = [];
nCF_all = [];
nGrC_ep = cell(1,4);
nCF_ep = cell(1,4);
nGrC_ep_chron = cell(1,4);
nCF_ep_chron = cell(1,4);
nGrC_chron = {[],[]};
nCF_chron = {[],[]};
eplist= ["day1","novice","mid","expert"];
for m = 1:nmice
    for d = 1:mice.ndays(m)
        nGrC_all(end+1) = mice.data{m}{d}.nIC_GrC;
        nCF_all(end+1) = mice.data{m}{d}.nIC_CF;
        curep = find(eplist==mice.days{m}.epoch(d));
        if mice.days{m}.chronDay(d)==min(mice.days{m}.chronDay), curep = 1; end
        nGrC_ep{curep}(end+1,1) = mice.data{m}{d}.nIC_GrC;
        nCF_ep{curep}(end+1,1) = mice.data{m}{d}.nIC_CF;
        if ~isnan(mice.days{m}.chronDay(d))
            nGrC_ep_chron{curep}(end+1,1) = mice.data{m}{d}.nIC_GrC;
            nCF_ep_chron{curep}(end+1,1) = mice.data{m}{d}.nIC_CF;
        end
        if ~isnan(mice.days{m}.chronDay(d))
            nGrC_chron{1}(end+1,:) = mice.data{m}{d}.nIC_GrC;
            nCF_chron{1}(end+1,:) = mice.data{m}{d}.nIC_CF;
        else
            nGrC_chron{2}(end+1,:) = mice.data{m}{d}.nIC_GrC;
            nCF_chron{2}(end+1,:) = mice.data{m}{d}.nIC_CF;
        end
    end
end
fprintf("all GrC/CF cell#s: ");ranksum_disp(nGrC_all,nCF_all);
fprintf("expt GrC cell#s: ");mean_stde_disp(nGrC_ep{4});
fprintf("expt CF cell#s: ");mean_stde_disp(nCF_ep{4});
ranksum_disp(nGrC_chron);
ranksum_disp(nCF_chron);
barwitherr_stde(nGrC_ep_chron,dots=true,overlap=1.3);
eplbls = ["Day 1","2-3","4-6","7-11"];
box off; xticklabels(eplbls); 
ylabel("nGrCs")
xlim([0.5 4.5])
barwitherr_stde(nCF_ep_chron,dots=true,overlap=1.3);
eplbls = ["Day 1","2-3","4-6","7-11"];
box off; xticklabels(eplbls); 
ylabel("nCFs")
xlim([0.5 4.5])
fprintf("GrC chron cell#s anova: ");anova_disp(nGrC_ep_chron);
fprintf("CF chron cell#s anova: ");anova_disp(nCF_ep_chron);

barwitherr_stde(nGrC_chron,dots=true,overlap=1.8);
eplbls = ["Day 1","2-3","4-6","7-11"];
box off; xticklabels(["chronic" "single"]); 
ylabel("nGrCs")
xlim([0.5 2.5])
ranksum_disp(nGrC_chron);

barwitherr_stde(nCF_chron,dots=true,overlap=1.5);
eplbls = ["Day 1","2-3","4-6","7-11"];
box off; xticklabels(["chronic" "single"]); 
ylabel("nCFs")
xlim([0.5 2.5])
fprintf("CF chron vs single cell#s: ");ranksum_disp(nCF_chron);
%% Fig S2O, Show single trials for some example GrCs
exGrCs = [5 7 52; 19 9 29; 17 7 115; 2 8 41];
figure(pos=[ 1771         473         275         810]);
for c = 1:4
    m = exGrCs(c,1); d=exGrCs(c,2); cc=exGrCs(c,3);
subaxis(4,1,c,'MT',0.05,'MB',0.05,'SV',0.01);
curd = mice.data{m}{d};
[tmpxUseCur,tmpxUseCur_s] = winfn(curd.tmpxCb,[-2 2]);
tr_rew = find(curd.mvlen>7 & curd.goodmvdir & curd.rewarded);
tr_o = find(curd.mvlen>7 & curd.goodmvdir & ~curd.rewarded);
nplto = min([length(tr_o),length(tr_rew),40]);
npltr = max([nplto,20]);
trrewuse = shufflev(tr_rew); trrewuse = trrewuse(1:npltr);
tr_ouse = shufflev(tr_o); tr_ouse = tr_ouse(1:nplto);
tmpsig = squeeze(curd.rewAlgn.sigFilt_GrC([trrewuse;tr_ouse],cc,tmpxUseCur));
tmpsigf = butter_filtfilt(2,0.3,tmpsig')';
ccplt = [repmat([0 0 1],npltr,1);repmat([1 0 0],nplto,1)];
rp = randperm(nplto+npltr);
for k = 1:(nplto+npltr)
plot(tmpxUseCur_s,tmpsigf(rp(k),:),col=ccplt(rp(k),:));
hold on;
end
hold on; plot([0 0],ylim,'k-');
hold on; plot(-[meanmvtime_all meanmvtime_all],ylim,'k-');
ylabel(strvec(exGrCs(c,:))+", nr="+npltr+", no="+nplto);
figureformat_forsaving
xlim tight;
end
%% GrC peak finding
for m = 1:nmice
    for d = 1:mice.ndays(m)
        curd = mice.data{m}{d};
        nc = curd.nIC_GrC;
        mice.data{m}{d}.grcpks = false(nc,curd.ntimCb);
        evtThresh = 2.96;
        for c = 1:nc
            [pkmags,pklocs,pkwidths,pkproms] = findpeaks(curd.sigFilt_GrC(c,:),1/curd.dtimCb,MinPeakHeight=evtThresh,MinPeakDistance=0.2,MinPeakProminence=0.5);
            pklocs = round(pklocs/curd.dtimCb);
            mice.data{m}{d}.grcpks(c,pklocs)=true;
        end
        fprintf("mouse %d day %d\n",m,d)
    end
end
%% S1F,G event rates and CV
frcv = cell(1,4);
frcvcat = cell(1,2); fc=[.8 .8 1; 1 .8 .8];
for m = 1:nmice
    for d = 1:mice.ndays(m)
        curd = mice.data{m}{d};
        curd.grcpks = curd.grcpks';
        pknms = ["grcpks","spMat_CF"];
        curmn = {}; curcv = {}; a=[];
        if m==3 && d==5
            figure; 
        end
        for c = 1:2
            curpks = curd.(pknms(c));
            curpks = butter_filtfilt(2,0.2,curpks);
            curmn{c} = (mean(curpks,1)./curd.dtimCb)';
            curcv{c} = curmn{c} ./ (std(curpks,[],1)./curd.dtimCb)';
            if m==3 && d==5
            scatter(curmn{c},curcv{c},30,markeredgecol=[0 0 0],markerfacecol=fc(c,:)); hold on;
            title(strvec([m d curd.nIC_GrC curd.nIC_CF]));
            end
            frcvcat{c} = cat(1,frcvcat{c},[curmn{c} curcv{c}]);
            frcv{(c-1)*2+1} = cat(1,frcv{(c-1)*2+1},mean(curmn{c}));
            frcv{(c-1)*2+2}=cat(1,frcv{(c-1)*2+2},mean(curcv{c}));
        end
    end
end
%
barwitherr_stde(frcv,dots=true,overlap=2)
%% repeat imaging stats
repeatmice = find(arrayfun(@(x)mice.days{x}.epoch(1)=="novice",1:nmice));
ndays_repeat = arrayfun(@(x)nnz(~isnan(mice.days{x}.trueDay)),repeatmice);
mean_stde_disp(ndays_repeat)
ntruedays_repeat = arrayfun(@(x)(max(mice.days{x}.trueDay)-min(mice.days{x}.trueDay)+1),repeatmice);
mean_stde_disp(ntruedays_repeat)
expertdays = arrayfun(@(x)mice.days{x}.chronDay(mice.days{x}.epoch=="expert"),1:nmice,'UniformOutput',false);
expertchronday = arrayfun(@(x)find(mice.days{x}.epoch=="expert"&mice.days{x}.chronDay==max(mice.days{x}.chronDay)),1:nmice,'UniformOutput',false);

numexpertsessions = arrayfun(@(x)nnz(mice.days{x}.epoch=="expert"&(mice.days{x}.chronDay==max(mice.days{x}.chronDay)|isnan(mice.days{x}.chronDay))),1:nmice);
mean_stde_disp(numexpertsessions);

total_days_repeatmice = arrayfun(@(x)mice.ndays(x),repeatmice);
mean_stde_disp(total_days_repeatmice)
%% Cross-correlations S1H,I
allxc = cell(1,4);
scnt=zeros(1,4);
for m = 1:nmice
    for d = 1:mice.ndays(m)
        curep = find(mice.days{m}.epoch(d)==eplist);
        if mice.days{m}.trueDay(d)==min(mice.days{m}.trueDay), curep=1; end
        if curep==4
            if isnan(mice.days{m}.chronDay(d)) || mice.days{m}.chronDay(d)==max(mice.days{m}.chronDay), skip = false;
            else, skip = true;
            end
        else, skip=false;
        end

        if ~skip
        curd = mice.data{m}{d};
        maxlag = round(1/curd.dtimCb);
        nc = curd.nIC_GrC;
        for c = 1:nc
            % [curxc,lags]=xcorr(mean(curd.spMat_CF,2)',curd.sigFilt_GrC(c,:),maxlag,'coeff');
            sig1 = mean(curd.sigFilt_CF,1)';
            sig1 = sig1-mean(sig1);
            sig2 = curd.sigFilt_GrC(c,:);
            sig2=sig2-mean(sig2);
            [curxc,lags]=xcorr(sig1,sig2,maxlag,'coeff');
            curlags = lags*curd.dtimCb;
            if m==1 & d==1, lagtmplt = curlags; end
            curxci = interp1(curlags,curxc,lagtmplt,"linear","extrap");
            allxc{curep} = cat(1,allxc{curep},curxci);
        end
        disp([m d]);
        scnt(curep)=scnt(curep)+1;
        end
    end
end
disp("scnt: "+scnt)
%
figure; errorbar_shadeSEM(lagtmplt,allxc{4}',[0 0 0])
hold on;plot([0 0],ylim,'k'); axis tight;
xlabel("Lag, GrC shift relative to CF (s)")
ylabel("Cross-correlation (r)")

maxes = cellfun(@(x)max(x,[],2),allxc([1 4]),unif=0);
histogram_cell(maxes(2),cdf=true);
xlabel("Peak cross-correlation (r)");
ylabel("Cumulative frac of GrCs")