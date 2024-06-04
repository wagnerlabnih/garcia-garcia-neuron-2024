%%
clear; clc; close all;
load("C:\Users\wagnermj\Desktop\data to share\learning_1s_to_2s_GrC_CF.mat");
%% group/compute quantities
rng(0);
ndirs = cellfun(@(x)length(x),groups);
ngrc =  cellfun(@(g)(cellfun(@(x)x.nIC_GrC,g)),groups,unif=false);
ng = length(groups);
titstrs = ["1-s expert","2-s novice / 1-s expert","2-s expert"];
RESTRICT_TO_INITIAL_NOVICE_TRIALS = true;
tmpxCb_tmplt = groups{3}{1}.tmpxCb;
dtimCb_tmplt = groups{3}{1}.dtimCb;
algnuse = "midAlgn";
xpltwin = [-1 4];
[xplt,xplt_s] = winfn(tmpxCb_tmplt,xpltwin);
nplt = length(xplt);
[rewgrcs, omitlick, omitgrcs, rewCFs, omitCFs, rewdels, s_ix] = deal(cell(1,ng));
for g= 1:ng
    [rewgrcs{g}, omitgrcs{g}, rewCFs{g}, omitCFs{g}] = deal(zeros(0,nplt));
for i = 1:ndirs(g)
    curd = groups{g}{i};
    curd.rewdel = (curd.rewtimes-curd.midpt)*curd.dtb;
    if g>1, rewdeluse = 1.75; else, rewdeluse = 0.75; end
    rewuse = find(curd.rewarded & curd.goodmvdir & curd.rewdel>rewdeluse);
    omituse = find(~curd.rewarded & curd.goodmvdir & curd.mvlen>6 & curd.rewdel>rewdeluse);

    test = filter(ones(300,1)/300,1,double(curd.rewAlgn.lick)')';
    testwin = winfn(curd.tmpxx,[-1.75 1.75]);
    valid_lick = find(max(test(:,testwin),[],2)<0.99);
    validtrials = valid_lick;
    omituse_lick = intersect(omituse,valid_lick);

    fprintf("g=%d i=%d %d rew %d omit\n",g,i,numel(rewuse),numel(omituse));
    if RESTRICT_TO_INITIAL_NOVICE_TRIALS
    if g==2,
        rewuse = rewuse(1:min(length(rewuse),50));
        omit = omituse(1:min(length(omituse),50));
    else
        rewuse = randsample(rewuse,min(length(rewuse),50));
        omit = randsample(omituse,min(length(omituse),15));
    end
    end

    [xplt_cur,xplt_cur_s] = winfn(curd.tmpxCb,xpltwin);
    [xplt_cur_b,xplt_cur_b_s] = winfn(curd.tmpxx,xpltwin);
    cur_grcrew = squeeze(mean(curd.(algnuse).sigFilt_GrC(rewuse,:,xplt_cur),1));
    cur_grcomit = squeeze(mean(curd.(algnuse).sigFilt_GrC(omituse,:,xplt_cur),1));
    tmpspvec = double(reshape(permute(curd.(algnuse).sp_CF,[3 1 2]),[],size(curd.(algnuse).sp_CF,2)));
    spkern = round(0.1 / curd.dtimCb);
    spkern_s = spkern*curd.dtimCb;
    tmpspvecf = filtfilt(ones(spkern,1)/spkern/spkern_s,1,tmpspvec);
    tmpspvecf = zscore(tmpspvecf,[],1);
    tmpspvecf = reshape(tmpspvecf,size(curd.(algnuse).sp_CF,3),size(curd.(algnuse).sp_CF,1),size(curd.(algnuse).sp_CF,2));
    curcfsig = permute(tmpspvecf,[2 3 1]);
    tmp3 = squeeze(mean(curcfsig(rewuse,:,xplt_cur),1));
    tmp4 = squeeze(mean(curcfsig(omituse,:,xplt_cur),1));
    if numel(xplt)~=numel(xplt_cur)
    cur_grcrew = interp1(xplt_cur_s,cur_grcrew',xplt_s,'linear','extrap')';
    cur_grcomit = interp1(xplt_cur_s,cur_grcomit',xplt_s,'linear','extrap')';
    tmp3 = interp1(xplt_cur_s,tmp3',xplt_s,'linear','extrap')';
    tmp4 = interp1(xplt_cur_s,tmp4',xplt_s,'linear','extrap')';
    end
    rewgrcs{g} = cat(1,rewgrcs{g},cur_grcrew);
    s_ix{g} = cat(1,s_ix{g},[(1:ngrc{g}(i))' repmat(i,ngrc{g}(i),1)]);
    filtwin = 60;
    lickraw = double(curd.midAlgn.lick);
    lickrate = (diff(lickraw,[],2)>0)./0.005;
    curlick = filter(ones(filtwin,1)/filtwin,1,lickrate')';
    omitlick{g} = cat(1,omitlick{g},curlick(omituse_lick,xplt_cur_b));
    rewdels{g} = cat(1,rewdels{g},curd.rewdel([rewuse; omituse]));
    omitgrcs{g} = cat(1,omitgrcs{g},cur_grcomit);
    rewCFs{g} = cat(1,rewCFs{g},tmp3);
end
end
%% Fig 4I-P
f1 = figure(pos=[172    68   400*ng   420]);
f2 = figure(pos=[172    68   400   420]);
nf2 = [];
f2cf = figure(pos=[172    68   400   420]);
ccplt = [0.2 0.2 0.2; 0.5 0.5 0; 0 0.5 0.5];
grc_scale = [-0.4 0.6];
cf_scale = [-0.4 0.4];
[q,q2,qt,long_active_sfrac] = deal(cell(ng,1));
hlines = []; hlines_o = [];
rewtimes = [1 2 2];
dtb = 0.005;
long_active = 1.6;
for g= 1:ng
    ix_2sdel = winfn(xplt_s,[0 2]);
    ix_2sdel_b = winfn(curd.tmpxx(xplt_cur_b),[0 2]);
    premv_s = [-1 0];
    ix_premv = winfn(xplt_s,premv_s);
    ix_premv_b = winfn(curd.tmpxx(xplt_cur_b),premv_s);
    ix_prerew = winfn(xplt_s,rewtimes(g)-[0.3 0.03]);
    ix_earlydel = winfn(xplt_s,[0.1 0.4]);
    delay_activ_duration = sum(rewgrcs{g}(:,ix_2sdel)-mean(rewgrcs{g}(:,ix_premv),2)>0,2)*dtimCb_tmplt;
    earlydelRise = mean(rewgrcs{g}(:,ix_prerew),2)-mean(rewgrcs{g}(:,ix_earlydel),2);
    delay_activ_cent = zeros(size(rewgrcs{g},1),1);
    for c= 1:size(rewgrcs{g},1)
        delay_activ_indeces = find(rewgrcs{g}(c,ix_2sdel)-mean(rewgrcs{g}(c,ix_premv),2)>0);
        tmp4 = rewgrcs{g}(c,ix_2sdel(delay_activ_indeces));
        delay_activ_cent(c) = sum(tmp4./sum(tmp4).*xplt_s(ix_2sdel(delay_activ_indeces)));
    end
    delay_activ_cent(delay_activ_cent==0)=-inf;
    [qt{g},idxsrt] = sort(delay_activ_cent);
    figure(f1); 
    subaxis(1,ng,g);
    imagesc(xplt_s,[],rewgrcs{g}(idxsrt,:),grc_scale);
    hold on; plot([0 0],ylim,'k',linew=1.25);
    plot(rewtimes(g)*[1 1],ylim,'r',linew=1.25); 
    if g==2
        plot([1 1],ylim,'r--',linew=1.5,col=[1 .5 .5]);
        xlabel("Time rel to movement (s)")
    elseif g==1
        plot([2 2],ylim,'r--',linew=1.5,col=[1 .5 .5]);
    end
    xlim([-1 4])
    title(titstrs(g)+sprintf(" %d GrCs %d sessions",length(rewgrcs{g}),length(groups{g})))

    ix_cf_postrew = winfn(xplt_s,[0 0.2]+rewtimes(g));
    ix_cf_prerew = winfn(xplt_s,[-0.3 -0.03]+rewtimes(g));
    ncf_postrewrise = mean(rewCFs{g}(:,ix_cf_postrew),2) - mean(rewCFs{g}(:,ix_cf_prerew),2);
    cfpostrew_activ = mean(rewCFs{g}(:,ix_cf_postrew),2);
    cfs_elevatedpostrew = find(ncf_postrewrise>0.1);
    [~,idxsrt_cf] = sort(ncf_postrewrise);

    q{g} = delay_activ_duration;
    long_active_bool = delay_activ_duration>long_active;
    dir_ix = arrayfun(@(x)find(s_ix{g}(:,2)==x),1:ndirs(g),UniformOutput=false);
    long_active_sfrac{g} = arrayfun(@(x)sum(long_active_bool(dir_ix{x}))/ngrc{g}(x),1:ndirs(g));
    q2{g} = long_active_bool;
    ixmean = find(earlydelRise>0);
    nf2(:,g) = [numel(ixmean)./numel(earlydelRise)*100; numel(ixmean)];
    nf2cf(:,g) = [numel(cfs_elevatedpostrew)./numel(ncf_postrewrise)*100; numel(cfs_elevatedpostrew)];
    figure(f2);
    hlines(g) = errorbar_shadeSEM(xplt_s,rewgrcs{g}(ixmean,:)',ccplt(g,:)); hold on;
    figure(f2cf);
    hlines_cf(g) = errorbar_shadeSEM(xplt_s,rewCFs{g}(cfs_elevatedpostrew,:)',ccplt(g,:)); hold on;
    ylabel("CF spike rate (zsc)")
end
figure(f1); colorbar_noresize;
fs = [f2 f2cf];
hl = {hlines,hlines_cf};
nfs = {nf2,nf2cf};
for fcur = 1:length(fs)
figure(fs(fcur));
axis tight;
hold on; plot([0 0],ylim,'k'); plot([2 2],ylim,'r--');
plot([1 1],ylim,'r--');
xlim([-1 4])
xlabel("Time rel to movement (s)")
uistack(hl{fcur},'top');
if fcur<2
legend(hl{fcur},cat(1,arrayfun(@(x)sprintf("%s: %1.2g%%, %d/%d GrCs",titstrs(x),nfs{fcur}(1,x),nfs{fcur}(2,x),size(rewgrcs{x},1)),1:ng)));
else
legend(hl{fcur},cat(1,arrayfun(@(x)sprintf("%s: %1.2g%%, %d/%d CFs",titstrs(x),nfs{fcur}(1,x),nfs{fcur}(2,x),size(rewCFs{x},1)),1:ng)));
end
end

barwitherr_stde(q(1:3),ranksumgrps={[1 3],[2 3]});
ylim tight;
% figure; boxplotcell(q);
ylabel("Duration of activity during [0,2] (s)")
xticklabels(titstrs)
guse = [1 2 3];
histogram_cell(q,[],cdf=true,col=ccplt);
xticks(0:.5:2); 
hold on; plot([1 1],ylim,'k--'); plot([2 2],ylim,'k')
ylabel("Fraction of GrCs")
xlabel("Duration of activity during [0,2] (s)"); box off;
legend(titstrs(guse));

% barwitherr_stde(q2(1:3));
barwitherr_stde(long_active_sfrac,dots=true,ranksumgrps={[1 3]});
ylim tight;
xticklabels(titstrs)
ylabel("Fraction of GrCs active for >1.75 s during [0,2]s")


%% decoding time
SCALE_TO_FIRST_SECOND = true;
twins = [0 2];
algnuse = "midAlgn";
tmpxCb_tmplt = groups{3}{1}.tmpxCb;
timedecode_coefs = cell(ng,1);
timedecode_stats_trs = cell(ng,1);
timedecode_stats_trs_m = cell(ng,1);
timedecode_stats_trs_m_01s= cell(ng,1);
timedecode_stats_trs_01s= cell(ng,1);
timedecode_stats_trs_m_12s= cell(ng,1);
timedecode_stats_trs_12s= cell(ng,1);
timedecode_kfl = cell(ng,1);
timedecode_kfl_m = cell(ng,1);
timedecode_output = cell(ng,1);
timedecode_output_trs = cell(ng,1);
twin_tmplt = tmpxCb_tmplt(winfn(twins,tmpxCb_tmplt));
twin_tmplt_ix = winfn(twins,tmpxCb_tmplt);
timedecode_timeax = twin_tmplt;
for g =1:ng
    nd = ndirs(g);
    timedecode_stats_trs_m{g} = cell(nd,1);
    timedecode_stats_trs_m_01s{g} = cell(nd,1);
    timedecode_stats_trs_m_12s{g} = cell(nd,1);
    timedecode_output{g} = cell(nd,1);
    timedecode_coefs{g} = cell(nd,1);
    timedecode_kfl_m{g} =  cell(nd,1);
    if g == 1, deluse = 0.75;
    else, deluse = 1.75; end
    for i = 1:nd
    cur_dtimCb = groups{g}{i}.dtimCb;
    curd = groups{g}{i};
    curd.rewdel = (curd.rewtimes-curd.trueend)*curd.dtb;
    cur_tmpxCb = curd.tmpxCb;
    rewuse = find(curd.rewarded & curd.goodmvdir & curd.rewdel>deluse);
    omituse = find(~curd.rewarded & curd.goodmvdir & curd.mvlen>6 & curd.rewdel>deluse);
    truse = [rewuse; omituse];
%     truse = rewuse;
    twin = winfn(twins,cur_tmpxCb);
    allcursigs = permute(curd.(algnuse).sigFilt_GrC(truse,:,twin),[3 1 2]);
    allcursigs = reshape(allcursigs,length(twin)*length(truse),[]);
    allcurtax = repmat(curd.tmpxCb(twin),[1 length(truse)])';
    [delay_activ_cent,mae,r2] = kfoldpred_fitlm(allcursigs,allcurtax);
    delay_activ_duration = reshape(delay_activ_cent,[length(twin),length(truse)])';

    if SCALE_TO_FIRST_SECOND
        firstsec = winfn(cur_tmpxCb(twin),[0 1]);
        twin_scale = firstsec;
    else
        twin_scale = 1:length(twinant);
    end
    tmpmn = mean(delay_activ_duration(:,twin_scale),1);
    tmpmn = min(tmpmn);
    tmpdec_scaled = delay_activ_duration-tmpmn; 
    bla = mean(tmpdec_scaled(:,twin_scale),1);
    bla = range(bla);
    if bla~=0,
    endpts = cur_tmpxCb(twin(twin_scale([1 end])));
    tmpdec_scaled = tmpdec_scaled./bla*range(endpts);
    delay_activ_duration = tmpdec_scaled;
    else
        pause;
    end

    curdecode_output_interp = interp1(cur_tmpxCb(twin),delay_activ_duration',twin_tmplt,"linear","extrap");
    timedecode_stats_trs_m{g}{i} = (corr(curd.tmpxCb(twin)',delay_activ_duration').^2)';
    timedecode_stats_trs{g} = cat(1,timedecode_stats_trs{g},(corr(curd.tmpxCb(twin)',delay_activ_duration').^2)');    
    if twins(2)==2
    twin_01s = winfn([0 1],cur_tmpxCb);
    twin_01s_s = cur_tmpxCb(twin_01s);
    twin_12s = winfn([1 2],cur_tmpxCb);
    twin_12s_s = cur_tmpxCb(twin_12s);
    delay_activ_cent_01s = delay_activ_cent(ismember(allcurtax,twin_01s_s));
    delay_activ_cent_12s = delay_activ_cent(ismember(allcurtax,twin_12s_s));
    delay_activ_duration_01s = reshape(delay_activ_cent_01s,[length(twin_01s),length(truse)])';
    delay_activ_duration_12s = reshape(delay_activ_cent_12s,[length(twin_12s),length(truse)])';    
    timedecode_stats_trs_m_01s{g}{i} = (corr(curd.tmpxCb(twin_01s)',delay_activ_duration_01s').^2)';
    timedecode_stats_trs_01s{g} = cat(1,timedecode_stats_trs_01s{g},(corr(curd.tmpxCb(twin_01s)',delay_activ_duration_01s').^2)');    
    timedecode_stats_trs_m_12s{g}{i} = (corr(curd.tmpxCb(twin_12s)',delay_activ_duration_12s').^2)';
    timedecode_stats_trs_12s{g} = cat(1,timedecode_stats_trs_12s{g},(corr(curd.tmpxCb(twin_12s)',delay_activ_duration_12s').^2)');    
    end
    timedecode_output{g}{i} = curdecode_output_interp;
    timedecode_output_trs{g} = cat(1,timedecode_output_trs{g},curdecode_output_interp');
    timedecode_kfl{g} = cat(1,abs(allcurtax-delay_activ_cent));
    timedecode_kfl_m{g}{i} = abs(allcurtax-delay_activ_cent);
    mdl2 = fitlm(allcursigs,allcurtax);
    timedecode_coefs{g}{i} = mdl2.Coefficients.Estimate(2:end);
    fprintf("%d %d\n",g,i);
    end
end
%% Fig. S4B-D
guse = 1:3;
ccplt = [0.2 0.2 0.2; 0.5 0.5 0; 0 0.5 0.5];
nguse = length(guse);
histogram_cell(timedecode_stats_trs_01s(guse),cdf=true,col=ccplt,pos=[867   447   377   341])
xlabel("[0,1] s R^2")
ylabel("Cumulative fraction of trials")
histogram_cell(timedecode_stats_trs_12s(guse),cdf=true,col=ccplt,pos=[867   447   377   341])
xlabel("[1,2]s R^2")
ylabel("Cumulative fraction of trials")
histogram_cell(timedecode_stats_trs(guse),cdf=true,col=ccplt,pos=[867   447   377   341])
xlabel("[0,2]s R^2")
ylabel("Cumulative fraction of trials")
% yticks(0:0.05:0.2);
figure(pos=[ 1164         471         203*nguse         193]);
mg = [2 5; 3 5];
ax = [];
nguse = 2; 
for g = 1:2
ax(g) = subaxis(1,nguse,g);
toplt = timedecode_output{mg(g,1)}{mg(g,2)};
rs = randsample(size(toplt,2),40);
size(toplt,2)
plot(twin_tmplt,toplt(:,rs),col=[0.7 0.7 0.7])
axis tight; box off;
hold on; plot(twins,twins,'k')
if g==1, plot([1 1],ylim,'r')
else
    plot([2 2],ylim,'r')
end
plot([0 0],ylim,'k')
if g==2, plot([1 1],ylim,'r--'), end
title(sprintf("mouse %d r2=%.3g",mg(g,1),mean(timedecode_stats_trs_m{mg(g,1)}{mg(g,2)})));
end
linkaxes(ax)

if twins(2)==2
figure;
ccplt = [0.2 0.2 0.2; 0.5 0.5 0; 0 0.5 0.5];
for i = 1:3
errorbar_shadeSEM(twin_tmplt,timedecode_output_trs{i}',ccplt(i,:));
hold on
end
xticks(0:.5:2); yticks([0:.5:2]); axis equal square; xlim([0 2]); ylim([0 2]); 
end
hold on; plot(xlim,ylim,'k');
plot([1 1],ylim,'r--')
plot([2 2],ylim,'r--')
title(sprintf("%d %d %d",arrayfun(@(x)size(timedecode_output_trs{x},1),1:3)))
%% LTD computation
reslts = cell(1,2);
twinant_tmplt = cell(1,2);
twin_sort_tmplt = cell(1,2);
for krun = 1:2
RESTRICT_TO_INITIAL_NOVICE_TRIALS = true;
SHOW_GRC_SORTING = false;
SHOW_WEIGHT_VECS = false;
SHOW_PKJ_OUTPUT = false;
ELEGWIN_S = [-0.15 -0.025];
grcsiguse_ltd = "sigFilt_GrC";
grcsiguse_decode = "sigFilt_GrC";
reslts{krun} = struct([]);
cnt = 1;
if krun==1
    algnuse = "midAlgn";
    SCALE_TO_FIRST_SECOND = true;
    twinant_s = [0 2.1];
    REW_CF_WIN = {[1 1.25],[2 2.25],[2 2.25]};
    twin_sort_s = [-1 3];
elseif krun==2
    algnuse = "rewAlgn";
    SCALE_TO_FIRST_SECOND = false;
    twinant_s = [-2 0];
    REW_CF_WIN = {[0 0.25],[0 0.25],[0 0.25]};
    twin_sort_s = [-2 2];
end
brkptspct = [0 0.05 0.1 0.2 0.4 0.6 0.8 1];
rng(0);
for m = 1:ng
    for d = 1:length(groups{m})
        curd = groups{m}{d};
        curd.rewdel = (curd.rewtimes-curd.midpt)*curd.dtb;
        tmpxCb = curd.tmpxCb;
        dtimCb = curd.dtimCb;
        twin = winfn(twin_sort_s,tmpxCb);
        twin_s = tmpxCb(twin);
        if m==1, deluse=0.75; else deluse=1.75; end
        tr_rew = find(curd.rewarded & curd.goodmvdir & curd.rewdel>deluse);
        tr_omit = find(~curd.rewarded & curd.goodmvdir & curd.mvlen>6 & curd.rewdel>deluse);
        tr_use_decode = tr_rew;
        if RESTRICT_TO_INITIAL_NOVICE_TRIALS
        if m==2
            tr_use_decode = tr_use_decode(1:min(length(tr_use_decode),50));
        else
            tr_use_decode = randsample(tr_use_decode,min(length(tr_use_decode),50));
        end
        end
        tr_use_ltd = tr_rew;
        curspcf = curd.(algnuse).sp_CF(tr_use_ltd,:,twin);
        curspcf_travg = squeeze(mean(curspcf,1));
        curspcf_r = curd.(algnuse).sp_CF(tr_rew,:,twin);
        curspcf_travg_r = squeeze(mean(curspcf_r,1));

        tmpwinpkj1 = winfn(REW_CF_WIN{m}-0.3,twin_s);
        tmpwinpkj2 = winfn(REW_CF_WIN{m},twin_s);
        tmprewwin = winfn(REW_CF_WIN{m},twin_s);
        pkjuse = find((mean(curspcf_travg_r(:,tmpwinpkj2),2)-mean(curspcf_travg(:,tmpwinpkj1),2))>0.0);
        
        if ~isempty(pkjuse)
        curspcf = curspcf(:,pkjuse,:);
        curspcf_travg = curspcf_travg(pkjuse,:);

        reslts{krun}(cnt).pkjuse = [length(pkjuse),curd.nIC_CF];
        reslts{krun}(cnt).day = d;
        reslts{krun}(cnt).travg_spcf = curspcf_travg;
        eleg_win = round(ELEGWIN_S(1)/dtimCb):round(ELEGWIN_S(2)/dtimCb);
        test2 = false(size(curspcf));
        test2_signed = zeros(size(curspcf));
        meanrewresp = prctile(max(mean(curspcf_r(:,pkjuse,tmpwinpkj2),2),[],3),75,1);
        for c = 1:size(curspcf,2)
            for tr = 1:size(curspcf,1)
                cursp = tmprewwin(find(squeeze(curspcf(tr,c,tmprewwin))))';
                test3 = false(length(curspcf(tr,c,:)),1);
                if ~isempty(cursp)
                tmpix = reshape(cursp+eleg_win,1,[]);
                tmpix(tmpix<1 | tmpix>length(test3)) = [];
                test3(tmpix) = true;
                end
                test2(tr,c,:) = test3;
            end
        end
        test4 = permute(test2,[1 4 3 2]);
        curgrc = curd.(algnuse).(grcsiguse_ltd)(tr_use_ltd,:,twin);
        curLTD = bsxfun(@times,curgrc,test4);
        logitrate = prctile(reshape(permute(curgrc,[3 1 2]),[],curd.nIC_GrC),95);
        curLTD(curLTD<0) = 0;
        curLTD = 1./(1+exp(-1./logitrate.*curLTD));
        curLTD_travg = squeeze(mean(curLTD,1));
        curLTD_trPkjavg = squeeze(mean(curLTD_travg,3));
        curLTD_trtimeavg = squeeze(mean(curLTD_travg,2));
        curLTD_trtimeavg = (curLTD_trtimeavg-mean(curLTD_trtimeavg,1))./sum(curLTD_trtimeavg,1);
        curLTD_trtimePkjavg = squeeze(mean(curLTD_trtimeavg,2));
        [~,LTDsort] = sort(-curLTD_trtimePkjavg);
        twin_sort = winfn(twin_sort_s,tmpxCb);
        twin_sort_tmplt{krun} = winfn(twin_sort_s,tmpxCb_tmplt);
        curgrc_tmpsort = curd.(algnuse).(grcsiguse_ltd)(tr_use_ltd,:,twin_sort);
        curgrcavg = squeeze(mean(curgrc_tmpsort,1));
        curgrcinterp = interp1(tmpxCb(twin_sort),curgrcavg',tmpxCb_tmplt(twin_sort_tmplt{krun}),"linear","extrap")';
        brkpts = round(curd.nIC_GrC*brkptspct);
        for k = 1:(length(brkpts)-1)
            curix = [(brkpts(k)+1) brkpts(k+1)];
            reslts{krun}(cnt).grcLTDquartile{k} = curgrcinterp(LTDsort(curix(1):curix(2)),:);
        end
        curLTD_trtimePkjavg = -curLTD_trtimePkjavg;
        reslts{krun}(cnt).LTD_trs_time_Pkj_avg = curLTD_trtimePkjavg;
        tmpgrc = permute(curd.(algnuse).(grcsiguse_decode),[3 1 2]);
        curkern = round(0.07/dtimCb);
        tmpgrc = filtfilt(ones(curkern,1)./curkern,1,double(tmpgrc));
        tmpgrc = permute(tmpgrc,[2 3 1]);
        twinant = winfn(twinant_s,tmpxCb);
        twinant_tmplt{krun} = winfn(twinant_s,tmpxCb_tmplt);
        twins = {twinant};
        finantix = {};
        for k = 1
        tmpgrc2 = squeeze(mean(tmpgrc(:,:,twins{k}),1));
        finantix{k} = zeros(size(tmpgrc2,1),1);
        for c = 1:size(tmpgrc2,1)
            hicut = max(mean(tmpgrc(:,c,:),1),[],3)*0.5;
            junk = find(tmpgrc2(c,:)>hicut,1,'last');
            if isempty(junk), junk = 1; end
            finantix{k}(c) = tmpxCb(twins{k}(junk));
        end
        [~,cursrt] = sort(curLTD_trtimePkjavg);
        [~,cursrt2] = sort(finantix{k}); if k == 1
        reslts{krun}(cnt).anticOfftimeLTD_corr = corr(-curLTD_trtimePkjavg,finantix{k});
        fprintf("antic off-time %d %d %.2g\n",m,d,reslts{krun}(cnt).anticOfftimeLTD_corr)
        if isnan(reslts{krun}(cnt).anticOfftimeLTD_corr)
             reslts{krun}(cnt).anticOfftimeLTD_corr = 0;
        end
        end
        tmpgrc = squeeze(mean(curd.(algnuse).(grcsiguse_decode)(tr_use_decode,:,twin),1));
        [~,LTDsort2] = sort(timedecode_coefs{m}{d},'descend');
        [reslts{krun}(cnt).regwveccorr,p] = corr(curLTD_trtimePkjavg,-timedecode_coefs{m}{d});
        reslts{krun}(cnt).ngrc = length(curLTD_trtimePkjavg);
        reslts{krun}(cnt).travg_grc = tmpgrc;
        test = sum(tmpgrc.*curLTD_trtimePkjavg,1);
        reslts{krun}(cnt).simplPkjCell_trAvg = test;

        tmpgrc = curd.(algnuse).(grcsiguse_decode)(tr_use_decode,:,twinant);
        tmpgrcf = filtfilt(curkern,1,reshape(permute(double(tmpgrc),[3 1 2]),[],size(tmpgrc,2)));
        tmpgrcf = permute(reshape(tmpgrcf,size(tmpgrc,3),size(tmpgrc,1),size(tmpgrc,2)),[2 3 1]);
        tmpgrc= tmpgrcf;
        tmpPkj = squeeze(sum(bsxfun(@times,tmpgrc,permute(curLTD_trtimePkjavg,[3 1 2])),2));
        if SCALE_TO_FIRST_SECOND
            firstsec = winfn(tmpxCb(twinant),[0 1]);
            twin_scale = firstsec;
        else
            twin_scale = 1:length(twinant);
        end
        tmpmn = mean(tmpPkj(:,twin_scale),1);
        tmpmn = max(tmpmn);
        tmpPkj_scaled = tmpPkj-tmpmn; 
        bla = mean(tmpPkj_scaled(:,twin_scale),1);
        bla = range(bla);
        endpts = tmpxCb(twinant(twin_scale([1 end])));
        tmpPkj_scaled = tmpPkj_scaled./bla*range(endpts);
        tmpPkj_travg = mean(tmpPkj,1);
        reslts{krun}(cnt).Pkj_singlTrs = tmpPkj_scaled;
        tmpPkj_interp = interp1(tmpxCb(twinant),tmpPkj_scaled',tmpxCb_tmplt(twinant_tmplt{krun}),"linear","extrap")';
        reslts{krun}(cnt).Pkj_interp = tmpPkj_interp;
        reslts{krun}(cnt).Pkj_r_interp = tmpPkj_interp;%(1:length(tr_rew),:);

        tmpTime = reshape(repmat(-(tmpxCb(twinant)),[size(tmpPkj,1) 1])',[],1);
        tmpPkj_rs = reshape(tmpPkj_scaled',[],1);
        reslts{krun}(cnt).PkjTime_corr_stScale = corr(tmpTime,tmpPkj_rs);
        reslts{krun}(cnt).PkjTime_corr_trs = corr(-(tmpxCb(twinant))',tmpPkj_scaled')';
        reslts{krun}(cnt).opt_corr_trs = timedecode_stats_trs_m{m}{d};
        reslts{krun}(cnt).PkjTime_corr_pctOfOptimal = ...
            min(reslts{krun}(cnt).PkjTime_corr_stScale^2/mean(timedecode_stats_trs_m{m}{d}),2);
        reslts{krun}(cnt).lineregreslts{krun} = sqrt(mean(timedecode_stats_trs_m{m}{d}));
        fprintf("%d %d Pkj-time corr st scaling %.2g, %% of optimal %.2g\n",...
            m,d,reslts{krun}(cnt).PkjTime_corr_stScale,reslts{krun}(cnt).PkjTime_corr_pctOfOptimal);
        reslts{krun}(cnt).ep=m;
        cnt = cnt+1;
        end
        end
    end
end
end
%% Fig 6F,G
guse = 1:3;
resltgrps = arrayfun(@(x)find(cat(1,reslts{1}.ep)==x),guse,'UniformOutput',false);
nrg = length(guse);
resgrps = cellfun(@(x)reslts{1}(x),resltgrps,'UniformOutput',false);
grpcorrs = cellfun(@(x)cat(1,x.PkjTime_corr_trs),resgrps,'UniformOutput',false);
grpr2 = cellfun(@(x)x.^2,grpcorrs,'UniformOutput',false);
grppkjoutputs = cellfun(@(x)cat(1,x.Pkj_r_interp),resgrps,'UniformOutput',false);
barwitherr_stde(grpr2,color=ccplt);
histogram_cell(grpr2,0:.001:1,cdf=true,color=ccplt)
ccplt = [0.2 0.2 0.2; 0.5 0.5 0; 0 0.5 0.5];
figure(pos=[1139         411         443         492]); 
for rg = 1:nrg
errorbar_shadeSEM(tmpxCb_tmplt(twinant_tmplt{1}),grppkjoutputs{rg}',ccplt(rg,:));
hold on; 
end
axis tight;
title(mean_stde_disp(grpcorrs));
plot([1 1],ylim,'r--'); 
plot([0 2],[0 -2],'k')
if tmpxCb_tmplt(twinant_tmplt{1}(end))>1.9, plot([2 2],ylim,'r--');  end
figureformat_forsaving
%% Fig 5M, S4N
twinall = [0,3];
twin02 = [0,2.1];
ixall = winfn(twinall,tmpxCb_tmplt(twin_sort_tmplt{1}));
ix02 = winfn(twin02,tmpxCb_tmplt(twin_sort_tmplt{1}));
resltgrps = arrayfun(@(x)find(cat(1,reslts{1}.ep)==x),guse,'UniformOutput',false);
fdat = figure(pos=[968   198   915   259]);
ffit = figure(pos=[968   198   915   259]);
cctmp = cool(length(brkptspct));
[hlines,hlf,herrs,ixmx,risefall] = deal([]);
titstr = "";
for k = 1:(length(brkptspct)-1)
    allqtile = arrayfun(@(x)cat(1,x.grcLTDquartile{k}),reslts{1}(resltgrps{3}),unif=false);
    allqtile = cat(1,allqtile{:});
    if any(isnan(allqtile(:))), pause; end
    allqtile_a = allqtile(:,ixall);
    tp = tmpxCb_tmplt(twin_sort_tmplt{1}(ix02));
    aqm_o = mean(allqtile(:,ix02),1);
    aqm = (aqm_o-min(aqm_o))/range(aqm_o);

    figure(fdat)
    hlines(end+1)=plot(tp,aqm,col=cctmp(k,:));
    hold on;
    xlim(tmpxCb_tmplt(twin_sort_tmplt{1}([ix02(1) ix02(end)])))
    ylim([-0.05 1.05])

    figure(ffit);
    xt = tmpxCb_tmplt(twin_sort_tmplt{1}(ixall))';
    yt = allqtile_a;
    hlf(end+1)=errorbar_shadeSEM(xt,yt',cctmp(k,:));
    hold on;
    titstr = titstr+sprintf("%1.2g-%1.2g, %d; ",brkptspct(k),brkptspct(k+1),size(allqtile,1));
    % pause;
end
figure(ffit)
xlabel("Time relative to movement (s)");
ylabel("Fluorescence (zsc)")
ylim tight;
plot([2 2],ylim,'k')
uistack(hlf,"top")
title(titstr)
figureformat_forsaving

figure(fdat)
xlabel("Time relative to movement (s)");
ylabel("Group activity (normalized)")
plot([2 2],ylim,'k')
ylim tight;
figureformat_forsaving
%% Novice vs expert 2s licking S2H,I
figure(pos=[1138          22         328         420]);
ccplt = [0.2 0.2 0.2; 0.5 0.5 0; 0 0.5 0.5];

errorbar_shadeSEM(xplt_cur_b_s,omitlick{2}',ccplt(2,:))
hold on;
errorbar_shadeSEM(xplt_cur_b_s,omitlick{3}',ccplt(3,:))
xlabel("Time rel to movement (s)"); ylabel("Lick rate (Hz)"); 
hold on; plot([2 2],ylim,'k')
yl=ylim;yl(1)=0;ylim(yl);
winearly = winfn([0.3,0.9],xplt_cur_b_s);
winlate = winfn([2,2.75],xplt_cur_b_s);
cq = cellfun(@(x)(mean(x(:,winlate),2)-mean(x(:,winearly),2)),omitlick([2 3]),unif=0);
histogram_cell(cq,title="S2G",cdf=true,color=ccplt(2:3,:))
hold on; plot([0 0],ylim,'k')