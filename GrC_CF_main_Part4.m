%% decoding time
twins_all = {[-1.1015 0.0],...
         [0 1],... % positive control
         [1 2]}; % negative control
ntw = length(twins_all);
[timedecode_coefs,timedecode_stats,timedecode_kfl,timedecode_output,...
    twin_tmplt,twin_tmplt_ix,]=deal(cell(ntw,1));
epoch_list= ["novice","mid","expert"];
for tw = 1:ntw
twins = twins_all{tw};
twin_tmplt{tw} = tmpxCb_tmplt(winfn(twins,tmpxCb_tmplt));
twin_tmplt_ix{tw} = winfn(twins,tmpxCb_tmplt);
timedecode_coefs{tw} = cell(nmice,1);
timedecode_stats{tw} = cell(nmice,1);
timedecode_kfl{tw} = cell(nmice,1);
timedecode_output{tw} = cell(nmice,1);
if tw==1, 
    timedecode_stats_trs = cell(nmice,1);
    timedecode_stats_wtf = cell(nmice,1); 
end
for m =1:nmice
    cur_dtimCb = mice.data{m}{1}.dtimCb;
    for d = 1:mice.ndays(m)
        curep = find(mice.days{m}.epoch(d)==epoch_list)+1;
        if mice.days{m}.trueDay(d)==min(mice.days{m}.trueDay), curep = 1; end
        curd = mice.data{m}{d};
        cur_tmpxCb = curd.tmpxCb;
        if all(twins==[0 1]) || all(twins==[1 2])
            truse = find(curd.rewarded);
        else
            truse = find(curd.rewarded | (~curd.rewarded & curd.mvlen>=7));
        end
        truse = intersect(truse,find(curd.goodmvdir));
        allcursigs = [];
        allcurtax = [];
        for t = 1:length(truse)
            twin = winfn(twins,cur_tmpxCb);
            curt = truse(t);
            cursigs = squeeze(curd.rewAlgn.sigFilt_GrC(curt,:,:));
            allcursigs = cat(1,allcursigs,cursigs(:,twin)');
            allcurtax = cat(1,allcurtax,curd.tmpxCb(twin)');
        end
        [tmp2,mae,r2] = kfoldpred_fitlm(allcursigs,allcurtax);
        tmp = reshape(tmp2,[length(twin),length(truse)])';
        curdecode_output_interp = interp1(cur_tmpxCb(twin),tmp',twin_tmplt{tw},"linear","extrap");
        timedecode_output{tw}{m}{d} = curdecode_output_interp;
        timedecode_stats{tw}{m}(d) = corr(allcurtax,tmp2)^2;
        if tw==1,
            junk = [corr(curd.tmpxCb(twin)',tmp')' truse];
            junk(isnan(junk)) = 0; timedecode_stats_trs{m}{d}=junk;
            assert(~isnan(mean(timedecode_stats_trs{m}{d}(:,1))));
            timedecode_stats_wtf{m}(d) = mean(timedecode_stats_trs{m}{d}(:,1));
        end
        timedecode_kfl{tw}{m}(d) = mean(abs(allcurtax-tmp2));
        mdl2 = fitlm(allcursigs,allcurtax);
        timedecode_coefs{tw}{m}{d} = mdl2.Coefficients.Estimate(2:end);
        fprintf("%s %d\n",mice.name(m),d);
    end
end
end
%% Fig 5E,F, S4A
[rsq_ep,kfl_ep,rsq_a,eps_a,kfl_a] = deal(cell(ntw,1));
for tw = 1:ntw
rsq_ep{tw} = cell(1,4);
if tw==1, rsq_ep_t = cell(1,4); end
kfl_ep{tw} = cell(1,4);
neps = length(epoch_list);
eplbls = ["Day 1","2-3","4-6","7-11"];
rsq_a{tw} = []; kfl_a{tw} = [];
eps_a{tw} = [];
for  m = 1:nmice
    for d = 1:mice.ndays(m)
        curep = find(mice.days{m}.epoch(d)==epoch_list)+1;
        if mice.days{m}.trueDay(d)==min(mice.days{m}.trueDay), curep = 1; end
        rsq_ep{tw}{curep} = cat(1,rsq_ep{tw}{curep},timedecode_stats{tw}{m}(d));
        if tw==1, 
            rsq_ep_t{curep} = cat(1,rsq_ep_t{curep},timedecode_stats_wtf{m}(d).^2); 
            assert(~isnan(timedecode_stats_wtf{m}(d).^2));
        end
        rsq_a{tw} = [rsq_a{tw} timedecode_stats{tw}{m}(d)];
        kfl_ep{tw}{curep} = cat(1,kfl_ep{tw}{curep},timedecode_kfl{tw}{m}(d));
        kfl_a{tw} = [kfl_a{tw} timedecode_kfl{tw}{m}(d)];
        eps_a{tw} = [eps_a{tw} curep];
    end
end
barwitherr_stde(rsq_ep{tw},dots=true);
arrayfun(@(x)fprintf("Epoch %d, %d sessions\n",x,length(rsq_ep{tw}{x})),1:4);
xticklabels(eplbls);
axis tight;
xlim([0.5 4.5]);
if tw==1, yl = ylim; yl(1) = -0.01; ylim(yl);
else, ylim(yl);
end
[~,~,s]=anova_disp(rsq_a{tw},y=eps_a{tw});
title("twin: ["+join(string(twins_all{tw})," ")+"]s, "+s)
ylabel("GrC time decoding accuracy (R^2)")
end
%% Fig 5BC
m = 5; d1  = 1; d2=7;
tw=1;
tmp = timedecode_output{tw}{m}{d1};
r2_1 = timedecode_stats{tw}{m}(d1);
tmp_2 = timedecode_output{tw}{m}{d2};
r2_2 = timedecode_stats{tw}{m}(d2);
nrs = min(size(tmp,2),size(tmp_2,2));
fprintf("%d trials\n",nrs)
fprintf("Day 1: %d, Day 2: %d cells\n",mice.data{m}{d1}.nIC_GrC,mice.data{m}{d2}.nIC_GrC)
ix_samp = randsample(size(tmp,2),nrs);
ix_samp_2 = randsample(size(tmp_2,2),nrs);
figure(pos=[1159          22         173         515]);
ax1=subaxis(2,1,1,'ML',0.15); plot(twin_tmplt{tw}',tmp(:,ix_samp),'color',[0.3 0.3 0.3])
hold on; plot(twin_tmplt{tw}',twin_tmplt{tw}','k'); axis tight;
text(-0.8,0,sprintf("R^2=%.3g",r2_1));
box off
ax2=subaxis(2,1,2); 
plot(twin_tmplt{tw}',tmp_2(:,ix_samp_2),'color',[0.3 0.3 0.3])
hold on; plot(twin_tmplt{tw}',twin_tmplt{tw}','k'); axis tight;
text(-0.8,0,sprintf("R^2=%.3g",r2_2));
box off
ylabel("GrC reward delay decoded time (s)")
linkaxes([ax1,ax2])
timax_str(5)
title(sprintf("mouse:%s d1:%d,%s d2:%d,%s; %d trials",mice.name(m),d1,mice.days{m}.date(d1),d2,mice.days{m}.date(d2),nrs))
%% Fig 5D, S4A
outputs_ep = cell(ntw,1);
for tw=[1 2]
twins = twins_all{tw};
herrs = [];
outputs_ep{tw} = cell(1,4);
for m = 1:nmice
    ccplt = [0.4667    0.6745    0.1882; 0.1020    0.1020    0.5020;...
         0.2510    0.2510    0.7490; 0.4941    0.1843    0.5569];
    for d = 1:mice.ndays(m)
    curep = find(mice.days{m}.epoch(d)==epoch_list)+1;
    if mice.days{m}.trueDay(d)==min(mice.days{m}.trueDay), curep = 1; end
    tmp = timedecode_output{tw}{m}{d}';
    tmp = tmp./mean(range(tmp,2));
    tmp = (tmp-mean(min(tmp,[],2)))*diff(twins)+min(twins);
    outputs_ep{tw}{curep} = cat(1,outputs_ep{tw}{curep},tmp);
    end
end
ccplt = [0.4667    0.6745    0.1882; 0.1020    0.1020    0.5020;...
     0.2510    0.2510    0.7490; 0.4941    0.1843    0.5569];
figure(pos=[1084         115         305         280]);
hlines = [];
for ep = 1:4
[h1,h2]=errorbar_shadeSEM(twin_tmplt{tw}',outputs_ep{tw}{ep}',ccplt(ep,:));
hlines = cat(1,hlines,h1);
hold on;
end
axis tight;
box off;
ylabel("twin: ["+join(string(twins_all{tw})," ")+"]s, GrC reward delay decoded time (s)")
timax_str(5)
title(cat(1,arrayfun(@(x)sprintf("%s, %d trials",eplist(x),length(outputs_ep{tw}{x})),1:4)))
legend(hlines,eplbls,box=0)
box off;
end
%% LTD computation
% Fig 5I,J, S4-J, S5D-F
n34all = [];
LTPtally = [];
for krun = 1:4
eplist = ["day 1","novice","mid","expert"];
USE_RANDOM_WEIGHTS = false;
RANDOM_GRC_PERMUTATION = false;
USE_UNIFORM_WEIGHTS = false;
ELEGWIN_S = [-0.15 -0.025];
REW_CF_WIN = [0 0.25];
USE_LTP = false;

% Compare readouts
switch krun
    case 2, USE_RANDOM_WEIGHTS = true;
    case 3, RANDOM_GRC_PERMUTATION = true;
    case 4, USE_UNIFORM_WEIGHTS = true;
end

% Test Plasticity windows
% switch krun
%     case 1, ELEGWIN_S = [-0.15 -0.025];
%     case 2, ELEGWIN_S = [-0.075 -0.075];
%     case 3, ELEGWIN_S = [-0.08 -0.025];
%     case 4, ELEGWIN_S = [0.075 0.2];
% end

% Test CF inclusion windows
% switch krun
%     case 2, REW_CF_WIN = [0 2];
%     case 3, REW_CF_WIN = [-2 2];
%     case 4, USE_UNIFORM_WEIGHTS = true;
% end

% Compare LTD to LTD+LTP
% switch krun
%     case 2, USE_LTP = true;
% end

USE_TR_SUBS = false;
DO_CSUBSETS = false;
PLOT_EXAMPLE_TRIALS = false;
SHOW_GRC_SORTING = false;
SHOW_WEIGHT_VECS = false;
MEAN_SUBTRACT = false;
LTP_CF_WIN = [-2 2];
LTP_CF_thresh = 0.5;
LTPscaleFac = 1;

grcsiguse_ltd = "sigFilt_GrC";
grcsiguse_decode = "sigFilt_GrC";
grcsiguse_ltd = "sigFilt_GrC";
reslts = struct([]);
cnt = 1;
cnt_ep = [1 1 1 1];
twinant_s = [-1.1015 0.0];
twin_sort_s = [-2 2];
brkptspct = [0 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.8 1];
rng(0);
for m = 1:nmice
    for d = 1:mice.ndays(m)
        curd = mice.data{m}{d};
        tmpxCb = curd.tmpxCb;
        dtimCb = curd.dtimCb;
        curep = find(eplist==mice.days{m}.epoch(d));
        if ~USE_TR_SUBS || curep==4
        if mice.days{m}.trueDay(d)==min(mice.days{m}.trueDay), curep = 1; end
        [twin,twin_s] = winfn([-3 2],tmpxCb);
        [twinall_tmplt,twinall_tmplt_s] = winfn([-3 2],tmpxCb_tmplt);
        [twinLTP_tmplt,twinLTP_tmplt_s] = winfn(LTP_CF_WIN,tmpxCb_tmplt);
        if nnz(curd.goodmvdir)/length(curd.goodmvdir)<0.5, fprintf("%s %d %s %% good: %.3g\n",mice.name(m),d,mice.days{m}.date(d),nnz(curd.goodmvdir)/length(curd.goodmvdir)*100); pause; end
        tr_rew = find(curd.rewarded & curd.mvlen>7 & curd.goodmvdir); %rewarded trials
        tr_omit = find(~curd.rewarded & curd.mvlen>7 & curd.goodmvdir); %omitted reward trials
        tr_fail = find(curd.mvlen<6);
        tr_use_decode = cat(1,tr_rew,tr_omit);
        tr_use_ltd = tr_rew;
        if numel(tr_use_decode)>1
        if USE_TR_SUBS
        ntruse = min(tri_use(tri)*trincrmnt,length(tr_use_ltd));
        tr_use_ltd = randsample(tr_use_ltd,ntruse);
        end
        curspcf = curd.rewAlgn.sp_CF(tr_use_ltd,:,twin);
        curspcf_travg = squeeze(mean(curspcf,1));
        curspcf_r = curd.rewAlgn.sp_CF(tr_rew,:,twin);
        curspcf_travg_r = squeeze(mean(curspcf_r,1));

        tmpwinpkj1 = winfn([-0.3 -0.025],tmpxCb(twin));
        tmpwinpkj2 = winfn([0 0.25],tmpxCb(twin));
        tmprewwin = winfn(REW_CF_WIN,tmpxCb(twin));
        pkjuse = find((mean(curspcf_travg_r(:,tmpwinpkj2),2)-mean(curspcf_travg(:,tmpwinpkj1),2))>0.0);
        if isempty(pkjuse), fprintf("no pkjuse %s %d ntr_rew: %d\n",mice.name(m),d,length(tr_rew)); end
        if ~isempty(pkjuse)
        curspcf = curspcf(:,pkjuse,:);
        curspcf_travg = curspcf_travg(pkjuse,:);

        reslts(cnt).pkjuse = [length(pkjuse),curd.nIC_CF];
        reslts(cnt).name = mice.name(m);
        reslts(cnt).day = d;
        reslts(cnt).ep = curep;
        reslts(cnt).travg_spcf = curspcf_travg;
        eleg_win = winfnz(ELEGWIN_S,dtimCb);
        test2 = false(size(curspcf));
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
        curgrc = curd.rewAlgn.(grcsiguse_ltd)(tr_use_ltd,:,twin);
        curgrc_omit = curd.rewAlgn.(grcsiguse_ltd)(tr_omit,:,twin);
        if RANDOM_GRC_PERMUTATION
        for i = 1:size(curgrc,1)
            for k = 1:size(curgrc,2)
                curgrc(i,k,:) = curgrc(i,k,randperm(size(curgrc,3)));
            end
        end
        end
        curLTD = bsxfun(@times,curgrc,test4);
        curLTD(curLTD<0) = 0;

        if USE_LTP
            curcfsig = double(curd.rewAlgn.sp_CF);
            tmpspvec = reshape(permute(curcfsig,[3 1 2]),[],size(curd.rewAlgn.sp_CF,2));
            spkern = round(2 / dtimCb);
            spkern_s = spkern*dtimCb;
            tmpspvecf = filtfilt(ones(spkern,1)/spkern/spkern_s,1,tmpspvec);
            tmpspvecf = zscore(tmpspvecf,[],1);
            tmpspvecf = reshape(tmpspvecf,size(curcfsig,3),size(curcfsig,1),size(curcfsig,2));
            tmpspvecf = permute(tmpspvecf,[2 3 1]);
            trgrps_LTP = {tr_use_ltd,tr_omit};
            timepoints_LTP = cell(1,2);
            junki=[];
            for trg = 1:2
            tmpspvecf2 = tmpspvecf(trgrps_LTP{trg},pkjuse,twin);
            timepoints_LTP{trg} = false(size(tmpspvecf2));
            ltpwincur = winfn(LTP_CF_WIN,tmpxCb(twin));
            for c = 1:size(tmpspvecf2,2)
                curth = prctile(reshape(tmpspvecf2(:,c,ltpwincur),1,[]),LTP_CF_thresh);
            for tr = 1:size(tmpspvecf2,1) 
                cursp = ltpwincur(find(squeeze(tmpspvecf2(tr,c,ltpwincur))<=curth))';
                cur_LTP = false(length(tmpspvecf2(tr,c,:)),1);
                if ~isempty(cursp)
                tmpix = reshape(cursp+eleg_win,1,[]);
                tmpix(tmpix<1 | tmpix>length(cur_LTP)) = [];
                cur_LTP(tmpix) = true;
                end
                timepoints_LTP{trg}(tr,c,:) = cur_LTP;
            end
            end
            timepoints_LTP{trg} = permute(timepoints_LTP{trg},[1 4 3 2]);
            if trg==1
            assert(nnz(timepoints_LTP{trg}&test4)<50)
            curLTP = bsxfun(@times,curgrc,timepoints_LTP{trg});
            curLTP(curLTP<0) = 0;
            curLTD = curLTD - LTPscaleFac*curLTP;
            n1=squeeze(sum(sum(test4,1),3));
            n2=squeeze(sum(sum(timepoints_LTP{trg},1),3));
            n3 = mean((n2-n1)./n1);
            n4 = median((n2-n1)./n1);
            n34all = cat(1,n34all,[n3 n4]);
            else
                curLTP = bsxfun(@times,curgrc_omit,timepoints_LTP{trg});
                curLTP(curLTP<0) = 0;
                curLTD = cat(1,curLTD,LTPscaleFac*curLTP);
            end
            junk = squeeze(sum(sum(sum(curLTP,1),2),4));
            junki = cat(2,junki,interp1(twin_s,junk,twinLTP_tmplt_s,"linear","extrap"));
            end
            LTPtally = cat(1,LTPtally,junki./sum(junki(:)));
        end
        tmpgrc = reshape(permute(curgrc,[3 1 2]),[],curd.nIC_GrC);
        logitrate = prctile(tmpgrc,95,1);
        curLTD_scaled = 1./(1+exp(-1./logitrate.*curLTD)); 

        curLTD = curLTD_scaled;
        
        curLTD_travg = squeeze(mean(curLTD,1));
        curLTD_trtimeavg = squeeze(mean(curLTD_travg,2));
        curLTD_trtimeavg = (curLTD_trtimeavg-mean(curLTD_trtimeavg,1))./sum(curLTD_trtimeavg,1);
        curLTD_trtimePkjavg = squeeze(mean(curLTD_trtimeavg,2));
        [~,LTDsort] = sort(-curLTD_trtimePkjavg);
        twin_sort = winfn(twin_sort_s,tmpxCb);
        twin_sort_tmplt = winfn(twin_sort_s,tmpxCb_tmplt);
        curgrc_tmpsort = curd.rewAlgn.(grcsiguse_ltd)(tr_use_ltd,:,twin_sort);
        twin_mid = winfn(tmpxCb,[0 4]);
        twin_mid_tmplt = winfn(tmpxCb_tmplt,[0 4]);
        curgrc_tmpsort_m = curd.midAlgn.(grcsiguse_ltd)(tr_use_ltd,:,twin_mid);
        curgrcavg = squeeze(mean(curgrc_tmpsort,1));
        curgrcavg_m = squeeze(mean(curgrc_tmpsort_m,1));
        curgrcinterp = interp1(tmpxCb(twin_sort),curgrcavg',tmpxCb_tmplt(twin_sort_tmplt),"linear","extrap")';
        curgrcinterp_m = interp1(tmpxCb(twin_mid),curgrcavg_m',tmpxCb_tmplt(twin_mid_tmplt),"linear","extrap")';
        brkpts = round(curd.nIC_GrC*brkptspct);
        for k = 1:(length(brkpts)-1)
            curix = [(brkpts(k)+1) brkpts(k+1)];
            reslts(cnt).grcLTDquartile{k} = curgrcinterp(LTDsort(curix(1):curix(2)),:);
            reslts(cnt).grcLTDquartile_m{k} = curgrcinterp_m(LTDsort(curix(1):curix(2)),:);
        end
        if USE_RANDOM_WEIGHTS
            curLTD_trtimePkjavg = shufflev(curLTD_trtimePkjavg);
        end
        if USE_UNIFORM_WEIGHTS
            curLTD_trtimePkjavg = 1./length(curLTD_trtimePkjavg)*ones(size(curLTD_trtimePkjavg));
        else
            curLTD_trtimePkjavg = -curLTD_trtimePkjavg;
        end
        reslts(cnt).LTD_trs_time_Pkj_avg = curLTD_trtimePkjavg;
        twinant = winfn(twinant_s,tmpxCb);
        tmpgrc2 = squeeze(mean(curd.rewAlgn.(grcsiguse_decode)(tr_use_decode,:,twinant),1));
        finantix = zeros(size(tmpgrc2,1),1);
        for c = 1:size(tmpgrc2,1)
            tmp = tmpgrc2(c,:);
            tmp = tmp-min(tmp);
            junk = tmp./sum(tmp);
            finantix(c)=sum(junk.*tmpxCb(twinant));
        end
        reslts(cnt).anticOfftimeLTD_corr = corr(-curLTD_trtimePkjavg,finantix);
        if isnan(reslts(cnt).anticOfftimeLTD_corr)
             reslts(cnt).anticOfftimeLTD_corr = 0;
        end
        [reslts(cnt).regwveccorr,p] = corr(curLTD_trtimePkjavg,-timedecode_coefs{1}{m}{d});
        if SHOW_GRC_SORTING && reslts(cnt).anticOfftimeLTD_corr>0.5 ... % mice.name(m)=="590" && d==8
                && ~USE_RANDOM_WEIGHTS && ~RANDOM_GRC_PERMUTATION && ~USE_UNIFORM_WEIGHTS
            pltsuse = finantix;
            figure(pos=[1203 63 700 604]);
            subaxis(2,2,1,1,1,2,'mr',0.15);
            LTDnormplot = -curLTD_trtimePkjavg;
            LTDnormplot = (LTDnormplot-min(LTDnormplot))./range(LTDnormplot);
            plot(pltsuse,LTDnormplot,'.');
            b = regress(LTDnormplot,[pltsuse ones(size(pltsuse))]);
            axis tight; xl = xlim; yl = ylim; xlim([xl(1)-0.05,xl(2)+0.05]);
            ylim([yl(1)-0.01*range(yl),yl(2)+0.01*range(yl)]);
            hold on; plot(xlim,(xlim*b(1))+b(2),'k')
            text(-0.95,0.85*range(yl)+yl(1),sprintf("r = %.3g\nn=%d",reslts(cnt).anticOfftimeLTD_corr,length(finantix)));
            xlim([-1.05 0]);
            ax = gca;
            tmptitstr = sprintf("%d %d %s %s",m,d,mice.name(m),mice.days{m}.date(d));
            title(tmptitstr)
            box off;
            subaxis(2,2,2,1,'mr',0.15);
            curcfsig = double(curd.rewAlgn.sp_CF);
            tmpspvec = reshape(permute(curcfsig,[3 1 2]),[],size(curd.rewAlgn.sp_CF,2));
            spkern = round(0.2 / dtimCb);
            spkern_s = spkern*dtimCb;
            tmpspvecf = filtfilt(ones(spkern,1)/spkern/spkern_s,1,tmpspvec);
            tmpspvecf = zscore(tmpspvecf,[],1);
            tmpspvecf = reshape(tmpspvecf,size(curcfsig,3),size(curcfsig,1),size(curcfsig,2));
            curcfsig = permute(tmpspvecf,[2 3 1]);
            curspcf_travg = squeeze(mean(curcfsig(tr_use_decode,pkjuse,twin),1));
            cclim = [prctile(curspcf_travg(:),0.1),prctile(curspcf_travg(:),99)];
            % imagesc(tmpxCb(twin),[],curspcf_travg,cclim)
            % h=gca; despos = h.Position; colorbar; h.Position=despos; 
            errorbar_shadeSEM(tmpxCb(twin),curspcf_travg',[0 0 0]); axis tight;
            title(sprintf("%d CFs, %d GrCs, %d trials",length(pkjuse),curd.nIC_GrC,length(tr_use_ltd)));
            subaxis(2,2,2,2,'mr',0.15);
            tmpgrc_plt = squeeze(mean(curd.rewAlgn.(grcsiguse_decode)(tr_rew,:,twin),1));
            tmpgrc_plt = tmpgrc_plt-prctile(tmpgrc_plt,1,2); tmpgrc_plt = tmpgrc_plt./prctile(tmpgrc_plt,95,2);
            cuts = [0.05 0.85];
            imagesc(tmpxCb(twin),[],tmpgrc_plt(flipud(LTDsort),:),cuts)
            h=gca; despos = h.Position; colorbar; h.Position=despos; 
            figureformat_forsaving
            savedir="C:\Users\wagnermj\OneDrive - National Institutes of Health\Lab share\Data\L5-GrC-CF\cohortfigs\science\sim\ex sortings normalized\figures";
            notdone=1;
            while notdone
                s = char(input("save? ","s"));
                if isempty(s), notdone=0;
                elseif s(1)=='y'
                    saveas(gcf,fullfile(savedir,tmptitstr),"svg"); notdone=0;
                else
                    try
                    subaxis(2,2,2,str2double(s(1)));
                    clix = str2double(s(2));
                    cl = clim;
                    cl(clix)=eval("cl(clix)"+s(3)+"0.1*range(cl)");
                    clim(cl);
                    end
                end
            end
            close;
        end
        if ~USE_RANDOM_WEIGHTS && ~RANDOM_GRC_PERMUTATION && ~USE_UNIFORM_WEIGHTS
        mice.data{m}{d}.regwveccorr = reslts(cnt).regwveccorr;
        end
        reslts(cnt).ngrc = length(curLTD_trtimePkjavg);
        if SHOW_WEIGHT_VECS && reslts(cnt).regwveccorr>0.5 && curep==4 ...
            && ~USE_RANDOM_WEIGHTS && ~RANDOM_GRC_PERMUTATION && ~USE_UNIFORM_WEIGHTS
            figure(pos=[716 328 851 491]);
            subaxis(1,2,1);
            junk2 = double(curd.ICmat_GrC>0); junk2 = junk2./sum(sum(junk2,1),2);
            junk3 = curLTD_trtimePkjavg; junk3 = junk3-min(junk3);
            allmasks = sum(bsxfun(@times,junk2,permute(junk3,[3 2 1])),3);
            cm = hot(200); cm = flipud(cm);
            imagesc(allmasks); axis image; colormap(cm);
            clim([prctile(allmasks(allmasks~=0),1),prctile(allmasks(allmasks~=0),99)])
            xticks([]); yticks([]);
            subaxis(1,2,2);
            junk = timedecode_coefs{m}{d};
            junk = ((junk-min(junk(:)))./range(junk(:)))*2;
            allmasks2 = sum(bsxfun(@times,junk2,permute(junk,[3 2 1])),3);
            imagesc(allmasks2); axis image; colormap(cm);
            clim([prctile(allmasks2(allmasks2~=0),1),prctile(allmasks2(allmasks2~=0),99)])
            xticks([]); yticks([]);
            h=gca; despos = h.Position; colorbar; h.Position=despos; 
            title(sprintf("r = %.3g, p = %.3g, %d cells, %d trials, %s %d %d",reslts(cnt).regwveccorr,p,length(curLTD_trtimePkjavg),length(tr_use_decode),mice.name(m),m,d))
            tmptitstr = sprintf("%d %d %s %s",m,d,mice.name(m),mice.days{m}.date(d));
            savedir="C:\Users\wagnermj\OneDrive - National Institutes of Health\Lab share\Data\L5-GrC-CF\cohortfigs\science\sim\ex weight corrs";
            notdone=1;
            while notdone
                s = char(input("save? ","s"));
                if isempty(s), notdone=0;
                elseif s(1)=='y'
                    saveas(gcf,fullfile(savedir,tmptitstr),"svg"); notdone=0;
                else
                    try
                    subaxis(1,2,str2double(s(1)));
                    clix = str2double(s(2));
                    cl = clim;
                    cl(clix)=eval("cl(clix)"+s(3)+"0.1*range(cl)");
                    clim(cl);
                    end
                end
            end
            close;
        end
        
        tmpgrc = curd.rewAlgn.(grcsiguse_decode)(tr_use_decode,:,twinant);
        curkern = round(0.07/dtimCb);
        tmpgrcf = filtfilt(curkern,1,reshape(permute(double(tmpgrc),[3 1 2]),[],size(tmpgrc,2)));
        tmpgrcf = permute(reshape(tmpgrcf,size(tmpgrc,3),size(tmpgrc,1),size(tmpgrc,2)),[2 3 1]);
        tmpgrc= tmpgrcf;
     
        tmpSums = bsxfun(@times,tmpgrc,permute(curLTD_trtimePkjavg,[3 1 2]));
        if MEAN_SUBTRACT
        tmpSums = tmpSums - mean(tmpgrc,2);
        end

        tmpPkj = squeeze(sum(tmpSums,2));
        tmpmn = mean(tmpPkj,2);
        tmpmn = mean(tmpmn);
        tmpPkj_scaled = tmpPkj-tmpmn; 
        bla = range(tmpPkj_scaled,2);
        junkix = find(bla==0);
        tmpPkj(junkix,:) = []; tmpPkj_scaled(junkix,:) = []; bla(junkix) = [];
        bla = mean(bla);
        tmpPkj_scaled = tmpPkj_scaled./bla + 0.5;
        tmpPkj_travg = mean(tmpPkj,1);
        reslts(cnt).Pkj_singlTrs = tmpPkj_scaled;
        twinant_tmplt = winfn(twinant_s,tmpxCb_tmplt);
        tmpPkj_interp = interp1(tmpxCb(twinant),tmpPkj_scaled',tmpxCb_tmplt(twinant_tmplt),"linear","extrap")';
        reslts(cnt).Pkj_interp = tmpPkj_interp;
        reslts(cnt).sID = repmat([cnt_ep(curep) curep m d],size(tmpPkj_interp,1),1);
        tmpTime = repmat(-(tmpxCb(twinant)),[size(tmpPkj,1) 1]);
        if ~USE_RANDOM_WEIGHTS && ~RANDOM_GRC_PERMUTATION && ~USE_UNIFORM_WEIGHTS
        mice.data{m}{d}.pkjtimeacc_trs = arrayfun(@(x)corr(tmpTime(x,:)',tmpPkj_scaled(x,:)'),1:size(tmpTime,1));
        mice.data{m}{d}.pkjtime_r2_trs = mice.data{m}{d}.pkjtimeacc_trs.^2;
        mice.data{m}{d}.pkjtime_trs_ix = tr_use_decode;
        end
        tmpTime = reshape(tmpTime',[],1);
        tmpPkj_rs = reshape(tmpPkj_scaled',[],1);
        if isnan(corr(tmpTime,tmpPkj_rs)), pause; end
        reslts(cnt).PkjTime_corr_stScale = single(corr(tmpTime,tmpPkj_rs));
        reslts(cnt).PkjTime_r2_stScale = corr(tmpTime,tmpPkj_rs).^2;
        if ~USE_RANDOM_WEIGHTS && ~RANDOM_GRC_PERMUTATION && ~USE_UNIFORM_WEIGHTS
        mice.data{m}{d}.pkjtimeacc = corr(tmpTime,tmpPkj_rs);
        mice.data{m}{d}.pkjtime_r2 = corr(tmpTime,tmpPkj_rs).^2;
        end
        reslts(cnt).PkjTime_corr_pctOfOptimal = ...
            min(reslts(cnt).PkjTime_corr_stScale^2/timedecode_stats{1}{m}(d),2);
        mice.data{m}{d}.pctopt = min(reslts(cnt).PkjTime_corr_stScale^2/timedecode_stats{1}{m}(d),1);
        reslts(cnt).lineregreslts = sqrt(timedecode_stats{1}{m}(d));
        fprintf("%s %d Pkj-time corr %.2g, %% of optimal %.2g\n",...
            mice.name(m),d,reslts(cnt).PkjTime_corr_stScale,reslts(cnt).PkjTime_corr_pctOfOptimal);
        cnt = cnt+1;
        cnt_ep(curep) = cnt_ep(curep)+1;
        end
        end
        end
    end
end
% Collating LTD results into a struct
clear L;
L.pkjuse = sum(cat(1,reslts.pkjuse));
L.ns = length(reslts);
fprintf("nsessions: %d, Npkj: %d, Npkjtot: %d\n",L.ns,L.pkjuse(1),L.pkjuse(2));
L.eps = cat(1,reslts.ep);
L.corrs_offtime = arrayfun(@(x)cat(1,reslts(L.eps==x).anticOfftimeLTD_corr),1:4,unif=false);
L.corrs_offtime_a = cat(1,reslts.anticOfftimeLTD_corr);
L.PkjTime_corr_stScale = arrayfun(@(x)cat(1,reslts(L.eps==x).PkjTime_corr_stScale),1:4,unif=false);
if isfield(reslts,"PkjTime_corr_cSubsets"),
    L.PkjTime_corr_cSubs = arrayfun(@(x)cat(1,reslts(L.eps==x).PkjTime_corr_cSubsets),1:4,unif=false);
end
L.PkjTime_corr_stScale_a = cat(1,reslts.PkjTime_corr_stScale);
L.linregaccuracy = cat(1,reslts.lineregreslts);
L.regwveccorr_e = arrayfun(@(x)cat(1,reslts(L.eps==x).regwveccorr),1:4,unif=false);
L.regwveccorr = cat(1,reslts.regwveccorr);
L.ngrcs = cat(1,reslts.ngrc);
L.linregaccuracy_e = arrayfun(@(x)cat(1,reslts(L.eps==x).lineregreslts),1:4,unif=false);
L.pctOfLinReg = arrayfun(@(x)100*cat(1,reslts(L.eps==x).PkjTime_corr_pctOfOptimal),1:4,unif=false);
L.pctOfLinReg_a = 100*cat(1,reslts.PkjTime_corr_pctOfOptimal);
L.PkjTrs = arrayfun(@(x)cat(1,reslts(L.eps==x).Pkj_interp),1:4,unif=false);
L.sID_Trs = arrayfun(@(x)cat(1,reslts(L.eps==x).sID),1:4,unif=false);
switch krun
    case 1, L_true = L; reslts_true = reslts;
    case 2, L_rand = L; reslts_rand = reslts;
    case 3, L_permuted = L;
    case 4, L_unif = L;
end
disp(krun);
end
% end
%% Fig. 5K, S5I
yl = [-0.03    0.65];
yl2 = [-0.03    0.65];
yls = {yl,yl2,yl,yl2,[-0.025 0.48],[-0.06 0.52]};

labls = ["peak time","peak time"];
labls = arrayfun(@(k)sprintf("Corr bw GrC %s and LTD",labls{k}),1:length(labls));
labls(end+1) = "Corr between GrC weights from CF-LTD vs regression (r)";
labls(end+1) = "Pkj time decoding (r)";
barwitherr_stde(L_true.corrs_offtime,dots=true); 
barwitherr_stde(L_permuted.corrs_offtime,fhandle=gcf,color=[0.7 0.7 0.7]);
% ylim(yls{k})
ylabel(labls{1})
[f,p,s]=anova_disp(L_true.corrs_offtime);
title(s);
xticklabels(["Day 1","2-3","4-6","7-11"]);

ccplt = [0.4667    0.6745    0.1882; 0.1020    0.1020    0.5020;...
     0.2510    0.2510    0.7490; 0.4941    0.1843    0.5569];
figure; colororder(ccplt);
h=gscatter(L_true.linregaccuracy.^2,L_true.PkjTime_corr_stScale_a.^2,L_true.eps); 
xl = xlim; yl = ylim; 
box off
b = regress(L_true.PkjTime_corr_stScale_a.^2,[L_true.linregaccuracy.^2 ones(size(L_true.linregaccuracy))]);
hold on; plot([0 0.8],[0 0.8]*b(1)+b(2),'k'); xlim(xl); ylim(yl);
hold on; plot([0 0.7],[0 0.7],'k--'); xlim(xl); ylim(yl);
axis equal square
fprintf("reg acc vs pkj acc (r): ");
[r,p]=corr_disp(L_true.linregaccuracy.^2,L_true.PkjTime_corr_stScale_a.^2);
text(0.1,0.5,sprintf("r = %.3g\np = %.3g\nn = %d",r,p,length(L_true.PkjTime_corr_stScale_a)))
legend(["Day 1","2-3","4-6","7-11"])
xlabel("CF-LTD timing accuracy (R2)")
ylabel("Optimal decodingaccuracy (R2)")
%% Fig 5L; S4K,L GrC avg profiles grouped by LTD percentile
twinquartileplts={[0 2.1],[0 3]};
normheight = [true false];
nplt = length(twinquartileplts);
for p = 1:nplt
twinquartileplt = twinquartileplts{p};

twincur = twin_mid_tmplt;
ixuse = winfn(twinquartileplt,tmpxCb_tmplt(twincur));

figure(pos=[968   198   915   259]);
cctmp = cool(length(brkptspct));
hlines = [];
herrs = [];
ixmx = [];
nstr="";
risefall=[];
for k = 1:(length(brkptspct)-1)
    allqtile = arrayfun(@(x)cat(1,x.grcLTDquartile_m{k}),reslts_true(ismember(cat(1,reslts_true.ep),[3 4])),unif=false);
    disp(length(reslts_true(ismember(cat(1,reslts_true.ep),[3 4]))))
    allqtile = cat(1,allqtile{:});
    if any(isnan(allqtile(:))), pause; end
    allqtile = allqtile(:,ixuse);
    
    tp = tmpxCb_tmplt(twincur(ixuse));
    aqm_o = mean(allqtile,1);
    aqm = (aqm_o-min(aqm_o))/range(aqm_o);
    aqm_a = (allqtile-min(aqm_o))/range(aqm_o);
    
    if normheight(p)
    hlines(end+1)=plot(tp,aqm,col=cctmp(k,:));
    hold on;
    xlim(tmpxCb_tmplt(twincur([ixuse(1) ixuse(end)])))
    ylim([-0.05 1.05])
    else
    xt = tmpxCb_tmplt(twincur(ixuse))';
    yt = allqtile';
    hlines(end+1) = errorbar_shadeSEM(xt,yt,cctmp(k,:));
    hold on;
    end
    nstr=nstr+" "+num2str(brkptspct(k+1))+":"+num2str(size(allqtile,1));
    hold on;
    fprintf("%1.2g-%1.2g, %d\n",brkptspct(k),brkptspct(k+1),size(allqtile,1))
end
uistack(hlines,"top")
disp(nstr)
xlabel("Time relative to movement (s)");
figureformat_forsaving
if ~normheight(p)
    legend(hlines,string(100*brkptspct(1:(end-1)))+"–"+string(100*brkptspct(2:end)))
    ylabel("Fluorescence (zsc)")
    ylim([-0.3584    2.7024]);
else
    ylabel("Normalized signal")
end
plot(meanmvtime_all*[1 1],ylim,'k');
xlim(twinquartileplt)
figureformat_forsaving
end
%% Fig 6D,E, S5G
% Comparing true, perumted-GrC LTD, and random weight stats
toplt = {L_true.pctOfLinReg{4},single(L_rand.pctOfLinReg{4}),single(L_unif.pctOfLinReg{4}),single(L_permuted.pctOfLinReg{4})};
toplt = toplt([1 2 4 3]);
rsg = {[1,2],[1,3],[1,4]};
barwitherr_stde(toplt,color=[.73 0.26 0.89; .77 .6 .42; .2 .2 .2; .7 .7 .7],...
    dots=true,ranksumgrps=rsg);
ylabel("pkj % of optim"); ylim([-2 100]); yticks([0:20:100]); 
figureformat_forsaving;

varuse = "PkjTime_corr_stScale";
grps = [L_true.(varuse),L_unif.(varuse)(4),L_rand.(varuse)(4)];
tmpcc = [repmat([.73 0.26 0.89],4,1);.2 .2 .2;.77 .6 .42];
ranksumgrps = {[4,5],[4,6]};
barwitherr_stde(grps,color=tmpcc,dots=true,...
    ranksumgrps=ranksumgrps,anovagrp=1:4); 
xlim([0.5 6.5])
ylim tight;
ylabel(labls{4})

grps = [L_true.regwveccorr_e,L_permuted.regwveccorr_e];
tmpcc=[repmat([.73 0.26 0.89],4,1);repmat([.75 .75 .75],4,1)];
barwitherr_stde(grps,color=tmpcc,dots=true,anovagrp=1:4,...
    dotgrp=1:4,xcoords=[1:4 1:4],ranksumgrps={[4,8]}); 
ylabel(labls{3}); yticks(0:.1:0.7)
%% Fig 6B,C, S5A-C; single session results
tmp3 = L_true.PkjTime_corr_stScale{4} - L_unif.PkjTime_corr_stScale{4};
[~,tmp4] = sort(tmp3,'desc');
tmpsid = L_true.sID_Trs{4}(:,1);
for i = 1:20
    curi = tmp4(i);
    curii = find(tmpsid==curi);
    m=L_true.sID_Trs{4}(curii(1),3);
    d=L_true.sID_Trs{4}(curii(1),4);
    curn = min(40,length(curii));
    notdone=true;
    curplt = randsample(curii,curn);
    figure; clear a;
    a(1)=subaxis(1,3,1);
    plot(tmpxCb_tmplt(twinant_tmplt),L_true.PkjTrs{4}(curplt,:),col=[0.7 0.7 0.7]);
    ylabel(mice.name(m)+" "+mice.days{m}.date(d))
    xlabel("CF-LTD");
    title(sprintf("r=%1.2g",L_true.PkjTime_corr_stScale{4}(curi))+...
        sprintf(" %dgrcs %dtrials",length(curii),mice.data{m}{d}.nIC_GrC))
    axis tight;
    xticks(-1:.2:0); axl = axis;
    hold on; plot(xlim,-xlim,'k'); axis(axl);
    a(2)=subaxis(1,3,3);
    plot(tmpxCb_tmplt(twinant_tmplt),L_unif.PkjTrs{4}(curplt,:),col=[0.7 0.7 0.7]);
    xlabel("Simple avg");
    title(sprintf("r=%1.2g",L_unif.PkjTime_corr_stScale{4}(curi)))
    axis tight;
    xticks(-1:.2:0);axl = axis;
    hold on; plot(xlim,-xlim,'k'); axis(axl);
    a(3)=subaxis(1,3,2);
    plot(tmpxCb_tmplt(twinant_tmplt),-timedecode_output{1}{m}{d}(:,curplt-curii(1)+1)',col=[0.7 0.7 0.7]);
    xlabel("Optimal");
    title(sprintf("r=%1.2g",sqrt(timedecode_stats{1}{m}(d))))
    axis tight;
    xticks(-1:.2:0);axl = axis;
    hold on; plot(xlim,-xlim,'k'); axis(axl);
    linkaxes(a);
    figureformat_forsaving;
    pause
    close;
end
%% Fig 6H; Compare decoding stats to behavior
[behavStat_m, neuralStats_m, behav_neural_corr, all_m_ix] = deal([]);
corrs_m = cell(1,4);
for m = 1:nmice
    for d = 1:mice.ndays(m)
        behavVar = "lickDiff";
        neuralVar = "pkjtime_r2";
        if isfield(mice.data{m}{d},behavVar) && isfield(mice.data{m}{d},neuralVar)
        curbehavstat = mice.data{m}{d}.(behavVar); 
        curd = mice.data{m}{d};
        alltrials = 1:length(curd.rewarded);
        if length(curbehavstat)>1
        alltrials_use = alltrials(~isnan(curbehavstat(alltrials,1))); %trials with valid lick-sensor readings
        [~,tr_use_joint_neur,~] = intersect(mice.data{m}{d}.pkjtime_trs_ix,alltrials_use);
        tr_use_joint = mice.data{m}{d}.pkjtime_trs_ix(tr_use_joint_neur); 
        assert(~any(isnan(curbehavstat(alltrials_use,1))));
        cur_neuralStat = mice.data{m}{d}.(neuralVar); 
        curep = find(mice.days{m}.epoch(d)==epoch_list)+1;
        if mice.days{m}.trueDay(d)==min(mice.days{m}.trueDay), curep = 1; end
        if length(curbehavstat)>1
            if length(cur_neuralStat)>1
            corrs_m{curep} = cat(1,corrs_m{curep},corr(curbehavstat(tr_use_joint),cur_neuralStat(tr_use_joint_neur)'));
            end
            curbehavstat = curbehavstat(alltrials_use,:); 
            curbehavstat = mean(curbehavstat);
            all_m_ix = cat(1,all_m_ix,[m d curep]);
            cur_neuralStat  = mean(cur_neuralStat);
            behavStat_m = cat(1,behavStat_m, curbehavstat);
            neuralStats_m = cat(1,neuralStats_m, cur_neuralStat);
        end
        end
        end
    end
end
[r,p]=corr_disp(behavStat_m,neuralStats_m,'spearman');
cmap1 = cool(4);
ccplt = [0.4667    0.6745    0.1882; 0.1020    0.1020    0.5020;...
     0.2510    0.2510    0.7490; 0.4941    0.1843    0.5569];
figure; colororder(ccplt);
xuse = behavStat_m(:,i);
yuse = neuralStats_m;
h=gscatter(xuse,yuse,all_m_ix(:,3),ccplt); box off
hold on;
xl = xlim; xl(1) = xl(1)-0.05*range(xl);
yl = ylim; yl(1) = yl(1)-0.05*range(yl);
xl(2) = xl(2)+0.05*range(xl); xlim(xl);
yl(2) = yl(2)+0.05*range(yl); ylim(yl);
b = regress(yuse,[xuse ones(length(behavStat_m(:,i)),1)]);
xl = xlim;
plot(xl,xl*b(1)+b(2),'k-')
title(sprintf("col=%d r=%.3g p=%.3g n=%d",i,r,p,length(behavStat_m)))
legend off;
p
figureformat_forsaving;
disp("mice:"+length(unique(all_m_ix(:,1)))+", sess:"+length(unique(all_m_ix(:,1:2),'rows')))
strjoin(arrayfun(@(x)"ep"+x+":"+length(unique(all_m_ix(all_m_ix(:,3)==x,1:2),'rows')),1:4))
%% Fig S5J,K; diff windows
varuse = "PkjTime_corr_stScale";
toplt = {L_true.(varuse){4},single(L_rand.(varuse){4}),single(L_permuted.(varuse){4}),single(L_unif.(varuse){4})};
barwitherr_stde(toplt)
signrank(toplt{1},toplt{4})
ylim([0 0.6])
%% S4O peak time versus distant activity levels
mice2s = load("C:\Users\wagnermj\Desktop\data to share\learning_1s_to_2s_GrC_CF.mat");
mice2s = mice2s.groups;
testtimes = [1.1 2.1];
grcpeaks_activity = zeros(0,3+length(testtimes)+3);
persescorr = zeros(0,2+3);
mdats = {mice.data,mice2s{3}};
for mm = 1:2
    nmice_cur = length(mdats{mm});
for m = 1:nmice_cur
    if mm==1, ndays_cur = mice.ndays(m);
    else, ndays_cur = 1;
    end
    for d = 1:ndays_cur
        if mm==1,
        curep = find(mice.days{m}.epoch(d)==eplist);
        if mice.days{m}.trueDay(d)==1, curep = 1; end
        curd = mice.data{m}{d};
        else, 
            curep=5; 
            curd = mice2s{3}{m};
        end
        trrew = find(curd.mvlen>7 & curd.goodmvdir & curd.rewarded);
        currewtime = mean(curd.rewtimes(trrew)-curd.midpt(trrew))*0.005;
        if mm==1, testtimes(1) = currewtime;
        else, testtimes(1)=1.1; end
        curgrcs = curd.midAlgn.sigFilt_GrC(trrew,:,:);
        if any(isinf(curgrcs(:))), disp([m d]); pause; end
        curgrcm = squeeze(mean(curgrcs,1));
        curgrcm(isinf(curgrcm))=0;
        curgrcmf = butter_filtfilt(2,0.1,curgrcm')';
        lateract= zeros(curd.nIC_GrC,2);
        pktimes= zeros(curd.nIC_GrC,2);
        r = [];
        for k = 1:length(testtimes)
        [windel,windel_s] = winfn([-1 testtimes(k)],curd.tmpxCb);
        curgrcmf = curgrcmf-min(curgrcmf,[],2);
        [pks,pktimes(:,k)] = max(curgrcmf(:,windel),[],2);
        risetimes = zeros(size(pktimes(:,k)));
        for c = 1:curd.nIC_GrC
            risetimes(c) = find(curgrcmf(c,windel)>0.7*pks(c),1);
        end
        if 1, pktimes(:,k)=risetimes; end
        pktimes(:,k) = windel_s(pktimes(:,k));
        [~,ix] = sort(pktimes(:,k));
        % [rnuse,rnuse_s] = winfn([-1 4],curd.tmpxCb);
        rnuse = windel; rnuse_s = windel_s;
        testtimes_ix = arrayfun(@(x)minix(abs(x-rnuse_s)),testtimes(k));
        if 1
        curgrcmnorm = curgrcm(:,rnuse)-min(curgrcm(:,rnuse),[],2);
        assert(all(all(curgrcmnorm>=0)))
        curgrcmnorm  = curgrcmnorm./range(curgrcm(:,rnuse),2);
        assert(all(all(curgrcmnorm <=1)))
        else, curgrcmnorm = curgrcm(:,rnuse);
        end
        lateract(:,k) = curgrcmnorm(:,testtimes_ix);
        r(k) = corr(pktimes(:,k),lateract(:,k));
        end
        if mm==-2
        figure; 
        subaxis(1,3,1);
        imagesc(curd.tmpxCb,[],curgrcm(ix,:),[-0.4 0.6]); 
        hold on; plot([0 0],ylim,'k');
        plot([1.1 1.1],ylim,'k'); 
        subaxis(1,3,2);
        scatter(pktimes(:,1),lateract(:,1))
        title(strvec([m d curep r1(1,2)]));
        subaxis(1,3,3);
        scatter(pktimes(:,2),lateract(:,2))
        title(r2(1,2))
        pause; close;
        end
        grcpeaks_activity = cat(1,grcpeaks_activity,[pktimes,pks,lateract,repmat([m d curep],curd.nIC_GrC,1)]);
        persescorr = cat(1,persescorr,[r m d curep]);
    end
end
end
epixa = arrayfun(@(x)find(grcpeaks_activity(:,end)==x),1:5,unif=0);
epix = arrayfun(@(x)find(persescorr(:,end)==x),1:5,unif=0);
sescorrep = arrayfun(@(x)persescorr((persescorr(:,end)==x),1),1:5,unif=0);
barwitherr_stde([sescorrep(1:4),{persescorr(epix{4},2)},{persescorr(epix{5},2)}],...
    dots=true,signrank=true,ranksumgrps={[4 5],[5,6]},...[4,6],
    color=cat(1,cc_ep,cc_ep(4,:),[0 .5 .5]),anovagrp=1:4);
xl=[-1.1 2.1];
figure; scatter(grcpeaks_activity(epixa{4},1),grcpeaks_activity(epixa{4},4),'filled',...
    'markeredgecolor',[0 0 0],'markerfacecolor',cc_ep(4,:));
xlim(xl);
[r2s,p2s]=corr(grcpeaks_activity(epixa{4},1),grcpeaks_activity(epixa{4},4));
title("1sexp at 1s: "+strvec([r2s p2s length(epixa{4})]))

figure; scatter(grcpeaks_activity(epixa{1},1),grcpeaks_activity(epixa{1},4),'filled',...
    'markeredgecolor',[0 0 0],'markerfacecolor',cc_ep(1,:));
xlim(xl);
[r2s,p2s]=corr(grcpeaks_activity(epixa{1},1),grcpeaks_activity(epixa{1},4));
title("1snov at 1s: "+strvec([r2s p2s length(epixa{1})]))

figure; scatter(grcpeaks_activity(epixa{4},2),grcpeaks_activity(epixa{4},5),'filled',...
    'markeredgecolor',[0 0 0],'markerfacecolor',cc_ep(4,:));
xlim(xl);
[r2s,p2s]=corr(grcpeaks_activity(epixa{4},2),grcpeaks_activity(epixa{4},5));
title("1sexp at 2s: "+strvec([r2s p2s length(epixa{4})]))

figure; scatter(grcpeaks_activity(epixa{5},2),grcpeaks_activity(epixa{5},5),'filled',...
    'markeredgecolor',[0 0 0],'markerfacecolor',[0 .5 .5]);
xlim(xl);
[r2s,p2s]=corr(grcpeaks_activity(epixa{5},2),grcpeaks_activity(epixa{5},5));
title("2sexp at 2s: "+strvec([r2s p2s length(epixa{5})]))
%% S5H global learning curves lick and grc
curtrs = cell(nmice,1);
[epcorrs,epcorrsxc,epcorrlag,epxc] = deal(cell(4,1));
for m = 1:nmice
    curtrs{m} = zeros(0,3);
for d = 1:mice.ndays(m)
    curd = mice.data{m}{d};
    if isfield(curd,'lickDiff')
    curep = find(mice.days{m}.epoch(d)==eplist);
    if mice.days{m}.trueDay(d)==min(mice.days{m}.trueDay), curep = 1; end
    tmplick = curd.lickPref;
    tmpgrc = timedecode_stats_trs{m}{d}(:,1);
    truse = intersect(find(~isnan(tmplick)),timedecode_stats_trs{m}{d}(:,2));
    [~,trused] = ismember(truse,timedecode_stats_trs{m}{d}(:,2));
    
    ntr = length(truse);
    curcurtrs = [tmplick(truse),tmpgrc(trused,1).^2,repmat(mice.days{m}.trueDay(d),ntr,1)];
    curtrs{m} = cat(1,curtrs{m},curcurtrs);
    test = filter(10,1,curcurtrs);
    [xc,lags] = xcorr(test(:,1)-mean(test(:,1)),test(:,2)-mean(test(:,2)),'coef',30);
    [mx,mxix] = min(xc); mxix=lags(mxix);
    epcorrs{curep} = cat(1,epcorrs{curep},corr(test(:,1),test(:,2)));
    epcorrsxc{curep} = cat(1,epcorrsxc{curep},mx);
    epcorrlag{curep} = cat(1,epcorrlag{curep},mxix);
    epxc{curep} = cat(1,epxc{curep},xc');
    end
end
test = curtrs{m};
end
lns = cellfun(@(x)length(x),curtrs);
lnst = find(lns>200);
hasd1 = find(cellfun(@(x)sum(x(:,3)<3,1)>1,curtrs));
hasd1 = intersect(hasd1,find(cellfun(@(x)sum(x(:,3)>=7 | isnan(x(:,3)),1)>1,curtrs)));
muse = intersect(lnst,hasd1);
curtrst = curtrs(muse);
lnmn = mean(lns(muse));
disp(lns(muse))
alllickinterp = arrayfun(@(x)interp1(1:lns(muse(x)),curtrst{x}(:,1),linspace(1,lns(muse(x)),lnmn),"linear","extrap")',1:length(muse),unif=false);
allgrcinterp = arrayfun(@(x)interp1(1:lns(muse(x)),curtrst{x}(:,2),linspace(1,lns(muse(x)),lnmn),"linear","extrap")',1:length(muse),unif=false);
alllickinterp = cat(2,alllickinterp{:});
allgrcinterp = cat(2,allgrcinterp{:});
nint=size(alllickinterp,1);
alllickinterpf = butter_filter(2,0.2,alllickinterp);
allgrcinterpf = butter_filter(2,0.2,allgrcinterp);

figure(pos=[1168         377         400         420]); 
h1=errorbar_shadeSEM(1:nint,alllickinterpf);
hold on; axis tight;
yyaxis right;
h2=errorbar_shadeSEM(1:nint,allgrcinterpf,[0 0 0]);
axis tight;
ax = gca; ax.YGrid=true;
xl= xlim; xl(1)=-10; xlim(xl)
%% LTP analysis
figure; h1=errorbar_shadeSEM(twinLTP_tmplt_s',LTPtally(:,1:length(twinLTP_tmplt))',[0 0 1])
hold on; h2=errorbar_shadeSEM(twinLTP_tmplt_s',LTPtally(:,(length(twinLTP_tmplt)+1):end)',[1 0 0])
xlabel("Time relative to reward (s)"); ylabel("LTP probability");
xlim([-2 1.95])
legend([h1,h2],["Rewarded trials","Omission"])
figureformat_forsaving
figure; histogram_cell({n34all(:,2)})
barwitherr_stde({L_true.pctOfLinReg{4},L_rand.pctOfLinReg{4}},dots=true)