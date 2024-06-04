%% behavior quantifications. 
% Fig2A,B,D; 
% Fig. S2A,B,C,D,E
epochs = ["day 1" "novice" "mid" {"expert","lowmag"}]; nep=4;
eplbls = ["Day 1","2-3","4-6","7-11"];
bquantnms = ["ITI","mvDur","ntrials","attempts_per_min","successes_per_min","success_rate","mvLen","failfrac","ITI_t",...
    "ITI_omit","mvDur_t","anticLick","STlick","STlick_o","STlick_raw",...
    "STlick_raw_o","STrate_prob","STreach","STrate","STrate_all","STrate_o",...
    "omitM_ix","preRewBodyStereotypy_coords","preRewBodyMovement_coords","STbody","STbody_speed","lick_mice"]; 
nbq = length(bquantnms);
filt_lick = true;
filtwin = 60;
algnuse = "rewAlgn";
firstbody=true;
COMPUTE_OMIT_OFFTIMES = false;
behaviorQuantsEpoch = table(Size=[nep,nbq],VariableTypes=repmat("cell",1,nbq));
behaviorQuantsEpoch.Properties.VariableNames = bquantnms;
behaviorQuantsEpoch.Row = epochs(1:4);
licksessioncount = zeros(1,4);
lick_mice  = [];
for m = 1:nmice
    for d = 1:mice.ndays(m)
        curd = mice.data{m}{d};
        tmpITI = [0; diff(curd.truestart)*curd.dtb];
        tmpomit = find(curd.mvlen>7 & ~curd.rewarded); tmpomit(tmpomit>=length(curd.mvlen))=[];
        tmpITI_omit = (curd.truestart(tmpomit+1)-curd.truestart(tmpomit))*curd.dtb;
        mice.data{m}{d}.ITI = tmpITI;
        mice.data{m}{d}.medITI = median(tmpITI);
        curep = mice.days{m}.epoch(d);
        tmpwin2 = winfn([-2 2],curd.tmpxx);
        if licksensor_valid{m}(d)~="n"
            valid_trials = intersect(curd.valid_lick_trials,find(curd.mvlen>7));
            valid_trials_r = intersect(valid_trials,find(curd.rewarded));
            valid_trials_o = intersect(valid_trials,find(~curd.rewarded));
            if licksensor_valid{m}(d)~="n"
            if ismember(algnuse,["startAlg","midAlgn","endAlgn","splitAlgn"]), win_prerew = meanmvtime - [0.3 0.05];
            elseif algnuse=="rewAlgn", win_prerew = [-0.2 0]; 
            end
            win_prerew = winfn(win_prerew,curd.tmpxx);
            pre_reward_win = winfn([-meanmvtime_all 0],curd.tmpxx);
            cur_rawlick = double(curd.(algnuse).lick);
            cur_rawrate = [zeros(size(cur_rawlick,1),1) double(diff(cur_rawlick,[],2)>0)]./curd.dtb;
            cur_filtlick = cur_rawlick;
            probnormsig = ones(size(cur_rawrate,1),1);
            probnormsig(valid_trials) = sum(cur_rawrate(valid_trials,pre_reward_win),2); % per-trial pre-reward sum
            cur_probrate = bsxfun(@times,cur_rawrate,1./probnormsig);
            cur_probrate(find(probnormsig==0),:) = 1./length(pre_reward_win);
            if filt_lick
                cur_filtlick = filtfilt(ones(filtwin,1)/filtwin,1,cur_filtlick')';
                cur_probrate = filter(ones(filtwin,1)/filtwin,1,cur_probrate')';
                cur_filtrate = filter(ones(filtwin,1)/filtwin,1,cur_rawrate')';
            end
            assert(~(any(isnan(cur_filtlick(:))) || any(isinf(cur_filtlick(:)))));
            assert(~(any(isnan(cur_probrate(:))) || any(isinf(cur_probrate(:)))));
            mice.data{m}{d}.anticLick = nan(size(curd.(algnuse).lick,1),1);
            mice.data{m}{d}.anticLick(valid_trials) = mean(cur_filtrate(valid_trials,win_prerew),2);
            mice.data{m}{d}.valid_trials_o = valid_trials_o;

            latewin = winfn([-0.2 0],curd.tmpxx);
            earlywin = winfn([-0.8 -0.6],curd.tmpxx);
            lickuse = cur_filtrate;
            lickuse_prob = cur_probrate;

            mice.data{m}{d}.lickDiff = nan(size(curd.(algnuse).lick,1),1);
            mice.data{m}{d}.lickPref = nan(size(curd.(algnuse).lick,1),1);
            mice.data{m}{d}.lickDiff(valid_trials) = mean(lickuse(valid_trials,latewin),2)-mean(lickuse(valid_trials,earlywin),2);
            assert(all(mean(lickuse(valid_trials,latewin),2)>=0) && all(mean(lickuse(valid_trials,earlywin),2)>=0))
            mice.data{m}{d}.lickPref(valid_trials) = mice.data{m}{d}.lickDiff(valid_trials)./(mean(lickuse(valid_trials,latewin),2)+mean(lickuse(valid_trials,earlywin),2));
            mice.data{m}{d}.lateLickFrac = nan(size(curd.(algnuse).lick,1),1);
            mice.data{m}{d}.lateLickFrac(valid_trials) = sum(lickuse_prob(valid_trials,latewin),2);

            end
        end
       
        if isfield(curd.midAlgn,"dlcnew")
            % pause
            windlc = winfn([0.3 1],curd.tmpxDLCnew);
            truse = find(curd.goodmvdir & curd.goodDLC & curd.mvlen>6);
            tmpdlc = nan(size(curd.midAlgn.dlcnew));
            for c= 1:8
                tmpdlc(truse,:,c) = squeeze(curd.midAlgn.dlcnew(truse,:,c));
            end
            if firstbody
                tmpxDLCnewtmplt = curd.tmpxDLCnew; firstbody=false;
            end
            trrew = intersect(truse,find(curd.rewarded));
            rewbodycur = tmpdlc(trrew,:,:);
            allbodycur = tmpdlc;
            rewbodycuri = zeros(size(rewbodycur,1),length(tmpxDLCnewtmplt),size(rewbodycur,3));
            for c = 1:8
                tmp = squeeze(rewbodycur(:,:,c));
                rewbodycur(:,:,c) = medfilt1(tmp,10,[],2);
                tmp = squeeze(allbodycur(:,:,c));
                tmptest = tmp(truse,:,:);
                assert(all(~isnan(tmptest(:))))
                allbodycur(truse,:,c) = zscore(medfilt1(tmp(truse,:,:),10,[],2),0,'all');
                tmptest = allbodycur(truse,:,:);
                assert(all(~isnan(tmptest(:))))
                % allbodycur(:,:,c) = medfilt1(tmp,10,[],2);
                rewbodycuri(:,:,c) = interp1(curd.tmpxDLCnew,squeeze(rewbodycur(:,:,c))',tmpxDLCnewtmplt,'linear','extrap')';
            end
            rewbodycuri_spd = zeros(size(rewbodycuri,1),size(rewbodycuri,2),size(rewbodycuri,3)/2);
            cuse = 1:2:7;
            for c = 1:length(cuse)
                cc = cuse(c);
                dtcam = mean(diff(tmpxDLCnewtmplt));
                tmp = sqrt(sum(diff(rewbodycuri(:,:,cc+(0:1)),[],2).^2,3))/dtcam;
                tmp = butter_filter(2,0.2,tmp')';
                rewbodycuri_spd(:,2:end,c) = zscore(tmp,0,'all');
            end
            [win22dlc_tmplt,win22dlc_tmplt_s] = winfn([-4 4],tmpxDLCnewtmplt);
            rewbodycuri_spd = rewbodycuri_spd(:,win22dlc_tmplt,:);
            for c = 1:8
                rewbodycuri(:,:,c) = zscore(rewbodycuri(:,:,c),0,'all');
            end
            rewbodycuri = rewbodycuri(:,win22dlc_tmplt,:);
        end
        if ismember(curep,epochs)
        curepn =find(curep==epochs);
        if curepn==5, curepn=4; end
        curep = curepn;
        if mice.days{m}.trueDay(d)==min(mice.days{m}.trueDay), curep = 1; end
        truse = find(curd.goodmvdir);
        trsucc = find(curd.goodmvdir & curd.mvlen>=7);
        behaviorQuantsEpoch.ntrials{curep} = cat(1,behaviorQuantsEpoch.ntrials{curep},length(find(curd.mvlen>=7)));
        behaviorQuantsEpoch.attempts_per_min{curep} = cat(1,behaviorQuantsEpoch.attempts_per_min{curep},length(curd.mvlen)/(curd.ntb*curd.dtb/60));
        behaviorQuantsEpoch.successes_per_min{curep} = cat(1,behaviorQuantsEpoch.successes_per_min{curep},length(find(curd.mvlen>=7))/(curd.ntb*curd.dtb/60));
        behaviorQuantsEpoch.success_rate{curep} = cat(1,behaviorQuantsEpoch.success_rate{curep},length(find(curd.mvlen>=7))/length(curd.mvlen));
        mice.data{m}{d}.success_rate = length(find(curd.mvlen>=7))/length(curd.mvlen);
        mice.data{m}{d}.successes_per_min = length(find(curd.mvlen>=7))/(curd.ntb*curd.dtb/60);
        behaviorQuantsEpoch.mvLen{curep} = cat(1,behaviorQuantsEpoch.mvLen{curep},curd.mvlen);
        behaviorQuantsEpoch.ITI{curep}=cat(1,behaviorQuantsEpoch.ITI{curep},mean(tmpITI));
        behaviorQuantsEpoch.ITI_t{curep}=cat(1,behaviorQuantsEpoch.ITI_t{curep},tmpITI);
        behaviorQuantsEpoch.ITI_omit{curep}=cat(1,behaviorQuantsEpoch.ITI_omit{curep},tmpITI_omit);
        mice.data{m}{d}.mvDur = (curd.trueend-curd.truestart)*curd.dtb;
        behaviorQuantsEpoch.mvDur{curep}=cat(1,behaviorQuantsEpoch.mvDur{curep},median(mice.data{m}{d}.mvDur(truse)));
        behaviorQuantsEpoch.mvDur_t{curep}=cat(1,behaviorQuantsEpoch.mvDur_t{curep},mice.data{m}{d}.mvDur(truse));

        behaviorQuantsEpoch.STreach{curep} = cat(1,behaviorQuantsEpoch.STreach{curep},curd.startAlgn.pos(:,tmpwin2,2));
        if licksensor_valid{m}(d)~="n"
            behaviorQuantsEpoch.anticLick{curep} = cat(1,behaviorQuantsEpoch.anticLick{curep},mice.data{m}{d}.anticLick(valid_trials_r));
            if algnuse=="rewAlgn",tmpwin2 = winfn([-2 2],curd.tmpxx);
            else, tmpwin2 = winfn([-1 4],curd.tmpxx);
            end
            % if isempty(valid_trials_o), fprintf("no omitted trials: %s %s %s\n",mice.name(m),mice.days{m}.date{d},mice.days{m}.epoch(d)); end
            lick_use = cur_filtlick(valid_trials_r,tmpwin2);
            behaviorQuantsEpoch.STlick{curep} = cat(1,behaviorQuantsEpoch.STlick{curep},lick_use);
            lick_use = cur_filtlick(valid_trials_o,tmpwin2);
            behaviorQuantsEpoch.STlick_o{curep} = cat(1,behaviorQuantsEpoch.STlick_o{curep},lick_use);
            
            lick_use = cur_filtrate(valid_trials,tmpwin2);
            behaviorQuantsEpoch.STrate_all{curep} = cat(1,behaviorQuantsEpoch.STrate_all{curep},lick_use);
            lick_use = cur_filtrate(valid_trials_r,tmpwin2);
            behaviorQuantsEpoch.STrate{curep} = cat(1,behaviorQuantsEpoch.STrate{curep},lick_use);
            lick_use = cur_filtrate(valid_trials_o,tmpwin2);
            behaviorQuantsEpoch.STrate_o{curep} = cat(1,behaviorQuantsEpoch.STrate_o{curep},lick_use);

            lick_use = cur_rawlick(valid_trials_r,tmpwin2);
            behaviorQuantsEpoch.STlick_raw{curep} = cat(1,behaviorQuantsEpoch.STlick_raw{curep},lick_use);
            lick_use = cur_rawlick(valid_trials_o,tmpwin2);
            behaviorQuantsEpoch.STlick_raw_o{curep} = cat(1,behaviorQuantsEpoch.STlick_raw_o{curep},lick_use);
            licksessioncount(curep) = licksessioncount(curep) +1;
        
            lick_use = cur_probrate(valid_trials,tmpwin2);
            assert(size(cur_filtrate,1)==size(cur_probrate,1));
            behaviorQuantsEpoch.STrate_prob{curep} = cat(1,behaviorQuantsEpoch.STrate_prob{curep},lick_use*100);

        behaviorQuantsEpoch.lick_mice{curep} = unique([behaviorQuantsEpoch.lick_mice{curep},mice.name(m)]);
        end
        if isfield(curd.midAlgn,"dlcnew")
            behaviorQuantsEpoch.STbody{curep}=cat(1,behaviorQuantsEpoch.STbody{curep},rewbodycuri);
            behaviorQuantsEpoch.STbody_speed{curep}=cat(1,behaviorQuantsEpoch.STbody_speed{curep},rewbodycuri_spd);
        end
        end
    end
end
%% Fig S2A-E
cqs = ["ITI","mvDur_t","attempts_per_min","successes_per_min","anticLick"];
titls = ["S2C","S2A","S2D","S2B","S2E"];
for cq = 1:length(cqs)
    curquant = cqs(cq);
fprintf("%s: ",curquant)
junk = cat(1,behaviorQuantsEpoch.(curquant){:});
cur_n = length(junk);
if cur_n<2000
barwitherr_stde(behaviorQuantsEpoch.(curquant),dots=true,overlap=1.2,title=titls(cq));
xticklabels(eplbls); 
ylabel(curquant); end
if cur_n>2000
    histogram_cell(behaviorQuantsEpoch.(curquant),[],cdf=true,title=titls(cq))
    legend(eplbls); 
    xlabel(curquant)
end
ranksum_disp(behaviorQuantsEpoch.(curquant){[1 4]})
alleps = arrayfun(@(x)x*ones(length(behaviorQuantsEpoch.(curquant){x}),1),1:4,unif=false);
alleps = cat(1,alleps{:});
anova_disp(cat(1,behaviorQuantsEpoch.(curquant){:}),y=alleps);
end

figure(pos=[ 584   512   257   206]); h1=errorbar_shadeSEM(curd.tmpxx(tmpwin2),behaviorQuantsEpoch.STrate{1}',[0.47 .67 .19]);
hold on; h2=errorbar_shadeSEM(curd.tmpxx(tmpwin2),behaviorQuantsEpoch.STrate{4}',[0.49 0.18 0.57]);
timax_str(algnuse);
ylabel("lick rate")
ylim tight;
yl = ylim; ylim([0 yl(2)]);
plot(-[1 1]*meanmvtime_all,ylim,'k');
plot([0 0],ylim,'k');
legend([h1,h2],["day 1","expert"])
s=sprintf("n: %d %d",size(behaviorQuantsEpoch.STrate{1},1),size(behaviorQuantsEpoch.STrate{4},1));
s="2A"+newline+s;
title(s);

figure(pos=[ 584   512   257   206]); 
h1=errorbar_shadeSEM(curd.tmpxx(tmpwin2),behaviorQuantsEpoch.STrate{4}',[0 0 1]);
hold on; h2=errorbar_shadeSEM(curd.tmpxx(tmpwin2),behaviorQuantsEpoch.STrate_o{4}',[1 0 0]);
timax_str(algnuse);
ylabel("lick rate (hz)")
ylim tight;
yl = ylim; ylim([0 yl(2)]);
plot(-[1 1]*meanmvtime_all,ylim,'k');
plot([0 0],ylim,'k');
legend([h1,h2],["expert rewarded","expert omission"])
s=sprintf("n: %d %d",size(behaviorQuantsEpoch.STrate{4},1),size(behaviorQuantsEpoch.STrate_o{4},1));
s="2D"+newline+s;
title(s);

figure(pos=[309    48   212   287]); 
ccplt = [0.4667    0.6745    0.1882; 0.1020    0.1020    0.5020;...
         0.2510    0.2510    0.7490; 0.4941    0.1843    0.5569];
h0=errorbar_shadeSEM(curd.tmpxx(tmpwin2),behaviorQuantsEpoch.STrate_prob{1}',ccplt(1,:));
hold on; 
h1=errorbar_shadeSEM(curd.tmpxx(tmpwin2),behaviorQuantsEpoch.STrate_prob{2}',ccplt(2,:));
h2=errorbar_shadeSEM(curd.tmpxx(tmpwin2),behaviorQuantsEpoch.STrate_prob{3}',ccplt(3,:));
h3=errorbar_shadeSEM(curd.tmpxx(tmpwin2),behaviorQuantsEpoch.STrate_prob{4}',ccplt(4,:));
uistack([h0,h1,h2,h3],"top");
timax_str(algnuse);
ylabel("relative licking probability");
xlim([-meanmvtime 0.0])
ylim tight;
yl = ylim; ylim([0 yl(2)]);
s=sprintf("n: %d %d %d %d",arrayfun(@(x)size(behaviorQuantsEpoch.STrate_prob{x},1),1:4));
s="2B"+newline+s;
title(s);


tmpwin_reach = winfn([-0.2 0.6],curd.tmpxx(tmpwin2));
figure(); errorbar_shadeSEM(curd.tmpxx(tmpwin2(tmpwin_reach)),behaviorQuantsEpoch.STreach{1}(:,tmpwin_reach)',[0.47 .67 .19]);
hold on; errorbar_shadeSEM(curd.tmpxx(tmpwin2(tmpwin_reach)),behaviorQuantsEpoch.STreach{4}(:,tmpwin_reach)',[0.49 0.18 0.57]);
timax_str("startAlgn");
ylabel("reach distance (mm)")
s=sprintf("day1 %d day7 %d",size(behaviorQuantsEpoch.STreach{1},1),size(behaviorQuantsEpoch.STreach{4},1));
s="S2A"+newline+s;
title(s);
axis tight;

fprintf("total sessions: %d\n",length(cat(1,behaviorQuantsEpoch.ITI{:})))
fprintf("lick data sessions: %d\n",licksessioncount)
fprintf("lick data sessions total: %d\n",sum(licksessioncount))
behaviorQuantsEpoch.lick_mice
cellfun(@(x)length(x),behaviorQuantsEpoch.lick_mice)
%% Fig 2C, S2E, 2D-inset
for i = 1
latemearly = cell(1,nep);
late = cell(1,nep);
latewin = winfn([-0.2 0],curd.tmpxx(tmpwin2));
earlywin = winfn([-0.8 -0.6],curd.tmpxx(tmpwin2));
lickuse = behaviorQuantsEpoch.STrate_all;
lickuse_prob = behaviorQuantsEpoch.STrate_prob;

for g = 1:nep
    latemearly{g} = mean(lickuse{g}(:,latewin),2)-mean(lickuse{g}(:,earlywin),2);
    % latemearly{g} = latemearly{g}./(mean(lickuse{g}(:,latewin),2)+mean(lickuse{g}(:,earlywin),2));
    % latemearly{g}(isnan(latemearly{g}))=0;
    late{g} = sum(lickuse_prob{g}(:,latewin),2);
end
histogram_cell(latemearly,cdf=true);
hold on; plot([0 0],ylim,'k')
legend(eplbls)
titstr=mean_stde_disp(latemearly([1 4]));
titstr=titstr+newline+"2C"; title(titstr);
signrank_disp(latemearly{4})
end
lickuses = {behaviorQuantsEpoch.STrate{4},behaviorQuantsEpoch.STrate_o{4}};
laterewwin = winfn([0.75 1],curd.tmpxx(tmpwin2)); laterewlick = {};
postrewofftimes = {[],[]};
trwin = winfn([-0.25 2],curd.tmpxx(tmpwin2));
mxwin = winfn([-0.25 0.5],curd.tmpxx(tmpwin2));
for i = 1:2
laterewlick{i} = mean(lickuses{i}(:,laterewwin),2);
for k = 1:size(lickuses{i},1)
    curtr = lickuses{i}(k,:); curtr = filtfilt(ones(50,1)/50,1,curtr')';
    [mx,mxix] = max(curtr(mxwin));
    firstfall = find(curtr(trwin(mxix:end))<mx*0.5,1,'first')+mxix-1;
    lasthigh = find(curtr(trwin(mxix:end))>mx*0.7,1,'last')+mxix-1;
    if isempty(firstfall), lasthigh=0;
    else, lasthigh = curd.tmpxx(tmpwin2(trwin(lasthigh)));
    end
    postrewofftimes{i}(end+1) = lasthigh;
end
end
histogram_cell(laterewlick,[],cdf=true,title="S2E");
histogram_cell(postrewofftimes,[],cdf=true,title="2D");
%% FIG. 3M and S2T
close all
miceuse = 1:nmice;
nmiceuse = length(miceuse);
cfsiguse = "sp_CF";
sigs = ["sigFilt_GrC",cfsiguse];
epoch_use = "expert";

ntg = 2;
grcsigs_session = cell(ntg,0);
cfsigs_session = cell(ntg,0);
% average across specific cells over specific windows or differences between windows per trial
quantheaders = ["name","algnuse","sigsUse","avgWins","cellSubset","quantFn","results","trTypes","session_ids"];
s_quantTyps = cell2table(cell(0,length(quantheaders)),VariableNames=quantheaders);
s_quantTyps = cat(1,s_quantTyps,{"CFpostRespTime","rewAlgn",{["sp_CF"]},{[0.02 1]},...
    {cfsigs_all_ix{2}{1}(proportions_ix_3gm{2},:)},"spiketime",...
    cell(1,1),cell(1,1),cell(1,1)}); %"valid_trials_o"
s_quantTyps = cat(1,s_quantTyps,{"GrCanticOffTime","rewAlgn",{["sigFilt_GrC"]},...
    {[-0.35 0; 0 1.5]},{grcsigs_all_ix{2}{1}(proportions_ix_3gm{1},:)},"offtime",cell(1,1),cell(1,1),cell(1,1)});
s_quantTyps.Properties.RowNames=s_quantTyps.name;

clear session_ids;
sessioncount=1;
grccount_total = 0;
cfcount_total = 0;
for mm = miceuse
    daysuse_cur = find(ismember(mice.days{mm}.epoch,epoch_use));
    for dd = 1:length(daysuse_cur)
        dcur = daysuse_cur(dd);
        curd = mice.data{mm}{dcur};
        zerptCb = minix(abs(curd.tmpxCb));
        dtimCb = curd.dtimCb;
        tmpxCb = curd.tmpxCb;
        clear trgrps;

        for q = 1:size(s_quantTyps,1)
            % rewarded / unrewarded
            trgrps{1} = find(curd.rewarded & curd.mvlen>7 & curd.goodmvdir);
            trgrps{2} = find(~curd.rewarded & curd.mvlen>7 & curd.goodmvdir);
            trgrps{1} = [trgrps{1};trgrps{2}];
            cur_rew_ix = find(curd.rewarded(trgrps{1}));
            cur_omit_ix = find(~curd.rewarded(trgrps{1}));
            if length(trgrps{1})>2
            for cg = 1:length(s_quantTyps.sigsUse{q})
                curtr_sigs = curd.(s_quantTyps.algnuse{q}).(s_quantTyps.sigsUse{q}(cg));
                curtr_sigs = curtr_sigs(trgrps{1},:,:);
                nccur = size(curtr_sigs,2);
                if cfsiguse=="sp_CF"
                    tmpspvec = reshape(permute(curtr_sigs,[3 1 2]),[],nccur);
                    spkern = round(0.15 / curd.dtimCb);
                    spkern_s = spkern*curd.dtimCb;
                    tmpspvecf = filtfilt(ones(spkern,1)/spkern/spkern_s,1,double(tmpspvec));
                    tmpspvecf = zscore(tmpspvecf,[],1);
                    tmpspvecf = reshape(tmpspvecf,size(curtr_sigs,3),size(curtr_sigs,1),nccur);
                    curtr_sigs = permute(tmpspvecf,[2 3 1]);
                end
                if ~all(all((ismissing(s_quantTyps.cellSubset{q}))))
                    cuse = find(s_quantTyps.cellSubset{q}(:,1)==mm & s_quantTyps.cellSubset{q}(:,2)==dcur);
                    if ~isempty(cuse)
                        cuse = s_quantTyps.cellSubset{q}(cuse,3);
                        curtr_sigs = squeeze(mean(curtr_sigs(:,cuse,:),2));
                    else
                        curtr_sigs = [];
                    end
                elseif ismissing(s_quantTyps.regweights{q}), curtr_sigs = squeeze(mean(curtr_sigs,2));
                else
                    if contains(s_quantTyps.sigsUse{q}(cg),"CF"), cgix = find(contains(regTyps.sigsUse{s_quantTyps.regweights{q}},"CF")); else, cgix=1; end
                    curweights = regTyps.regCoefs{s_quantTyps.regweights{q}}{mm}{dcur}{cgix};
                    curtr_sigs = squeeze(sum(bsxfun(@times,permute(curweights(2:end),[3 1 2]),curtr_sigs),2));
                end
                if ~isempty(curtr_sigs)
                win1 = winfn(s_quantTyps.avgWins{q}(1,:),tmpxCb);
                tr_sig = {};
                switch s_quantTyps.quantFn{q}
                case "offtime"
                    cur_kern = ceil(0.2/dtimCb);
                    win2 = winfn(s_quantTyps.avgWins{q}(2,:),tmpxCb);
                    cur_smoothed_sigs = filtfilt(ones(cur_kern,1)/cur_kern,1,curtr_sigs')';
                    tmptrid = zeros(size(curtr_sigs,1),1); tmptrid(cur_rew_ix)=1;
                    s_quantTyps.trTypes{q}{cg}{sessioncount} = tmptrid;
                    if 1
                    cur_max = max(cur_smoothed_sigs(:,win1),[],2);
                    cur_smoothed_sigs=cur_smoothed_sigs(:,win2);
                    junk_antic = curtr_sigs(cur_rew_ix,winfn([-0.5 -0.05],tmpxCb));
                    threshval = max(0.2,cur_max*0.5);
                    tr_sig{1} = zeros(size(curtr_sigs,1),1);
                    cur_kern = ceil(0.4/dtimCb);
                    cur_raw_sigs = filtfilt(ones(cur_kern,1)/cur_kern,1,curtr_sigs')';
                    cur_raw_sigs = cur_raw_sigs(:,win2);
                    for ii = 1:size(curtr_sigs,1)
                        cur_ix_off = find(cur_raw_sigs(ii,:)<threshval(ii),1);
                        if isempty(cur_ix_off) | threshval(ii)<prctile(junk_antic(:),5)
                            tr_sig{1}(ii) = tmpxCb(win2(end));
                        else
                            tr_sig{1}(ii) = tmpxCb(win2(cur_ix_off));
                        end
                    end                       
                    end
                    fprintf("grc fail frac: %.3g, rew off times: %.3g omit off times: %.3g\n",sum(isnan(tr_sig{1}))/length(tr_sig{1}),nanmean(tr_sig{1}(cur_rew_ix)),nanmean(tr_sig{1}(cur_omit_ix)))
                case "spiketime"
                    tmpsig = curtr_sigs(:,win1);
                    tr_sig{1} = zeros(size(tmpsig,1),1);
                    failcount = 0;
                    tmpjunk = curtr_sigs(cur_rew_ix,winfn([0 0.2],tmpxCb));
                    junk = maxix(mean(tmpjunk));
                    curthresh = max(0,prctile(tmpjunk(:,junk),20));
                    for ii = 1:size(tmpsig,1)
                        warning off;
                        [pkmags,tmppks] = findpeaks(tmpsig(ii,:),MinPeakHeight=curthresh);
                        warning on;
                        if isempty(tmppks)
                            tmpfind = maxix(tmpsig(ii,:));
                            failcount = failcount+1;
                            tr_sig{1}(ii) = nan;
                        else
                            tmpfind = tmppks(1);
                            tr_sig{1}(ii) = tmpxCb(win1(tmpfind));
                        end
                    end
                    fprintf("CF threshold %.3g find fail faction: %.3g, ",curthresh,failcount/size(tmpsig,1))
                    fprintf("rew peak times: %.3g omit peak times: %.3g\n",nanmean(tr_sig{1}(cur_rew_ix)),nanmean(tr_sig{1}(cur_omit_ix)))
                    [~,tmpsortix] = sort(tr_sig{1});
                    if s_quantTyps.name(q)=="CFpostRespTime"
                        tmp_plt_cf_im_srt = curtr_sigs(:,winfn([-2 2],tmpxCb));
                        tmp_plt_ix_tr_srt = tmpsortix;
                        tmp_plt_ix_tr_vals = tr_sig{1};
                    end
                end
                clear cur_reslt;
                cur_reslt = tr_sig{1};
                s_quantTyps.results{q}{cg}{sessioncount} = cur_reslt;

                s_quantTyps.session_ids{q} = cat(1,s_quantTyps.session_ids{q},[mice.name{mm} dcur mice.days{mm}{dcur,["date" "epoch"]}]);
            end
            end
        end
        end
        session_ids(sessioncount,:) = [mice.name{mm} dcur mice.days{mm}{dcur,["date" "epoch"]}];
        sessioncount = sessioncount+1;
    end
end
% session quantification correlations
q1 = s_quantTyps.results{"GrCanticOffTime"}{1};
sq1 = s_quantTyps.session_ids{"GrCanticOffTime"};
q2 = s_quantTyps.results{"CFpostRespTime"}{1};
sq2 = s_quantTyps.session_ids{"CFpostRespTime"};
trTyps = s_quantTyps.trTypes{"GrCanticOffTime"}{1};
ixuse = intersect(find(cellfun(@(x)~isempty(x),q1)),find(cellfun(@(x)~isempty(x),q2)));
sess_str = sprintf("nsessions q1: %d, q2: %d, both: %d\n",sum(cellfun(@(x)~isempty(x),q1)),sum(cellfun(@(x)~isempty(x),q2)),length(ixuse));
sess_str = sess_str+sprintf("nmice q1: %d, q2: %d, both: %d\n",...
    length(unique(sq1(:,1))),length(unique(sq2(:,1))),length(intersect(unique(sq2(:,1)),unique(sq1(:,1)))));
disp(sess_str);
ns = length(ixuse);
allcorrs_sessions = zeros(ns,1);
allcorrs_rew = zeros(ns,1);
allnvals_sessions = zeros(ns,1);
allomit_means = zeros(ns,2);
allrew_means = zeros(ns,2);
alltr_use = cell(ns,1);
for s = 1:ns
    tr_use = find(~isnan(q1{ixuse(s)}) & ~isnan(q2{ixuse(s)}));
    alltr_use{s} = tr_use;
    allnvals_sessions(s,1) = length(tr_use);
    cur_omit = find(trTyps{ixuse(s)}(tr_use)==0);
    cur_rew = find(trTyps{ixuse(s)}(tr_use)==1);
    allnvals_sessions(s,2:3) = [length(cur_rew) length(cur_omit)];
    allomit_means(s,:) = mean([q1{ixuse(s)}(tr_use(cur_omit)),q2{ixuse(s)}(tr_use(cur_omit))],1);
    allrew_means(s,:) = mean([q1{ixuse(s)}(tr_use(cur_rew)),q2{ixuse(s)}(tr_use(cur_rew))],1);
    junk=corr([q1{ixuse(s)}(tr_use(cur_rew)),q2{ixuse(s)}(tr_use(cur_rew))]);
    allcorrs_rew(s,:) = junk(1,2);
    junk=corr([q1{ixuse(s)}(tr_use([cur_rew; cur_omit])),q2{ixuse(s)}(tr_use([cur_rew; cur_omit]))]);
    allcorrs_sessions(s,:) = junk(1,2);
end
fprintf("# trials per session: "),mean_stde_disp(allnvals_sessions);
fprintf("num mice: %d\n",length(unique(session_ids(ixuse,1))));

figure;
plot([allrew_means(:,1) allomit_means(:,1)]',[allrew_means(:,2) allomit_means(:,2)]','color',[0.8 0.8 0.8])
hold on;
plot([allrew_means(:,1)]',[allrew_means(:,2)]','b.',MarkerSi=15)
plot([allomit_means(:,1)]',[allomit_means(:,2)]','r.',MarkerSi=15)
axis square; box off;
titstr = sprintf("%d sessions\n",size(allrew_means,1))+sess_str;
titstr = [titstr; "Fig3M"];
title(titstr);
xlabel("Anticipation GrC Off-Time (s)")
ylabel("First CF population spike time (s)")

s1=signrank_disp(allomit_means(:,1),allrew_means(:,1));
s2=signrank_disp(allomit_means(:,2),allrew_means(:,2));

barwitherr_stde({allrew_means(:,1),allomit_means(:,1),allrew_means(:,2),allomit_means(:,2)},...
    color=[0 0 1; 1 0 0; 0 0 1; 1 0 0]);
titstr = sprintf("GrCs: %sCFs: %s",s1,s2)+sess_str;
titstr = [titstr; "Fig3M"];
title(titstr);
ylabel("GrC off-time/CF spike time (s)")
xticks([1.5 3.5]); xticklabels(["GrCs","CFs"])

barwitherr_stde({allcorrs_sessions(:,1),allcorrs_rew(:,1)},dots=true);
titstr = signrank_disp(allcorrs_sessions)+" "+signrank_disp(allcorrs_rew)+sess_str+"trial-by-trial correlations";
titstr = [titstr; "Fig3M"];
title(titstr);
xticklabels(["All trials","Rewarded trials"])
%% Fig S3 lick triggered average
jitters = [0 1];
[lick_trig_all,grcs_lick_trig_all,grcs_all,grcs_all_omit,...
    grcs_all_ix] = deal({});
firstrun = true;
for k_run = 1:2
do_jitter = jitters(k_run);
[lick_trig_all{k_run}, grcs_lick_trig_all{k_run}, grcs_all{k_run},...
    grcs_all_omit{k_run},grcs_all_ix{k_run}] = deal([]);
s=1;
for m = 1:nmice
    for d = 1:mice.ndays(m)
        if licksensor_valid{m}(d) ~= "n" && ...
                mice.days{m}.epoch(d)=="expert" && ...
                (isnan(mice.days{m}.chronDay(d)) || mice.days{m}.chronDay(d)==max(mice.days{m}.chronDay(d)))
        curd = mice.data{m}{d};
        valid_trials = mice.data{m}{d}.valid_lick_trials;
        assert(all(curd.rewtimes(valid_trials)-curd.midpt(valid_trials)<2.5/curd.dtb))
        tr_rew = find(curd.rewarded & curd.mvlen>7 & curd.goodmvdir); %rewarded trials
        tr_omit = find(~curd.rewarded & curd.mvlen>7 & curd.goodmvdir); %omitted reward trials
        tr_rew = intersect(tr_rew,valid_trials);
        tr_omit = intersect(tr_omit,valid_trials);
        trtyp = cat(1,zeros(size(tr_rew)),ones(size(tr_omit)));
        tr_use = cat(1,tr_rew,tr_omit);
        ixuse = winfn(curd.tmpxx,[-1 0]);
        lickuse = curd.rewAlgn.lick(tr_use,ixuse);
        trigwinsz = 0.4;
        trigwinlen = round(trigwinsz/curd.dtb);
        trigwinlen_cb = round(trigwinsz/curd.dtimCb);
        [ix_all,ix_all_s] = winfn([-2 2],curd.tmpxCb);
        trigax_cur = (-trigwinlen_cb:trigwinlen_cb)*curd.dtimCb;
        trigax_b = (-trigwinlen:trigwinlen)*curd.dtb;
        if firstrun,
            trigax_tmplt = trigax_cur;
            ix_all_tmplt = ix_all_s;
            firstrun=0;
        end
        lick_trig = [];
        grcs_lick_trig = [];
        for k = 1:length(tr_use)
            ixuse_cur = ixuse;
            lickonsets = find(diff(lickuse(k,:))>0)';
            lickonsets(lickonsets<trigwinlen | lickonsets>(length(ixuse_cur)-trigwinlen)) = [];
            if ~isempty(lickonsets)
                % pause
            lickonsets_s = curd.tmpxx(ixuse_cur(lickonsets))';
            if do_jitter, lickonsets_s = lickonsets_s+(rand(size(lickonsets_s))-0.5)*0.3; end
            lickonsets_cb = round(lickonsets_s/curd.dtimCb)+minix(abs(curd.tmpxCb));
            lickonsets_b = round(lickonsets_s/curd.dtb)+minix(abs(curd.tmpxx));
            lickwins_b = lickonsets_b + (-trigwinlen :trigwinlen);
            lickwins = lickonsets_cb + (-trigwinlen_cb :trigwinlen_cb);
            grcveclen = size(curd.rewAlgn.sigFilt_GrC,3);
            lickveclen = size(curd.rewAlgn.lick,2);
            if (nnz(lickwins(lickwins>grcveclen))>0), disp(nnz(lickwins(lickwins>grcveclen))>0); end
            if (nnz(lickwins_b(lickwins_b>lickveclen))>0), disp(nnz(lickwins_b(lickwins_b>lickveclen))); end
            lickwins(lickwins<1) = 1; lickwins(lickwins>grcveclen) = grcveclen;
            lickwins_b(lickwins_b<1) = 1; lickwins_b(lickwins_b>lickveclen) = lickveclen;
            tmptrsig = squeeze(curd.rewAlgn.sigFilt_GrC(tr_use(k),:,:));
            for j = 1:size(lickwins,1)
                trigsig_interp = interp1(trigax_cur,tmptrsig(:,lickwins(j,:))',trigax_tmplt,"linear","extrap")';
                grcs_lick_trig = cat(3,grcs_lick_trig,trigsig_interp);
                lick_trig = cat(1,lick_trig,curd.rewAlgn.lick(tr_use(k),lickwins_b(j,:)));
            end
            end
        end
        if 1
        grcs_all_cur = squeeze(mean(curd.rewAlgn.sigFilt_GrC(tr_rew,:,ix_all),1));
        ix_b = winfn([-2 2],curd.tmpxx);
        lick_trig_all{k_run} = cat(1,lick_trig_all{k_run},mean(lick_trig,1));
        grcs_lick_trig_all{k_run} = cat(1,grcs_lick_trig_all{k_run},mean(grcs_lick_trig,3));
        grcs_all{k_run} = cat(1,grcs_all{k_run},interp1(ix_all_s,grcs_all_cur',ix_all_tmplt,"linear","extrap")');
        grcs_all_ix{k_run} = cat(1,grcs_all_ix{k_run},[repmat([m d s],curd.nIC_GrC,1) (1:curd.nIC_GrC)']);
        s=s+1;
        end
        % curstr=sprintf('m:%d d:%d',m,d);
        % fprintf(curstr);
        end
    end
end
end

[~,ixant] = ismember(ixantic_ss{1},grcs_all_ix{1}(:,[1 2 4]),'rows');
ixant = ixant(ixant~=0);
ixant_true = ixant;
ixcorr_b = winfn(trigax_b,[-0.2 0.2]);
ixcorr_cb = winfn(trigax_tmplt,[-0.2 0.2]);
[tmpr,tmpp] = corr(grcs_lick_trig_all{1}(:,ixcorr_cb)',...
    interp1(trigax_b(ixcorr_b),lick_trig_all{1}(grcs_all_ix{1}(:,3),ixcorr_b)',...
    trigax_tmplt(ixcorr_cb),'linear','extrap'));
tmpr = diag(tmpr); tmpp = diag(tmpp);
ixant_neg = find(tmpr<0 & tmpp<0.05);
ixant_pos = find(tmpr>0 & tmpp<0.05);
ix_posneg = [ixant_neg; ixant_pos];
length(ixant_true)*length(ix_posneg)/size(grcs_all_ix{1},1).^2;
length(intersect(ixant_true,ix_posneg))/size(grcs_all_ix{1},1);

ixant_srt = mean(grcs_lick_trig_all{1}(ixant,round(end/2)+[-1 1]),2);
[~,ixant_srt] = sort(ixant_srt);
rng(0);
rpp = rand(1,3);
rpn = rand(1,3);
figure; h1=errorbar_shadeSEM(ix_all_tmplt,grcs_all{1}(ixant_pos,:)',rpp); axis tight;
hold on; h2=errorbar_shadeSEM(ix_all_tmplt,grcs_all{1}(ixant_true,:)',[0 0 1]); axis tight;
hold on; h3=errorbar_shadeSEM(ix_all_tmplt,grcs_all{1}(ixant_neg,:)',rpn); axis tight;
hold on; plot([0 0],ylim,'k'); plot(-meanmvtime_all*[1 1],ylim,'k')
legend([h1 h2 h3],[sprintf("%d %0.2g",length(ixant_pos),length(ixant_pos)./size(grcs_all_ix{1},1)),...
    sprintf("%d %.2g",length(ixant_neg),length(ixant_neg)./size(grcs_all_ix{1},1)),...
    sprintf("%d %0.2g",length(ixant_true),length(ixant_true)./size(grcs_all_ix{1},1))])
figureformat_forsaving
figure(pos=[1504         322         171         431]); 
subaxis(3,1,1,'ML',0.2,'MB',0.2,'SV',0.07); 
hl2=errorbar_shadeSEM(trigax_b,lick_trig_all{2}',[0.5 0.5 0.5],linewidth=1,linespec='--');
hl1=errorbar_shadeSEM(trigax_b,lick_trig_all{1}',[0 0 0]);
uistack([hl1,hl2],'top');
axis tight;
ylabel("Lick probability");
hold on; subaxis(3,1,1,3); 

hl2=errorbar_shadeSEM(trigax_tmplt,grcs_lick_trig_all{1}(ixant_pos,:)',rpp);
hl1=errorbar_shadeSEM(trigax_tmplt,grcs_lick_trig_all{1}(ixant_neg,:)',rpn);
uistack([hl1,hl2],'top');
axis tight;
ylabel("GrC fluorescence (zsc)");
xlabel(["Time relative to lick onset","for delay licks (ms)"])
title(sprintf("%d mice %d sessions %d neg %d pos",length(unique(grcs_all_ix{1}([ixant_pos;ixant_neg],1))),...
    length(unique(grcs_all_ix{1}([ixant_pos;ixant_neg],3))),size(ixant_pos,1),size(ixant_neg,1)))

hold on; subaxis(3,1,1,2); 
hl2=errorbar_shadeSEM(trigax_tmplt,grcs_lick_trig_all{2}(ixant_true,:)',[0.5 0.5 0.5],linewidth=1,linespec='--');
hl1=errorbar_shadeSEM(trigax_tmplt,grcs_lick_trig_all{1}(ixant_true,:)',[0 0 1]);
uistack([hl1,hl2],'top');
axis tight;
ylabel("GrC fluorescence (zsc)");
title(sprintf("%d mice %d sessions %d cells",length(unique(grcs_all_ix{1}(ixant_true,1))),length(unique(grcs_all_ix{1}(ixant_true,3))),size(ixant_true,1)))

figureformat_forsaving
%% DLC Fig S3L,M
dlccoorduse = [1:2 5:8];
cc = rand(length(dlccoorduse),3)*0.8;
tmp = behaviorQuantsEpoch.STbody{"expert"};
figure(pos=[1063         150         425         710]);
[windlccur,windlccur_s] = winfn(win22dlc_tmplt_s,[-3 3]);
for k = 1:length(dlccoorduse)
    kk = dlccoorduse(k);
    errorbar_shadeSEM(windlccur_s,(squeeze(tmp(:,windlccur,kk)))'+(k-1),cc(k,:));
    hold on;
end
yticks(0:(length(dlccoorduse)-1))
yticklabels(dlclabls(dlccoorduse))
ylim tight;
hold on; plot([0 0],ylim,'k');
hold on; plot([1 1]*meanmvsttime_all,ylim,'k--');
title(sprintf("%d trials",size(tmp,1)));

tmp = behaviorQuantsEpoch.STbody_speed{"expert"};
figure(pos=[1063         150         425         710]);
dlccoorduse = [1 3 4];
for k = 1:length(dlccoorduse)
    kk = dlccoorduse(k);
    errorbar_shadeSEM(windlccur_s,(squeeze(tmp(:,windlccur,kk)))'+(k-1)/6,cc(k,:));
    hold on;
end
yticks(0:(1/6):(length(dlccoorduse)-1))
yticklabels(dlclablshalf(dlccoorduse))
ylim tight;
hold on; plot([0 0],ylim,'k');
hold on; plot([1 1]*meanmvsttime_all,ylim,'k--');
title(sprintf("%d trials",size(tmp,1)));