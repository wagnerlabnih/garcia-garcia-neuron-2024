%% Load in pieces and metadata
clear; clc;
wd = "C:\Users\wagnermj\Desktop\data to share\";
mice = load(fullfile(wd,"main_metadat.mat")); mice = mice.mice2;
mdall = {};
for i = 1:4
    fn = fullfile(wd,"main_part"+num2str(i)+".mat");
    load(fn);
    mdall = cat(2,mdall,eval("md"+num2str(i)));
    fprintf("loaded part %d\n",i);
end
for i = 1:height(mice)
    for j = 1:mice.ndays(i)
        mice.data{i}{j} = mdall{i}{j};
    end
end
clear mdall;
%% initialize vars
% make split alignment
% concatenate data
% compute single cell parameters
% split lick sensor data
nmice = size(mice,1);
chronicmice = mice.name(arrayfun(@(x)~isempty(mice.chronAlgnGrC{x}),1:nmice));
eplist= ["day1","novice","mid","expert"];
cc_ep = [119 172 48; 26 26 128; 64 64 191; 126 47 142]./255;
% matlab data file saving issue produces a handful of nan frames
for m = 1:nmice
for k = 1:mice.ndays(m)
for algnuse=["startAlgn","midAlgn","rewAlgn"]
for sig = ["sigFilt_GrC","sigFilt_CF","sp_CF"]
    if any(isnan(mice.data{m}{k}.(algnuse).(sig)(:)))
        fprintf("m%d d%d %s replacing %d nans with 0\n",m,k,sig,...
            nnz(isnan(mice.data{m}{k}.(algnuse).(sig)(:)))); 
        mice.data{m}{k}.(algnuse).(sig)(isnan(mice.data{m}{k}.(algnuse).(sig)(:)))=0;
    end
end
end
end
end
for m = 1:nmice
    mice.nGrCtot(m) = 0;
    mice.nCFtot(m) = 0;
    tmpdaysuse = 1:mice.ndays(m);
    for kk = 1:length(tmpdaysuse)
        k = tmpdaysuse(kk);
        mice.nGrCtot(m) = mice.nGrCtot(m) + mice.data{m}{k}.nIC_GrC;
        mice.nCFtot(m) = mice.nCFtot(m) + mice.data{m}{k}.nIC_CF;
        mice.days{m}(kk,"nGrC") = {mice.data{m}{k}.nIC_GrC};
        mice.days{m}(kk,"nCF") = {mice.data{m}{k}.nIC_CF};
    end
end
nCFtot = sum(cat(1,mice.nCFtot));
nGrCtot = sum(cat(1,mice.nGrCtot));

for m = 1:nmice
    for k = 1:mice.ndays(m)
    mvstart = minix(abs(mice.data{m}{k}.tmpxx));
    mice.data{m}{k}.mvlen = sqrt(sum(mice.data{m}{k}.rewAlgn.pos(:,mvstart,:).^2,3));
    lengrc= size(mice.data{m}{k}.sigFilt_GrC,2);
    lencf = size(mice.data{m}{k}.sigFilt_CF,2);
    if lengrc ~= mice.data{m}{k}.ntimCb || lencf ~= mice.data{m}{k}.ntimCb
        fprintf("%s %s grc sig length: %d, cf sig length: %d\n",mice.name(m),mice.days{m}.date(k),lengrc,lencf)
        s = input(sprintf("'y' to crop longer one to %d",mice.data{m}{k}.ntimCb),'s');
        if s=="y", 
            mice.data{m}{k}.sigFilt_GrC = mice.data{m}{k}.sigFilt_GrC(:,1:mice.data{m}{k}.ntimCb);
            mice.data{m}{k}.sigFilt_CF = mice.data{m}{k}.sigFilt_CF(:,1:mice.data{m}{k}.ntimCb);
        end
    end
    end
end
clear tmpdaysuse mvstart lencf lengrc kk k j i

% define windows, labels
ccs = ["GrC" "CF"];
scrsz = get(0,'screensize');
curd = mice.data{"R40"}{end};
dtimCb_tmplt = curd.dtimCb;
dtb = 0.005;
tmpxCb = curd.tmpxCb;
tmpxx = curd.tmpxx;
tpltCb = minix(abs(tmpxCb+2)):minix(abs(tmpxCb-2));
tpltCb_s = tmpxCb(tpltCb);
tpltB = minix(abs(tmpxx+2)):minix(abs(tmpxx-2));
tpltB_s = tmpxx(tpltB);
zerptCb_tmplt = minix(abs(tmpxCb));
tmpxCb_tmplt = mice.data{"R40"}{end}.tmpxCb;
zerptB = minix(abs(tmpxx));
winfnz = @(x,y)((round(x(1)/y):round(x(2)/y)));

algnlbls = ["Movement onset";"Movement midpoint";"Movement end";"Reward";"Reward"];
alignnms = table(algnlbls,RowNames=["startAlgn"; "midAlgn"; "endAlgn"; "rewAlgn"; "splitAlgn"]);
timax_str = @(x)(xlabel("Time relative to "+alignnms.algnlbls(x)+" (s)"));

siglabls = ["Granule cells","Climbing fibers"];
ncg = length(siglabls);

dlclabls = ["nose x","nose y","tongue x","tongue y","l forelimb x","l forelimb y","r forelimb x","r forelimb y"];
dlclablshalf = ["nose","tongue","l forelimb","r forelimb"];

% make split alignment
meanmvtime_all = zeros(0,1);
meanmvsttime_all = zeros(0,1);
for m = 1:nmice
    for d = 1:mice.ndays(m)
        curd = mice.data{m}{d};
        meanmvtime = median((curd.rewtimes-curd.midpt)*dtb);
        meanmvtime_all = cat(1,meanmvtime_all,meanmvtime);
        meanmvsttime = median((curd.rewtimes-curd.truestart)*dtb);
        meanmvsttime_all = cat(1,meanmvsttime_all,meanmvsttime);
    end
end
meanmvtime_all = mean(meanmvtime_all);
meanmvsttime_all = mean(meanmvsttime_all);
for m = 1:nmice
    for d = 1:mice.ndays(m)
        curd = mice.data{m}{d};
        zerptCb = minix(abs(curd.tmpxCb));
        mice.data{m}{d}.timesharefac = round(curd.frameCb(end)/curd.ntimCb);
tmp_CFsps_mv = permute(filtfilt([1 1],1,permute(double(curd.midAlgn.sp_CF),[3 1 2])),[2 3 1]);
tmp_CFsigs_mv = curd.midAlgn.sigFilt_CF;
tmp_GrCsigs_mv = curd.midAlgn.sigFilt_GrC;
tmp_CFsps_r = permute(filtfilt([1 1],1,permute(double(curd.rewAlgn.sp_CF) ...
    ,[3 1 2])),[2 3 1]);
tmp_CFsigs_r = curd.rewAlgn.sigFilt_CF;
tmp_GrCsigs_r = curd.rewAlgn.sigFilt_GrC;
tmp_CFsps = zeros(size(tmp_CFsps_mv));
tmp_CFsigs = zeros(size(tmp_CFsigs_mv));
tmp_GrCsigs = zeros(size(tmp_GrCsigs_mv));
tmpxCb = curd.tmpxCb;
meanmvtime = meanmvtime_all;
meanmvix = minix(abs(tmpxCb+meanmvtime));
meanmvixrev = zerptCb-minix(abs(tmpxCb+meanmvtime));
meanmvix_b = minix(abs(tmpxx+meanmvtime));
meanmvixrev_b = zerptB-minix(abs(tmpxx+meanmvtime));
brktime = meanmvtime-0.3;
brkpt = minix(abs(tmpxCb+brktime));
brkpto = minix(abs(tmpxCb-0.3));
brkpt_b = minix(abs(tmpxx+brktime));
brkpto_b = minix(abs(tmpxx-0.3));
win1 = 1:brkpt;
win1o = brkpto+(-(length(win1)-1):0);
win2 = brkpt:length(tmpxCb);
tmp_CFsigs(:,:,win1) = tmp_CFsigs_mv(:,:,win1o);
tmp_CFsps(:,:,win1) = tmp_CFsps_mv(:,:,win1o);
tmp_CFsigs(:,:,win2) = tmp_CFsigs_r(:,:,win2);
tmp_CFsps(:,:,win2) = tmp_CFsps_r(:,:,win2);
tmp_GrCsigs(:,:,win1) = tmp_GrCsigs_mv(:,:,win1o);
tmp_GrCsigs(:,:,win2) = tmp_GrCsigs_r(:,:,win2);
mice.data{m}{d}.splitAlgn.sigFilt_CF = tmp_CFsigs;
mice.data{m}{d}.splitAlgn.sp_CF = tmp_CFsps;
mice.data{m}{d}.splitAlgn.sigFilt_GrC = tmp_GrCsigs;
% mice.data{m}{d}.trialTyp = mice.data{m}{d}.mvtyp(mice.data{m}{d}.rewtimes);
win1 = 1:brkpt_b;
win1o = brkpto_b+(-(length(win1)-1):0);
win2 = brkpt_b:length(tmpxx);
tmp_lick_mv = curd.midAlgn.lick;
tmp_lick_r = curd.rewAlgn.lick;
mice.data{m}{d}.splitAlgn.lick(:,win1) = tmp_lick_mv(:,win1o);
mice.data{m}{d}.splitAlgn.lick(:,win2) = tmp_lick_r(:,win2);
    end
end

%%% within day concatenation

miceuse = 1:nmice;
nmiceuse = length(miceuse);
cfsiguse = "sp_CF";
sigs = ["sigFilt_GrC",cfsiguse];
epochs_use = ["novice","expert"]; % leave as is
nep = length(epochs_use);
tgrplbls = ["Rewarded","Omitted reward","rew half","rew half2"];
tgrpalgn = ["splitAlgn","splitAlgn","splitAlgn","splitAlgn"];

ntg = 4;
[grcsigs_all,cfsigs_all,grcsigs_all_ix,cfsigs_all_ix] = deal(cell(1,nep));
[grcsigs_all{:},cfsigs_all{:},grcsigs_all_ix{:},cfsigs_all_ix{:}] = deal(cell(ntg,1));
[grcsigs_all_s,cfsigs_all_s] = deal(cell(1,nep));
[grcsigs_all_s{:},cfsigs_all_s{:}] = deal(cell(ntg,0));
scnt=ones(length(epochs_use),1);
session_ids = {};
for mm = miceuse
    for ep = 1:nep
        epoch_use = epochs_use(ep);
        if epoch_use=="novice"
            % this condition limits to the FIRST day for "novice"
            daysuse_cur = find(mice.days{mm}.epoch==epoch_use & mice.days{mm}.trueDay==min(mice.days{mm}.trueDay));
        elseif ismember(epoch_use,"expert")
            % this condition omits "expert" days where same FOV was imaged >once
            daysuse_cur = find(mice.days{mm}.epoch==epoch_use & ...
                (isnan(mice.days{mm}.chronDay) | mice.days{mm}.chronDay==max(mice.days{mm}.chronDay)));
        end
    for dcur = 1:length(daysuse_cur)
        curd = mice.data{mm}{daysuse_cur(dcur)};
        tmpxCb = curd.tmpxCb;
        zerptCb = minix(abs(curd.tmpxCb));
        clear trgrps;
        
        % rewarded / unrewarded
        trgrps{1} = find(curd.rewarded & curd.mvlen>7 & curd.goodmvdir);
        trgrps{2} = find(~curd.rewarded & curd.mvlen>7 & curd.goodmvdir);
        
        % random halves of rewarded trials
        rp = shufflev(trgrps{1});
        trgrps{3} = rp(1:round(end/2));
        trgrps{4} = rp(round(end/2)+1:end);

        ntg = length(trgrps);
        grcsigs_cur = cell(ntg,1);
        cfsigs_cur = cell(ntg,1);
        for tg = 1:ntg
            if ~isempty(trgrps{tg})
            tmp = squeeze(mean(curd.(tgrpalgn(tg)).sigFilt_GrC(trgrps{tg},:,:),1));
            if length(tmpxCb)~=length(tmpxCb_tmplt),
                tmp = interp1(tmpxCb,tmp',tmpxCb_tmplt,"linear","extrap")';
                assert(~any(isnan(tmp(:))));
            end
            grcsigs_all_s{ep}{tg,scnt(ep)} = tmp;
            grcsigs_all{ep}{tg} = cat(1,grcsigs_all{ep}{tg},tmp);
            grcsigs_all_ix{ep}{tg} = cat(1,grcsigs_all_ix{ep}{tg},[repmat([mm,daysuse_cur(dcur)],curd.nIC_GrC,1),(1:curd.nIC_GrC)']);
            
            curtr_sigs = curd.(tgrpalgn(tg)).(cfsiguse)(trgrps{tg},:,:);
            nccur = size(curtr_sigs,2);
            if cfsiguse=="sp_CF"
                tmpspvec = reshape(permute(curtr_sigs,[3 1 2]),[],nccur);
                spkern = round(0.15 / curd.dtimCb);
                spkern_s = spkern*curd.dtimCb;
                tmpspvecf = filter(ones(spkern,1)/spkern/spkern_s,1,tmpspvec);
                tmpspvecf = zscore(tmpspvecf,[],1);
                tmpspvecf = reshape(tmpspvecf,size(curtr_sigs,3),size(curtr_sigs,1),nccur);
                curtr_sigs = permute(tmpspvecf,[2 3 1]);
            end

            tmp = squeeze(mean(curtr_sigs,1));
            if length(tmpxCb)~=length(tmpxCb_tmplt),
                tmp = interp1(tmpxCb,tmp',tmpxCb_tmplt,"linear","extrap")';
                assert(~any(isnan(tmp(:))))
            end
            cfsigs_all{ep}{tg} = cat(1,cfsigs_all{ep}{tg},tmp);
            cfsigs_all_ix{ep}{tg} = cat(1,cfsigs_all_ix{ep}{tg},[repmat([mm,daysuse_cur(dcur)],curd.nIC_CF,1),(1:curd.nIC_CF)']);
            cfsigs_all_s{ep}{tg,scnt(ep)} = tmp;
            end
        end
        session_ids{ep}(scnt(ep),:) = [mice.name{mm} daysuse_cur(dcur) mice.days{mm}{daysuse_cur(dcur),["date" "epoch"]}];
        scnt(ep)=scnt(ep)+1;
        end
    end
end
clear  tg tgrpalgn tgvec  


%% single cell significance
miceuse = 1:nmice;
cfsiguse = "sp_CF";
sigs = ["sigFilt_GrC",cfsiguse];
epochs_use = ["novice","expert"]; % leave as is
neps= length(epochs_use);

% collect cell quantifications and pvals
quantheaders = ["name","algnuse","sigsUse","avgWins","winDiff","trialGrps","results","sinfo"];
quantTyps = cell2table(cell(0,length(quantheaders)),VariableNames=quantheaders);

quantTyps = cat(1,quantTyps,{"delaySuppCF","splitAlgn",{["sp_CF"]},{[-0.75 0]},{missing},{missing},cell(1,1),cell(1,1)});
quantTyps = cat(1,quantTyps,{"anticipatoryCF","splitAlgn",{["sp_CF"]},{[-0.3 -0.03]},{missing},{missing},cell(1,1),cell(1,1)});
quantTyps = cat(1,quantTyps,{"anticipatoryGrC","splitAlgn",{["sigFilt_GrC"]},{[-0.3 -0.03]},{missing},{missing},cell(1,1),cell(1,1)});

quantTyps = cat(1,quantTyps,{"anticipatoryRiseGrC","splitAlgn",{["sigFilt_GrC"]},{[-0.3 -0.03; -1.5 -1.3]},{[1 2]},{missing},cell(1,1),cell(1,1)});
quantTyps = cat(1,quantTyps,{"GrCrewReduction","splitAlgn",{["sigFilt_GrC"]},{[-0.3 -0.03; 0.3 0.5]},{[1 2]},{[1 1]},cell(1,1),cell(1,1)});

quantTyps = cat(1,quantTyps,{"CFrewElevation","splitAlgn",{["sp_CF"]},{[0 0.2; -0.3 -0.03]},{[1 2]},{[1 1]},cell(1,1),cell(1,1)});
quantTyps = cat(1,quantTyps,{"GrCrewElevation","splitAlgn",{["sigFilt_GrC"]},{[0 0.2; -0.3 -0.03]},{[1 2]},{[1 1]},cell(1,1),cell(1,1)});
quantTyps = cat(1,quantTyps,{"CFrewResponse","splitAlgn",{["sp_CF"]},{[0 0.2]},{missing},{[1]},cell(1,1),cell(1,1)});

quantTyps = cat(1,quantTyps,{"pre-movementGrC","startAlgn",{["sigFilt_GrC"]},{[-0.5 -0.25]},{missing},{[1]},cell(1,1),cell(1,1)});
quantTyps = cat(1,quantTyps,{"pre-movementCF","startAlgn",{["sp_CF"]},{[-0.5 -0.25]},{missing},{[1]},cell(1,1),cell(1,1)});
quantTyps = cat(1,quantTyps,{"GrCmoveElevation","startAlgn",{["sigFilt_GrC"]},{[0 0.2; -0.3 -0.03]},{[1,2]},{[1]},cell(1,1),cell(1,1)});
quantTyps = cat(1,quantTyps,{"CFmoveElevation","startAlgn",{["sp_CF"]},{[0 0.2; -0.3 -0.03]},{[1 2]},{[1]},cell(1,1),cell(1,1)});

quantTyps.Properties.RowNames=quantTyps.name;

clear session_ids;
sessioncount= ones(length(epochs_use),1);
session_ids = cell(length(epochs_use),1);
for ep = 1:neps
    epoch_use = epochs_use(ep);
for mm = miceuse
    if epoch_use=="novice"
        % this condition limits to the FIRST day for "novice"
        daysuse_cur = find(mice.days{mm}.epoch==epoch_use & mice.days{mm}.trueDay==min(mice.days{mm}.trueDay));
    elseif ismember(epoch_use,"expert")
        % this condition omits "expert" days where same FOV was imaged >once
        daysuse_cur = find(mice.days{mm}.epoch==epoch_use & ...
            (isnan(mice.days{mm}.chronDay) | mice.days{mm}.chronDay==max(mice.days{mm}.chronDay)));
    end

    for dd = 1:length(daysuse_cur)
        dcur = daysuse_cur(dd);
        curd = mice.data{mm}{dcur};
        cnums = [curd.nIC_GrC curd.nIC_CF];
        zerptCb = minix(abs(curd.tmpxCb));
        dtimCb = curd.dtimCb;
        tmpxCb = curd.tmpxCb;
        clear trgrps;
        
        % rewarded / unrewarded
        trgrps{1} = find(curd.rewarded & curd.mvlen>7 & curd.goodmvdir);
        trgrps{2} = find(~curd.rewarded & curd.mvlen>7 & curd.goodmvdir);
        tgrplbls = ["Rewarded","Omitted reward"];
        tgrpalgn = ["rewAlgn","rewAlgn"];

        ntg = length(trgrps);
        grcsigs_cur = cell(ntg,1);
        cfsigs_cur = cell(ntg,1);
        
        ntr_cur = curd.nc;
        for q = 1:size(quantTyps,1)
            for cg = 1:length(quantTyps.sigsUse{q})
                curtr_sigs = curd.(quantTyps.algnuse{q}).(quantTyps.sigsUse{q}(cg));
                nccur = size(curtr_sigs,2);
                
                if quantTyps.sigsUse{q}(cg)=="sp_CF"
                    tmpspvec = reshape(permute(curtr_sigs,[3 1 2]),[],nccur);
                    spkern = round(0.2 / dtimCb);
                    spkern_s = spkern*dtimCb;
                    tmpspvecf = filtfilt(ones(spkern,1)/spkern,1,double(tmpspvec));
                    tmpspvecf = zscore(tmpspvecf,[],1);
                    tmpspvecf = reshape(tmpspvecf,size(curtr_sigs,3),size(curtr_sigs,1),nccur);
                    tmpspvecf = permute(tmpspvecf,[2 3 1]);
                    curtr_sigs = tmpspvecf;
                end

                if ~ismissing(quantTyps.trialGrps{q}), 
                    curtrgrp = trgrps{quantTyps.trialGrps{q}(1)};
                else curtrgrp = 1:ntr_cur;
                end
                win1 = winfn(quantTyps.avgWins{q}(1,:),tmpxCb);
                tr_sig = {};
                mnval = squeeze(mean(curtr_sigs(curtrgrp,:,win1),3));
                pval = arrayfun(@(x)signrank(mnval(:,x)),1:nccur);
                mnval = squeeze(mean(mnval,1));
                
                if quantTyps.sigsUse{q}(cg)=="sigFilt_GrC"
                tmpwin = winfn([-0.75 0],tmpxCb);
                varval = squeeze(std(curtr_sigs(curtrgrp,:,tmpwin),[],1));
                mnval2 = squeeze(mean(curtr_sigs(curtrgrp,:,tmpwin),1));
                mnval2 = abs(mnval2);
                vsrval = mean(varval./mnval2,2);
                tmpwin2 = winfn([-2 2],tmpxCb);
                curhist = squeeze(mean(curtr_sigs(curtrgrp,:,tmpwin2),1));
                curhist = abs(curhist);
                curent = [];
                for i = 1:size(curhist,1)
                    [~,curix] = max(curhist(i,:));
                    tmpix = curix+(-1:1); tmpix(tmpix<1) = []; tmpix(tmpix>length(tmpwin2))=length(tmpwin2);
                    tmpcurtr_sigs = curtr_sigs(curtrgrp,i,tmpwin2);
                    tmpcurtr_sigs = abs(tmpcurtr_sigs);
                    tmp = mean(tmpcurtr_sigs(:,tmpix),2);
                    curent(end+1,1) = squeeze(mean(std(curtr_sigs(curtrgrp,i,tmpwin2),[],3)./tmp,1))';
                end

                tr_sig{1} = [mnval' pval' vsrval curent];
                else
                    tr_sig{1} = [mnval' pval'];
                end
                if any(~ismissing(quantTyps.winDiff{q}))
                    winuse = quantTyps.winDiff{q}(2);
                    win2 = winfn(quantTyps.avgWins{q}(winuse,:),tmpxCb);
                    if length(quantTyps.trialGrps{q})>1
                        curtrgrp2 = trgrps{quantTyps.trialGrps{q}(2)};
                    else, curtrgrp2 = curtrgrp;
                    end
                    if ~isempty(curtrgrp2)
                    mnval1 = squeeze(mean(curtr_sigs(curtrgrp,:,win1),3));
                    mnval2 = squeeze(mean(curtr_sigs(curtrgrp2,:,win2),3));
                    if size(mnval1,1)==size(mnval2,1)
                    pval = arrayfun(@(x)signrank(mnval1(:,x),mnval2(:,x)),1:nccur);
                    mnval = squeeze(mean(mnval1-mnval2,1));
                    else
                        pval = arrayfun(@(x)ranksum(mnval1(:,x),mnval2(:,x)),1:nccur);
                        mnval = squeeze(mean(mnval1)-mean(mnval2,1));
                    end
                    if quantTyps.sigsUse{q}(cg)=="sigFilt_GrC"
                    cur_reslt = [mnval' pval' vsrval curent];
                    else
                        cur_reslt = [mnval' pval'];
                    end
                    else
                        cur_reslt = [];
                    end
                else
                    cur_reslt = tr_sig{1};
                end
                quantTyps.results{q}{ep}{cg}{sessioncount(ep)} = cur_reslt;
                assert(size(cur_reslt,1)==nccur);
                quantTyps.sinfo{q}{ep}{cg}{sessioncount(ep)} = [repmat([mm dcur],nccur,1) (1:nccur)']; 
            end
        end
        session_ids{ep}(sessioncount(ep),:) = [mice.name{mm} dcur mice.days{mm}{dcur,["date" "epoch"]}];
        sessioncount(ep) = sessioncount(ep)+1;
    end
end
end
% identify sessions and trials with functioning lick sensor
algnuse = "rewAlgn";
licksensor_valid = cell(nmice,1);
for muse=string(mice.name)'
muse_ix = find(mice.name==muse);
licksensor_valid{muse_ix} = strings(size(mice.days{muse},1),1);
for d_use = 1:mice.ndays(muse_ix)
curdate = mice.days{muse}.date(d_use);
curd = mice.data{muse}{d_use};
rewtrs = find(curd.rewarded);
tmpsig = mean(curd.(algnuse).lick(rewtrs,:),2);
tmp = find(tmpsig>0.9 | tmpsig<0.01);
licksensor_valid{muse_ix}(d_use) = "y";
tmprat = length(tmp)/length(rewtrs);
if  tmprat>0.05, licksensor_valid{muse_ix}(d_use) = "p"; end
if  tmprat>0.9, licksensor_valid{muse_ix}(d_use) = "n"; 
    mice.data{muse_ix}{d_use}.valid_lick_trials = [];
end
fprintf("%s %s lick sensor bad trial ratio: %1.2f, valid: %s\n",mice.name(muse_ix),curdate,tmprat,licksensor_valid{muse_ix}(d_use));
tmpwin2 = winfn([-2 2],curd.tmpxx);
if licksensor_valid{muse_ix}(d_use)~="n"
    test = filter(ones(300,1)/300,1,curd.rewAlgn.lick(:,tmpwin2)')';
    % valid_lick = find(max(test,[],2)<0.95);
    valid_lick = find(max(test,[],2)<0.9);
    valid_lick_trials = intersect(valid_lick, find(curd.goodmvdir));
    mice.data{muse_ix}{d_use}.valid_lick_trials = valid_lick_trials;
    if length(valid_lick_trials)/size(test,1)<=0.2
    % if length(valid_lick_trials)/size(test,1)<0.7
        licksensor_valid{muse_ix}(d_use) = "n"; 
        mice.data{muse_ix}{d_use}.valid_lick_trials=[];
    end
    if ismember([muse_ix d_use],[4 6; 4 7; 5 6]) % additional datasets with corrupted sensor
        licksensor_valid{muse_ix}(d_use) = "n"; 
        mice.data{muse_ix}{d_use}.valid_lick_trials=[];
    end
    fprintf("valid trial frac: %2.1g\n",length(valid_lick_trials)/size(test,1));
end
end
end
clear val_len;
val_len_ep = cell(4,1);
for i = 1:nmice
    for d = 1:mice.ndays(i)
        curd = mice.data{i}{d};
        val_len{i}(d) = length(curd.valid_lick_trials);
        ncur = curd.nc;
        curep = find(eplist==mice.days{i}.epoch(d));
        if mice.days{i}.chronDay(d)==min(mice.days{i}.chronDay), curep = 1; end
        if val_len{i}(d)<10 && licksensor_valid{i}(d)~="n", 
            fprintf("m=%d d=%d val=%d tot=%d frac=%2.1g\n",...
                i,d,val_len{i}(d),ncur,double(val_len{i}(d))/ncur);
            licksensor_valid{muse_ix}(d_use) = "n"; 
            mice.data{muse_ix}{d_use}.valid_lick_trials=[];
        end
        val_len_ep{curep} = cat(1,val_len_ep{curep},val_len{i}(d));
    end
end
clear lickstats;
lickstats.tottrials = sum(cellfun(@(x)sum(x),val_len));
lickstats.trialsPerEp = cellfun(@(x)sum(x),val_len_ep)';
lickstats.totSessions = sum(cellfun(@(x)sum(x~=0),val_len));
lickstats.sessionsPerEp = cellfun(@(x)sum(x~=0),val_len_ep)';
lickstats.trialsPerMouse = cellfun(@(x)sum(x),val_len);
lickstats.nMice = sum(lickstats.trialsPerMouse~=0);
disp(lickstats)
%% within day plotting, 3A-D
ntg = 2;
plotheaders = ["name","algnuse","sortbycrit","showomit","sortwins","celltyps","results"];
plottyps = cell2table(cell(0,length(plotheaders)),VariableNames=plotheaders);
algnuse = "splitAlgn";
% FIG 3a,b
plottyps = cat(1,plottyps,{"3A,B sorted by peak time",algnuse,"max_ix",false,{missing},{["GrCs","CFs"]},cell(2,ntg)});
% FIG 3c
plottyps = cat(1,plottyps,{"3C pre-reward response",algnuse,"mean",false,{[-0.35 -0.05]},{"GrCs"},cell(2,ntg)});
% FIG 3d
plottyps = cat(1,plottyps,{"3D post-reward response",algnuse,"mean",false,{[0.01 0.25]},{"CFs"},cell(2,ntg)});

for p = 1:size(plottyps,1)
sigsuse = {};
if ismember("GrCs",plottyps.celltyps{p}), sigsuse=cat(2,sigsuse,{grcsigs_all}); end
if ismember("CFs",plottyps.celltyps{p}), sigsuse=cat(2,sigsuse,{cfsigs_all}); end
ep = 2;
figure(pos=[610    82   708   444*length(sigsuse)])

crng1 = [-0.4 0.6];
if cfsiguse=="sigFilt_CF", crng2 = [-0.5 1.3];
elseif cfsiguse=="sp_CF", crng2 = [-0.3 0.5]; end
crng = {crng1,crng2};
toquant = cell(2,ntg);
ix_neg_pos = cell(ncg,ntg,2);
tgvec = [1]; if plottyps.showomit(p), tgvec = [tgvec 2]; end
for s = 1:size(sigsuse{1}{ep},2)
for cg = 1:length(sigsuse)
    curtgvec=tgvec;
    for tg = curtgvec
        ax_x = tg;
        subaxis(length(sigsuse),ntg,ax_x,cg);
        cur = sigsuse{cg}{ep}{tg,s};
        ncur = size(cur,1);
        if ~all(ismissing(plottyps.sortwins{p}))
            clear sortrng;
            for i = 1:size(plottyps.sortwins{p},1)
                sortrng{i} = winfn(plottyps.sortwins{p}(i,:),tmpxCb_tmplt);
            end
        tmpsiguse = mean(sigsuse{cg}{ep}{tg,s}(:,sortrng{1}),2);
        end
        if plottyps.sortbycrit(p)=="max_ix"% sort by max time
        [~,sortcrit_cur] = max(cur(:,tpltCb,:),[],2);
        elseif plottyps.sortbycrit(p)=="mean"
        sortcrit_cur = tmpsiguse;
        toquant{cg,tg} = sortcrit_cur;
%             elseif sortbycrit==rew
        end

        [~,tmpix] = sort(sortcrit_cur);
        ix_neg_pos{cg,tg,1} = find(sortcrit_cur<-0.);
        ix_neg_pos{cg,tg,2} = find(sortcrit_cur>0.);
        
        imagesc(tpltCb_s,1:ncur,cur(tmpix,tpltCb),crng{cg});
        hold on; plot([0 0],ylim,'k');
        hold on; plot(-[meanmvtime_all meanmvtime_all],ylim,'k');
        yl = ylim;
        hold on; plot(-brktime,yl(2),'r.','markersize',5)
        if cg==length(sigsuse), timax_str(plottyps.algnuse(p)); end
        if ax_x == 1, 
            ylabel(plottyps.celltyps{p}(cg)+sprintf("%d cells, %d sessions %d mice",...
                ncur,size(session_ids{ep},1),length(unique(session_ids{ep}(:,1))))); 
        end
        if cg == 1, title(tgrplbls(tg)); end
        if cg==1 && ax_x==max(curtgvec), h=gca; despos = h.Position; colorbar; h.Position=despos; end
        if cg==2 && ax_x==max(curtgvec) && cfsiguse=="sp_CF", h=gca; despos = h.Position; colorbar; h.Position=despos; end
        if ax_x==1
            cur_ax = gca;
            axpos = cur_ax.Position;
            if ~ismissing(plottyps.sortwins{p})
            annotation('line',(axpos(1)+axpos(3)/2)+(plottyps.sortwins{p}(1,:))/(diff(tpltCb_s([1 end])))*axpos(3),[axpos(2) axpos(2)]-0.003,col='r',linew=2);
            annotation('line',(axpos(1)+axpos(3)/2)+(plottyps.sortwins{p}(1,:))/(diff(tpltCb_s([1 end])))*axpos(3),[1 1]*(axpos(2)+axpos(4)+0.003),col='r',linew=2);
            end
        end
    end
end
plottyps.results{p} = toquant;
title(plottyps.name(p));
pause; if s<size(sigsuse{1}{ep},2), clf; end
end
end
%% display single cell quantifications. 3E-L, 4A-D, S2N-S
clear fp;
pausebw = false;
fp(1).name = "3E";
fp(1).pthresh = 0.05;
fp(1).magthr = 0.2;
fp(1).useomit = true;
fp(1).sigtyps = ["GrCrewElevation";"CFrewElevation"];
fp(1).epuse = "expert";
fp(1).signchks = ["pos","neg"];
fp(1).comptyp = "cell";
fp(1).plottyps = 1; % 1=prevalence and 2=magnitude
fp(1).showhist = true; % 
fp(1).showtraces = false; %

fp(2) = fp(1);
fp(2).name = "3F";
fp(2).sigtyps = ["anticipatoryGrC";"anticipatoryCF"];

fp(3) = fp(1);
fp(3).name = "3G-M";
fp(3).sigtyps = ["anticipatoryRiseGrC","GrCrewReduction";"CFrewResponse","CFrewElevation"];
fp(3).signchks = "pos";
fp(3).pthresh = 1;
fp(3).magthr = 0.1;
fp(3).plottyps = 1;
fp(3).showhist = false;
fp(3).showtraces = true;

fp(4) = fp(1);
fp(4).name = "S2P-pre-move";
fp(4).sigtyps = ["pre-movementGrC";"pre-movementCF"];
fp(4).signchks = "pos";
fp(4).plottyps = [];
fp(4).showhist=false;

fp(5) = fp(4);
fp(5).name = "S2P-post-move";
fp(5).sigtyps = ["GrCmoveElevation";"CFmoveElevation"];

fp(6) = fp(1);
fp(6).name = "4a,c";
fp(6).pthresh = 1;
fp(6).magthr = 0;
fp(6).useomit = false;
fp(6).sigtyps = ["anticipatoryRiseGrC","GrCrewReduction"];
fp(6).comptyp = "epoch";
fp(6).epuse = ["novice","expert"];
fp(6).signchks = "pos";
fp(6).plottyps = 2;
fp(6).showhist = true; % 
fp(6).showtraces = true; %

fp(7) = fp(6);
fp(7).name = "4b,d";
fp(7).sigtyps = ["CFrewResponse"];

fp(8) = fp(6);
fp(8).name = "4c-inset";
fp(8).pthresh = 0.05;
fp(8).magthr = 0.2;
fp(8).showtraces = false;
fp(8).plottyps = 1;

fp(9) = fp(8);
fp(9).name = "4d-inset";
fp(9).sigtyps = ["CFrewResponse"];

fp(10) = fp(6);
fp(10).name = "S2Q-cfsupp";
fp(10).sigtyps = ["delaySuppCF"];
fp(10).signchks = "neg";
fp(10).showtraces = false; %

fp(11) = fp(8);
fp(11).name = "S2Q-cfsupp-inset";
fp(11).signchks = "neg";
fp(11).sigtyps = ["delaySuppCF"];

for f = 1:length(fp)
[~,epsusen] = ismember(fp(f).epuse,epochs_use);

thresh = fp(f).magthr;
pthresh =fp(f).pthresh;
lp = 1; hp = 99;
inc = 0.1;
clear h hinf;
if fp(f).useomit
    trtypes = [1 2];
else
    trtypes = 1;
end
for e = 1:length(epsusen)
    ep = epsusen(e);
    for c = 1:size(fp(f).sigtyps,1)
        for cond = 1:size(fp(f).sigtyps,2)
            h{e}{c,cond} = cat(1,quantTyps.results{fp(f).sigtyps(c,cond)}{ep}{1}{:});
            hinf{e}{c,cond} = cat(1,quantTyps.sinfo{fp(f).sigtyps(c,cond)}{ep}{1}{:});
        end
    end
end
nc = sum(double([contains(strjoin(fp(f).sigtyps(:)),"GrC"),contains(strjoin(fp(f).sigtyps(:)),"CF")]));
eps = repmat(epsusen,nc,1); ctyps = repmat(1:nc,length(fp(f).epuse),1)';
if fp(f).comptyp == "cell"
    eps = eps'; 
    ctyps = ctyps';
end
if fp(f).comptyp == "cell"
for i = 1:size(eps,1)
    e = eps(i,:); c = ctyps(i,:);
    [~,e] = ismember(e,epsusen);
    sig2 = h{e(2)}{c(2),1}(:,1);
    sig1 = h{e(1)}{c(1),1}(:,1);
    if fp(f).showhist
    histogram_cell({sig1,sig2},cdf=true)
    [~,s]=ranksum_disp(sig1,sig2,ks=true);
    titstr = sprintf("%s %s %s",fp(f).epuse(e(1))+fp(f).sigtyps(c(1)),fp(f).epuse(e(2))+fp(f).sigtyps(c(2)),s);
    titstr = [titstr;fp(f).name];
    title(titstr);
    yl = ylim;
    hold on; plot([0 0],yl,'k')
    xlabel(fp(f).sigtyps(1)+" (zsc)");
    ylabel("Fraction of cells");
    legend(fp(f).epuse(e(1))+fp(f).sigtyps(c(1)),fp(f).epuse(e(2))+fp(f).sigtyps(c(2)),box=0)
    box off;
    end
end
end
proportions = {};
proportions_persess = {};
proportions_ix = {};
sinf_selected = {};
allinf = {};
sigsigs = {};
for e = 1:length(epsusen)
for c = 1:nc
    curconds = fp(f).sigtyps(c,:);
    for s = 1:length(fp(f).signchks)
        condreslts = true(size(h{e}{c,cond}(:,1)));
        cursign = fp(f).signchks(s);
        if cursign=="pos",curfn = @gt; curfac = 1; else, curfn=@lt; curfac = -1; end
        for cond = 1:length(curconds)
            curthresh = thresh;
            curreslt = h{e}{c,cond}(:,2)<pthresh & curfn(h{e}{c,cond}(:,1),curthresh*curfac);
            condreslts = and(condreslts,curreslt);
        end
        proportions{e,c,s} = condreslts;
        proportions_ix{e,c,s} = find(condreslts);
        sinf_selected{e,c,s} = hinf{e}{c,cond}(proportions_ix{e,c,s},:);
        sigsigs{e,c,s} = h{e}{c,1}(proportions_ix{e,c,s},1);
        allinf{e,c,s} = hinf{e}{c,cond};
        allsess = unique(allinf{e,c,s}(:,1:2),'rows');
        proportions_persess{e,c,s} = zeros(length(allsess),1);
        for sess = 1:length(allsess)
            cursess = find(all(allinf{e,c,s}(:,1:2)==allsess(sess,:),2));
            proportions_persess{e,c,s}(sess) = mean(proportions{e,c,s}(cursess,:));
        end

    end
end
end
if fp(f).name=="3G-M"
    cursinf = sinf_selected(1,1,:); cursinf = squeeze(cursinf(:));
    ixantic_ss = cursinf; 
end
tmppropsmag = {proportions_persess,sigsigs};
if fp(f).name=="S2P-pre-move", premvprops=proportions_persess; end
if fp(f).name=="S2P-post-move", postmvprops=proportions_persess; end
if contains(fp(f).sigtyps(1,1),"GrC")
tmpsigsuse  = {grcsigs_all,cfsigs_all};
else
tmpsigsuse  = {cfsigs_all};
end
for fcur = 1:size(eps,1)
    ep = unique(eps(fcur,:));
    trueep = ep;
    [~,ep] = ismember(ep,epsusen);
    curc = unique(ctyps(fcur,:));
    for curplt = fp(f).plottyps
    cursig = tmppropsmag{curplt};
    cursinf = sinf_selected(ep,curc,:); cursinf = squeeze(cursinf(:));
    cursinf_all = allinf(ep,curc,:); cursinf_all = squeeze(cursinf_all(:));
    tmpcel = cursig(ep,curc,:); tmpcel = squeeze(tmpcel(:));
    if curplt==1
    barwitherr_stde(tmpcel,dots=true);
    else
        histogram_cell(tmpcel,cdf=true); f2=gcf;
    end
    ranksum_disp(tmpcel{1},tmpcel{2});
    if length(tmpcel)>2
    ranksum_disp(tmpcel{3},tmpcel{4});
    end
    hold on;
    tmpsignstr = repelem(fp(f).signchks,1,length(tmpcel)/length(fp(f).signchks));
    tmpcstr = repmat(fp(f).sigtyps(ctyps(fcur,curc)),1,length(tmpcel)/length(fp(f).sigtyps(ctyps(fcur,curc))));
    [~,tmpep] = ismember(eps(fcur,:),epsusen);
    epsstr = fp(f).epuse(tmpep);
    epsstr = repmat(epsstr,1,length(tmpcel)/length(epsstr));
    lablstr = arrayfun(@(x)sprintf("%s-%s-%s",tmpsignstr(x),...
        tmpcstr(x),epsstr(x)),1:length(tmpcel));
    if curplt==1
    set(gca,"Units","Pixels")
    despos=get(gca,"position");
    xticklabels(lablstr);
    xlim([0.5 length(tmpcel)+0.5])
    fprintf("Proportion significance differences:\n")
    ylabel("% of cells")
    elseif curplt==2
    figure(f2);
    legend(lablstr);
    fprintf("Magnitude significance differences:\n")
    ylabel("signal (zsc)")
    end
    titstr = "";
    for i = 1:length(tmpcel)/2
        tmpstr = sprintf("%s-%s-%s\nvs %s-%s-%s\n ",tmpsignstr((i-1)*2+1),tmpcstr((i-1)*2+1),...
            epsstr((i-1)*2+1),tmpsignstr((i-1)*2+2),tmpcstr((i-1)*2+2),...
            epsstr((i-1)*2+2));
        [~,tmpstr2]=ranksum_disp(tmpcel((i-1)*2+(1:2)));
        titstr = titstr + tmpstr+tmpstr2;
        uniq_s_1 = unique(cursinf{(i-1)*2+1}(:,1:2),'rows');
        uniq_s_2 = unique(cursinf{(i-1)*2+2}(:,1:2),'rows');
        uniq_m_1 = unique(cursinf{(i-1)*2+1}(:,1),'rows');
        uniq_m_2 = unique(cursinf{(i-1)*2+2}(:,1),'rows');
        uniq_s_a = unique(cursinf_all{(i-1)*2+1}(:,1:2),'rows');
        uniq_m_a = unique(cursinf_all{(i-1)*2+1}(:,1),'rows');
        uniq_s_a_2 = unique(cursinf_all{(i-1)*2+2}(:,1:2),'rows');
        uniq_m_a_2 = unique(cursinf_all{(i-1)*2+2}(:,1),'rows');
        sessstr = sprintf("nsess: %d/%d and %d/%d\n",length(uniq_s_1),length(uniq_s_a),length(uniq_s_2),length(uniq_s_a_2));
        micestr = sprintf("nmice: %d/%d and %d/%d\n",length(uniq_m_1),length(uniq_m_a),length(uniq_m_2),length(uniq_m_a_2));
        disp(sessstr);
        disp(micestr);
    end
    titstr = [titstr;fp(f).name];
    title(titstr);
    if curplt==2, figure(f2); title(titstr); 
    elseif curplt==1, restoreAxSz(despos); end
    end
    if 1
    tmpxUse = winfn([-2 2],tmpxCb_tmplt);
    tmpxUse_s = tmpxCb_tmplt(tmpxUse);
    clear tmpsig2 tmpsig1
    for c = 1:length(curc)
    tmpsigs = tmpsigsuse{curc(c)};
    % for epoch comparisons, tmpsig2=rewarded expert,tmpsig1=rewarded novice
    % for expert comparison, tmpsig2=omitted expert, tmpsig1=rewarded expert
    tmpsig2{c} = tmpsigs{trueep(end)}{trtypes(end)}(proportions_ix{ep(end),curc(c),1},tmpxUse);
    tmpsig1{c} = tmpsigs{trueep(1)}{trtypes(1)}(proportions_ix{ep(1),curc(c),1},tmpxUse);
    end
    end
end
%
if 0
tmpwin = winfn(tmpxUse_s,[-0.35,-0.05]);
barwitherr_stde({mean(tmpsig1{c}(:,tmpwin),2),mean(tmpsig2{c}(:,tmpwin),2)})
end

if fp(f).showtraces
% cell averages from significance tests above
climsCF = [-0.25 0.5];
climsGrC = [-0.5 1];
figure; 
errorbar_shadeSEM(tmpxUse_s,tmpsig2{1}',[ 1 0 0]);
hold on; errorbar_shadeSEM(tmpxUse_s,tmpsig1{1}',[0 0 1]);
titstr = sprintf("n1:%d, n2:%d",size(tmpsig1{1},1),size(tmpsig2{1},1));
titstr = [titstr;fp(f).name];
title(titstr);
axis tight;
hold on; plot([0 0],ylim,'k-');
hold on; plot(-[meanmvtime_all meanmvtime_all],ylim,'k-');
if length(tmpsig1)>1
figure;
errorbar_shadeSEM(tmpxUse_s,tmpsig1{2}',[0 0 1]);
hold on; 
errorbar_shadeSEM(tmpxUse_s,tmpsig2{2}',[1 0 0]);
axis tight;
titstr = sprintf("n1:%d, n2:%d",size(tmpsig1{2},1),size(tmpsig2{2},1));
titstr = [titstr;fp(f).name];
title(titstr);
hold on; plot([0 0],ylim,'k-');
hold on; plot(-[meanmvtime_all meanmvtime_all],ylim,'k-');
end

if fp(f).useomit
tmpgrcrew_f = filtfilt(ones(3,1)/3,1,double(tmpsig1{1})')';
tmpgrcomit_f = filtfilt(ones(3,1)/3,1,double(tmpsig2{1}'))';
tmpmaxix = minix(abs(tmpxUse_s));
tmpmax = tmpgrcrew_f(:,tmpmaxix);
tmpgrcrew_halfmax = zeros(size(tmpmax));
tmpgrcomit_halfmax = zeros(size(tmpmax));
for i = 1:size(tmpgrcrew_f,1)
    junk = find(tmpgrcrew_f(i,tmpmaxix:end)<0.5*tmpmax(i),1);
    if isempty(junk)
        tmpgrcrew_halfmax(i) = length(tmpxUse); 
    else
        tmpgrcrew_halfmax(i) = junk+tmpmaxix-1;
    end
    junk = find(tmpgrcomit_f(i,tmpmaxix:end)<0.5*tmpmax(i),1);
    if isempty(junk)
        tmpgrcomit_halfmax(i) = length(tmpxUse); 
    else
        tmpgrcomit_halfmax(i) = junk+tmpmaxix-1;
    end
end
if fp(f).name=="3G-M"
tmpcfrew_f = tmpsig1{2};
tmpcfomit_f = tmpsig2{2};
tmpminix = minix(abs(tmpxUse_s+0.2));
tmpminwin = winfn([-0.35 0.05],tmpxUse_s);
tmpmin = mean(tmpcfrew_f(:,tmpminwin),2);
tmpcfrew_halfmin = arrayfun(@(x)find(tmpcfrew_f(x,tmpmaxix:end)>0.5*tmpmin(x),1),1:size(tmpcfrew_f,1))+tmpmaxix-1;
tmpcfomit_halfmin = arrayfun(@(x)find(tmpcfomit_f(x,tmpmaxix:end)>0.5*tmpmin(x),1),1:size(tmpcfomit_f,1))+tmpmaxix-1;
fprintf("\nGrC offtimes rew vs omit %s\n",signrank_disp(tmpxUse_s(tmpgrcrew_halfmax),tmpxUse_s(tmpgrcomit_halfmax)));
[meanmax,meanmaxix] = max(mean(tmpsig1{1}),[],2);
meanhalfmax_rew = tmpxUse_s(find(mean(tmpsig1{1}(:,meanmaxix:end))<0.5*meanmax,1)+meanmaxix-1);
meanhalfmax_omit = tmpxUse_s(find(mean(tmpsig2{1}(:,meanmaxix:end))<0.5*meanmax,1)+meanmaxix-1);
fprintf("\nmean GrC off-time rewarded/omitted %2.1g/%2.1g\n",meanhalfmax_rew,meanhalfmax_omit)
fprintf("\nCF ontime rew vs omit %s\n",signrank_disp(tmpxUse_s(tmpcfrew_halfmin),tmpxUse_s(tmpcfomit_halfmin)))
meanmin = mean(mean(tmpsig1{2}(:,tmpminwin),2));
meanhalfmax_rew = tmpxUse_s(find(mean(tmpsig1{2}(:,tmpmaxix:end))>0.5*meanmin,1)+tmpmaxix-1);
meanhalfmax_omit = tmpxUse_s(find(mean(tmpsig2{2}(:,tmpmaxix:end))>0.5*meanmin,1)+tmpmaxix-1);
fprintf("\nmean CF ontime rewarded/omit %2.1g/%2.1g\n",meanhalfmax_rew,meanhalfmax_omit)

ngc = size(tmpsig1{1},1);
delay_activ_cent = zeros(ngc,1);
ix_1sdel = winfn(tmpxUse_s,[-1.1 0]);
ix_premv = winfn(tmpxUse_s,[-1.9 -1.4]);
for c= 1:ngc
    delay_activ_indeces = find(tmpsig1{1}(c,ix_1sdel)-mean(tmpsig1{1}(c,ix_premv),2)>0);
    tmp4 = tmpsig1{1}(c,ix_1sdel(delay_activ_indeces));
    delay_activ_cent(c) = sum(tmp4./sum(tmp4).*tmpxUse_s(ix_1sdel(delay_activ_indeces)));
end
delay_activ_cent(delay_activ_cent==0)=-inf;
[~,idxsrt] = sort(delay_activ_cent); rs = idxsrt;
figure; subaxis(1,2,1); imagesc(tmpxUse_s,[],tmpsig1{1}(rs,:),climsGrC); 
hold on; plot([0 0],ylim,'k-');
hold on; plot(-[meanmvtime_all meanmvtime_all],ylim,'k-');
subaxis(1,2,2); imagesc(tmpxUse_s,[],tmpsig2{1}(rs,:),climsGrC);
h=gca; despos = h.Position; colorbar; h.Position=despos; 
hold on; plot([0 0],ylim,'k-');
hold on; plot(-[meanmvtime_all meanmvtime_all],ylim,'k-');
title(fp(f).name)

ncf = size(tmpsig1{2},1); 
rs_cf = randperm(ncf);
figure; subaxis(1,2,1); imagesc(tmpxUse_s,[],tmpsig1{2}(rs_cf,:),climsCF); 
hold on; plot([0 0],ylim,'k-');
hold on; plot(-[meanmvtime_all meanmvtime_all],ylim,'k-');
subaxis(1,2,2); imagesc(tmpxUse_s,[],tmpsig2{2}(rs_cf,:),climsCF);
h=gca; despos = h.Position; colorbar; h.Position=despos; 
hold on; plot([0 0],ylim,'k-');
hold on; plot(-[meanmvtime_all meanmvtime_all],ylim,'k-');
title(fp(f).name)

% show same anti GrC sorting with averages from random halves of trials
if 0
figure; 
subaxis(1,2,1); imagesc(tmpxUse_s,[],grcsigs_all{2}{3}(proportions_ix{1,1,1}(rs),tmpxUse),climsGrC); 
hold on; plot([0 0],ylim,'k-');
hold on; plot(-[meanmvtime_all meanmvtime_all],ylim,'k-');
subaxis(1,2,2); imagesc(tmpxUse_s,[],grcsigs_all{2}{4}(proportions_ix{1,1,1}(rs),tmpxUse),climsGrC);
h=gca; despos = h.Position; colorbar; h.Position=despos; 
hold on; plot([0 0],ylim,'k-');
hold on; plot(-[meanmvtime_all meanmvtime_all],ylim,'k-');
title(fp(f).name)
end

end
end
end
if fp(f).name=="3G-M", proportions_ix_3gm = proportions_ix; end
if pausebw, pause; end
end
mvprops = {premvprops{1},postmvprops{1},premvprops{2},postmvprops{2}};
barwitherr_stde(mvprops,dots=true,dotcalc="90")
titstr = "S2P Premove";
[~,s]=ranksum_disp(mvprops(1:2));
titstr = titstr+newline+s;
[~,s]=ranksum_disp(mvprops(3:4));
titstr = titstr+newline+s;
title(titstr);
xticklabels(["GrCs-pre","GrCs-post","CFs-pre","CFs-post"])
figureformat_forsaving