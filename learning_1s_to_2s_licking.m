%%
clear; clc;
load("C:\Users\wagnermj\Desktop\data to share\learning_1s_to_2s_licking");
%% Initialize
ngl=2;
[rewlick, omitlick, omitrewalgn, rewdels, s_ix] = deal(cell(1,ngl));
for g = 1:ngl
    for d = 1:length(groupsl{g})
        curd= groupsl{g}{d};
        curd.rewdel = (curd.rewtimes-curd.midpt)*curd.dtb;
        if g>1, rewdeluse = 1.75; else, rewdeluse = 0.75; end
        rewuse = find(curd.rewarded & curd.goodmvdir & curd.rewdel>rewdeluse);
        [rewwin,rewwin_s] = winfn(curd.tmpxx,[-2 2]);
        omituse = find(~curd.rewarded & curd.goodmvdir & curd.mvlen>6 & curd.rewdel>rewdeluse);
        if isfield(curd,"laser"), 
            rewuse = rewuse(sum(curd.rewAlgn.laser(rewuse,rewwin),2)==0);
            omituse = omituse(sum(curd.rewAlgn.laser(omituse,rewwin),2)==0);
        end
        test = filter(ones(300,1)/300,1,double(curd.rewAlgn.lick)')';
        testwin = winfn(curd.tmpxx,[-1.75 1.75]);
        valid_lick = find(max(test(:,testwin),[],2)<0.99);
        invalid_lick = find(max(test(:,testwin),[],2)>=0.99);
        validtrials = valid_lick;
        rewuse_lick = intersect(rewuse,valid_lick);
        omituse_lick = intersect(omituse,valid_lick);

        filtwin = 60;
        lickraw = double(curd.midAlgn.lick);
        lickrate = (diff(lickraw,[],2)>0)./0.005;
        curlick = filter(ones(filtwin,1)/filtwin,1,lickrate')';
        rewlick{g} = cat(1,rewlick{g},curlick(rewuse_lick,:));
        omitlick{g} = cat(1,omitlick{g},curlick(omituse_lick,:));
        rewdels{g} = cat(1,rewdels{g},curd.rewdel([rewuse; omituse]));

        filtwin = 60;
        lickraw = double(curd.rewAlgn.lick);
        lickrate = (diff(lickraw,[],2)>0)./0.005;
        curlick = filter(ones(filtwin,1)/filtwin,1,lickrate')';
        omitrewalgn{g} = cat(1,omitrewalgn{g},curlick(omituse_lick,:));
    end
end
%% Fig. 2F, S2G, 2G,H, S3H,I, S3G
close all
xpltwin = [-1 4];
[xplt_cur_b,xplt_cur_b_s] = winfn(curd.tmpxx,xpltwin);
fo=figure(pos=[1158         299         345         420]); ho=[];
fr=figure(pos=[1158         299         345         420]); hr=[];
guse = [1 2]; ngu=length(guse); ccplt = [0 0 0; 0 .5 .5]; ccplt = ccplt(guse,:);
titstro = "";
titstrr = "";
delay_lick_duration = cell(ngu,1);
for g = 1:ngu
gg=guse(g);
cur = omitlick{gg};
figure(fo);
ho(end+1)=errorbar_shadeSEM(xplt_cur_b_s,cur(:,xplt_cur_b)',ccplt(g,:));
titstro=titstro+sprintf(" n-"+g+"s:%d, ",size(omitlick{gg},1));
cur = rewlick{gg};
figure(fr);
hr(end+1)=errorbar_shadeSEM(xplt_cur_b_s,cur(:,xplt_cur_b)',ccplt(g,:));
titstrr=titstrr+sprintf(" n-"+g+"s:%d, ",size(rewlick{gg},1));
hold on;
end
titstrs = [titstrr,titstro];
fs = [fr fo]; hs={hr,ho}; yls=[8 7];
for i = 1:2
    figure(fs(i));
    uistack(hs{i},"top"); 
    ylim([0 yls(i)])
    plot([0 0],ylim,'k');
    plot([1 1],ylim,'k');
    plot(2*[1 1],ylim,'k');
    xlabel("Time rel to movement (s)")
    ylabel("Lick rate (Hz)")
    title(titstrs(i))
end

pkTime = cell(1,ngu);
pkTimeRA = cell(1,ngu);
winomit2 = winfn(curd.tmpxx,[2.1 2.6]);
lickratepost2 = cellfun(@(x)mean(x(:,winomit2),2),omitlick(guse),unif=0);
winbase = winfn(curd.tmpxx,[0.3 0.8]);
lickratebase = cellfun(@(x)mean(x(:,winbase),2),omitlick(guse),unif=0);
lickratepost2mbase = cellfun(@(x)(mean(x(:,winomit2),2)-mean(x(:,winbase),2)),...
    omitlick(guse),unif=0);
for g = 1:ngu
    gg=guse(g);
    cur = omitlick{gg};
    [xtest,xtest_s] = winfn(curd.tmpxx,[0 3]);
    smoothed_lick = filtfilt(ones(40,1)/40,1,cur')';
    smoothed_lick = smoothed_lick(:,xtest);
    for t = 1:size(smoothed_lick,1)
        pkTime{g} = cat(1,pkTime{g},xtest_s(maxix(smoothed_lick(t,:))));
    end
    
    cur = omitrewalgn{gg};
    [xtest,xtest_s] = winfn(curd.tmpxx,[-4 4]);
    smoothed_lick = filtfilt(ones(40,1)/40,1,cur')';
    smoothed_lick = smoothed_lick(:,xtest);
    for t = 1:size(smoothed_lick,1)
        pkTimeRA{g} = cat(1,pkTimeRA{g},xtest_s(maxix(smoothed_lick(t,:))));
    end

    ix_2sdel_b = winfn(curd.tmpxx,[0 2]);
    premv_s = [-1 0];
    dtb = 0.005;
    ix_premv_b = winfn(curd.tmpxx,premv_s);
    delay_lick_duration{g} = sum(rewlick{g}(:,ix_2sdel_b)-mean(rewlick{g}(:,ix_premv_b),2)>0,2)*dtb;
end

histogram_cell(lickratepost2mbase,cdf=true,pos=[898   291   262   262],col=ccplt)
hold on; plot([0 0],ylim,'k')

histogram_cell(pkTime,cdf=true,pos=[898   291   262   262],col=ccplt,medians=true)
hold on; plot([1 1],ylim,'k')
hold on; plot([2 2],ylim,'k')

histogram_cell(delay_lick_duration,cdf=true,col=ccplt);
xticks(0:.5:2); 
hold on; plot([1 1],ylim,'k--'); plot([2 2],ylim,'k')
ylabel("Fraction of trials")
xlabel("Duration of licking during [0,2] (s)"); box off;
barwitherr_stde(delay_lick_duration,color=ccplt);
ylim tight;
ylabel("Duration of licking during [0,2] (s)")

figure(pos=[922   401   171   265]);
[xdel_b,xdel_b_s]=winfn(curd.tmpxx,[-0.5 2.5]);
errorbar_shadeSEM(xdel_b_s,rewlick{1}(:,xdel_b)',ccplt(1,:));
hold on;
errorbar_shadeSEM(xdel_b_s,rewlick{2}(:,xdel_b)',ccplt(2,:))
ylim([0 8]);
hold on; plot([1 1],ylim,'k--'); plot([2 2],ylim,'k--'); plot([0 0],ylim,'k')