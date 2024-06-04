function [preds,mae,r2] = kfoldpred_fitlm(X,y)

ns = length(y);
% ixp = randperm(ns);
ixp = 1:ns;
npf = floor(ns/10);
preds = zeros(size(y));
for i = 1:10
    if i<10
        curvalix = ((i-1)*npf+1):(i*npf);
    else
        curvalix = ((i-1)*npf+1):ns;
    end
    val_fold = ixp(curvalix);
    tr_fold = setdiff(ixp,val_fold);
    b = fitlm(X(tr_fold,:),y(tr_fold));
    b = b.Coefficients.Estimate;
    preds(val_fold) = X(val_fold,:)*b(2:end)+b(1);
end
% pause
mae = mean(abs(preds-y));
r2 = corr(preds,y).^2;