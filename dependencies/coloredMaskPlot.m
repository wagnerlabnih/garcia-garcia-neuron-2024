function [h,ccmult]=coloredMaskPlot(newfilters,o)
arguments
    newfilters
    o.meanim = zeros(size(newfilters(:,:,1)))
    o.colors = rand(size(newfilters,3),3)*0.7 + 0.15
    o.hicut = 1
    o.curfig = false
    o.textlabls = false;
end
newfilters = double(newfilters);
o.hicut = max(max(o.meanim(:)),1);
nc = size(newfilters,3);
flatfilts = reshape(newfilters,[],nc);
tmpsd = 1./(sum(flatfilts,1)./sum(flatfilts~=0,1));
filtersplot2 = permute(bsxfun(@times,permute(newfilters,[3 1 2]),tmpsd'),[2 3 1]);
filtersplotc = permute(filtersplot2,[1 2 4 3]);
filtersplotc = repmat(filtersplotc,[1 1 3 1]);
ccmult_p = permute(o.colors,[3 4 2 1]);
filtersplotc = bsxfun(@times,filtersplotc,ccmult_p);
cents = computeMaskCentroids(newfilters);
clear filtersplot;
filtsumc = squeeze(sum(filtersplotc,4));
meanimrep = repmat(o.meanim,[1 1 3]); meanimrep = meanimrep./max(meanimrep(:));
meanimrep(meanimrep>o.hicut) = o.hicut; meanimrep(meanimrep<0) = 0;
meanimrep = meanimrep./max(meanimrep(:));
meanimrep(filtsumc>0) = 0;
if ~o.curfig, figure; end
imagesc(gca,filtsumc+meanimrep); axis image;
% imagesc(gca,filtsumc); axis image; shg;
% imagesc(gca,meanimrep); axis image; shg;
text(cents(:,1),cents(:,2),cellstr(num2str((1:nc)')));
h = gcf;
ccmult = o.colors;