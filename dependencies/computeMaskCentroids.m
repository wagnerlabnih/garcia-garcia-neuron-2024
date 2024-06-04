function cents = computeMaskCentroids(newfilters)
nc = size(newfilters,3);
cents = zeros(nc,2);
% compute binary masks and centroids
thresh=2;
masksb = zeros(size(newfilters),'logical');
for k = 1:nc
    curc = newfilters(:,:,k);
    masksb(:,:,k) = curc>thresh*std(curc(:));
    stats = regionprops(squeeze(masksb(:,:,k)));
    if ~isempty(stats)
        cents(k,:) = stats.Centroid;
    end
end