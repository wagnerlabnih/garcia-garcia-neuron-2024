function filtsout = clearMaskBackground(newfilters,o)

arguments
    newfilters
    o.smoothedfilts = imgaussfilt(newfilters,1)
    o.thresh = 2
    o.cellszthresh = 50 %200
    o.dispresults=true
end

mins = min(min(newfilters,[],1),[],2);
newfilters = bsxfun(@plus,newfilters,-mins);
newfilters = bsxfun(@times,newfilters,1./sum(sum(newfilters,1),2));
nc = size(newfilters,3);
stds = std(reshape(o.smoothedfilts,[],nc),[],1);
means = mean(reshape(o.smoothedfilts,[],nc),1);
if o.dispresults
f=figure('position',[617    477	1289	220]);
f.WindowStyle='docked';
end
for k = 1:nc
    curc = newfilters(:,:,k); 
    curcsmoothed = o.smoothedfilts(:,:,k);
    curc(curcsmoothed<means(k)+o.thresh*stds(k)) = 0;
    curcsmoothed(curcsmoothed<means(k)+o.thresh*stds(k)) = 0;
    % delete crud
    curcbw = curcsmoothed>0;
    regs = bwconncomp(curcbw);
    rp = regionprops(regs,'Area'); rpa = cat(1,rp.Area);
    regixdel = find(rpa<o.cellszthresh);
    rppix = cat(1,regs.PixelIdxList);
    rppixdel = cat(1,rppix{regixdel});
    curc(rppixdel)=0;
    if o.dispresults
    subaxis(1,2,1); imagesc(newfilters(:,:,k)); axis image;
    end
    newfilters(:,:,k) = curc;
    if o.dispresults
    subaxis(1,2,2); imagesc(newfilters(:,:,k)); axis image;
    pause; clf
    end
end
% delete empty cells
remempt = [];
for k = 1:nc
    if nnz(newfilters(:,:,k))==0, remempt = cat(1,remempt,k); end
end
newfilters(:,:,remempt) = []; 
filtsout = newfilters;
if o.dispresults, close(f); end
end