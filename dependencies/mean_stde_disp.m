function s = mean_stde_disp(x,o)
arguments
    x
    o.dim = max(size(x));
    o.meds = false;
end

if ~iscell(x)
assert(length(size(x))<=2);
if nargin==1
    [~,o.dim] = max(size(x));
end
mu = squeeze(mean(x,o.dim));
tmpsd = squeeze(std_e(x,o.dim));
tmpnvals = repmat(size(x,o.dim),1,length(mu));
else
if nargin==1
    [~,o.dim] = max(size(x{1}));
end
if o.meds
    mu = cellfun(@(y)median(y,o.dim),x);
    tmpsd = cellfun(@(y)mad(y,o.dim),x);
else
mu = cellfun(@(y)mean(y,o.dim),x);
tmpsd = cellfun(@(y)std_e(y,o.dim),x);
end
tmpnvals = cellfun(@(y)size(y,o.dim),x);
end
% fprintf("n=%d\n",max(size(x)))
s = "";
for i = 1:length(mu)
    if i>1, s=s+newline; end
s = s+sprintf("m%d = %.3g+/-%.3g, n=%d",i,mu(i),tmpsd(i),tmpnvals(i));
end