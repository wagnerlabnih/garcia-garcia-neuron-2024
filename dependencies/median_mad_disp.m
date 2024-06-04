function s = median_mad_disp(x,dim)

assert(length(size(x))<=2);
if nargin==1
    [~,dim] = max(size(x));
end
mu = squeeze(median(x,dim));
tmpsd = squeeze(mad(x,1,dim));
s="";
s = s+newline+sprintf("n=%d\n",max(size(x)));
for i = 1:length(mu)
s= s+newline+sprintf("mu%d = %.3g+/-%.3g,n=%d",i,mu(i),tmpsd(i),size(x,dim));
end