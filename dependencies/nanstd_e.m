function se = nanstd_e(x,d)

if nargin==1, d = 1; end

if isrow(x), x = x'; end

se = nanstd(x,[],d)./sqrt(size(x,d));