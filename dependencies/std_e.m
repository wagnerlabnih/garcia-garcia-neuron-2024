function se = std_e(x,d,isr)

if nargin<2 || isempty(d), d = 1; end
if nargin<3, isr = false; end

if isr && isrow(x), x = x'; end

se = std(x,[],d)./sqrt(size(x,d));