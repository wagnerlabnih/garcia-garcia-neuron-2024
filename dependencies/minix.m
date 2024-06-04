function ix = minix(x,y,d)
if isrow(x), x = x'; end
if nargin == 1, y = []; end
if nargin < 3, d = 1; end
[~,ix] = min(x,y,d);