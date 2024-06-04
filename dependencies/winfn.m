function [out,winout] = winfn(x,win)
% x = time axis variable (must be a monotonically increasing vector of numbers)
% win must be a 2-element vector containing a desired start and end time
% out is a vector of indeces into x such that out(1) is the element of x
% closest to win(1) and out(end) is the element of x closest to win(2)
if numel(x)==2 && numel(win)>2, tmpx = x; x = win; win = tmpx; clear tmpx; end
out = minix(abs(x-win(1))):minix(abs(x-win(2)));
if nargin==2, winout = x(out); end
end