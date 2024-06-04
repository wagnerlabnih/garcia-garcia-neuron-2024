function out = shufflev(x,ns,n_subs)
if nargin<2, ns=1; end
x = x(:);
out = x(randperms(ns,length(x)));
if nargin==3 && ~isempty(n_subs)
    out = out(1:n_subs,:);
end