function [p,chi2,s] = chisq_disp(x,o)
arguments 
    x
    o.y = [];
    o.disp = true;
end
if isempty(o.y) && iscell(x)
    o.y = arrayfun(@(z)z*ones(length(x{z}),1),1:length(x),unif=false);
    o.y = cat(1,o.y{:});
    dimuse = find(size(x{1})>1);
    x = cat(dimuse,x{:});
end
[~,chi2,p] = crosstab(x,o.y);
s = sprintf("p = %.3g, chi2=%2.1g, n=%d\n",p,chi2,length(o.y));
if o.disp
disp(s);
end