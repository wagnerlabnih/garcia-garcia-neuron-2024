function [f,p,s] = anova_disp(x,o)
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
tbl = anova(fitlm(x,o.y));
f = tbl.F(1);
p = tbl.pValue(1);
s = sprintf("f = %.3g, p = %.3g, n=%d\n",f,p,length(o.y));
s = s + sprintf("n=%d, ",arrayfun(@(xx)nnz(o.y==xx),1:length(unique(o.y))));
if o.disp
disp(s);
end