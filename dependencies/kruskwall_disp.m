function [p,s] = kruskwall_disp(x,o)
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
p = kruskalwallis(x,o.y,'off');
s = sprintf("p = %.3g, n=%d\n",p,length(o.y));
if o.disp
disp(s);
end