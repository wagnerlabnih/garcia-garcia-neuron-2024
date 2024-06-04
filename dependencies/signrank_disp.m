function [s,p] = signrank_disp(x,y,flag)

if nargin<2 | isempty(y),
    y = [];
end
if nargin<3
    flag = 0;
end
if isrow(x), x = x'; end
if isrow(y), y = y'; end
p = zeros(size(x,2),1);
for c = 1:size(x,2)
if isempty(y)
p(c) = signrank(x(:,c));
else
p(c) = signrank(x(:,c),y(:,c));
end
s = sprintf("p = %.3g, n=%d",p,length(x));
if flag==0, 
    s = s+newline+mean_stde_disp(x);
    if ~isempty(y)
        s = s+newline+mean_stde_disp(y);
        s = s+newline+"difference: "+mean_stde_disp(x-y);
    end
else
    median_mad_disp(x);
    if ~isempty(y)
        s = s+newline+median_mad_disp(y);
    end
end
end