function [p,strout] = ranksum_disp(x,y,o)
arguments
    x
    y=[];
    o.ks = false;
    o.meds = false;
end

if isempty(y),
    if min(size(x))==1
    y = x{2};
    x = x{1};
    else
        y = x(:,2);
        x = x(:,1);
    end
end
if isrow(x), x = x'; end; if isrow(y), y = y'; end
if ~iscell(x), x = {x}; y = {y}; end
for i = 1:length(x)
    xc = x{i}; yc = y{i};
if all(ismember(unique([xc(:); yc(:)]),[0 1])),
    xtab = cellfun(@(x)[nnz(x) length(find(x==0))],{xc,yc},unif=0);
    xtab = cat(1,xtab{:});
    [~,p,~] = fishertest(xtab);%,'Tail','left');
else 
    if ~o.ks
    p = ranksum(xc,yc);
    else
    [~,p] = kstest2(xc,yc);
    end
end
if o.meds, fn1 = @median; fn2 = @mad; else, fn1=@mean; fn2=@std_e; end
strout = sprintf("m1 = %1.2g+/-%1.2g\nm2 = %1.2g+/-%1.2g\np = %1.2g, n1=%d, n2=%d\n",fn1(xc),fn2(xc),fn1(yc),fn2(yc),p,length(xc),length(yc));
if nargout<2
disp(strout);
end
end