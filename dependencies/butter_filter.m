function out = butter_filter(n,wn,x)
[b,a] = butter(n,wn);
out = filter(b,a,double(x));
end