function out = butter_filtfilt(n,wn,x)
[b,a] = butter(n,wn);
out = filtfilt(b,a,double(x));
end