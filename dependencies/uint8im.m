function out = uint8im(im,range)

range = double(range);
im2 = max(min(double(im),range(2)),range(1));
im2 = uint8(round((im2-range(1))./(range(2)-range(1))*255));
out = im2;
