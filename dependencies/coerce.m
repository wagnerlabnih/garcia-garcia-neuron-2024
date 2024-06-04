function out=coerce(in,outrange)
in(in<outrange(1)) = outrange(1);
in(in>outrange(2)) = outrange(2);
out = in;
end