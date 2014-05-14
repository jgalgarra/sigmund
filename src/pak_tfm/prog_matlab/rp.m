function [ y ] = rp( ran,per )
    %a = double(log(1+ran));
    %b = double(exp(a/per));
    %y = double(b-1);
    y = (exp(ran/per)-1)
end

