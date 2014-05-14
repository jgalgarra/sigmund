%function rperiod( r, period )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    %(1-((1-abs(r))^(1/period)));
    %r
%end

function y = rperiod(x)

 y =x.^2 + (1 - exp(x)).^2 - 4;
end