function [ step ] = getStep( min_el, max_el )
%GETSTEP Summary of this function goes here
%   Detailed explanation goes here

if ((max_el - min_el) <= 100)
    step = 1;
else 
    step  = ceil((max_el - min_el)/50);
end;

return;


