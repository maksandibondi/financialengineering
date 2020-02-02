function [out] = nonUniTransform(y, ymax, y0)


out = 0.5*log((ymax+y)/(ymax-y))+y0;