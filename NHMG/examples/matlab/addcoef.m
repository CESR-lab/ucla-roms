function [iT,jT,cT,v] = addcoef(iT,jT,cT,v,coef,idx,jdx); 
 
   nd = length(coef(:));
   iT(v:v+nd-1) = idx;
   jT(v:v+nd-1) = jdx;
   cT(v:v+nd-1) = coef(:);
   v = v + nd;
  
