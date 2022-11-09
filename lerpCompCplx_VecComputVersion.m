function magnitude = lerpCompCplx_VecComputVersion(idx,magSqrt,chmaskPtr,u,v,x,y)
    global msr
   
    a = chmaskPtr((idx-1)*21*21+1:idx*21*21); 
  
    gf = ((1.0 - x).*(1.0 - y).*a(u * 21 + v) + x .* (1.0 - y).*a((u + 1) * 21 + v) + (1.0 - x).*y.*a(u * 21 + v + 1) + x .* y.*a((u + 1) * 21 + v + 1));
    difGf = gf - msr.smoothedFunc(idx,:); 
   
    idx0 = difGf > msr.smoothing;
    gf(idx0) = msr.smoothedFunc(idx,idx0) + msr.smoothing; 
%     
    idx0 = difGf < -msr.smoothing;
    gf(idx0) = msr.smoothedFunc(idx,idx0) - msr.smoothing; 
   
    msr.smoothedFunc(idx,:) = gf; 
    magnitude = magSqrt .* gf;

end