function   LLPAMSProcessNPR_VecComputVersion
global msr
global MAX_OUTPUT_BUFFERS
global mode
global FLT_EPSILON            
global DBL_EPSILON
  % find dimensions
        freqs = msr.halfLen + 2; % two-sided + nyquist point
%       TempLBuffer = zeros(1,msr.frameLen);
%       TempRBuffer = zeros(1,msr.frameLen);
%           for i = 0: msr.frameLen-1
%               k = i + msr.mInputPos;
%               if (k >= msr.frameLen)
%                     k = k - msr.frameLen;
%               end
%               w = msr.analysisWnd(i+1);
%               TempLBuffer(i+1) = msr.mInput(1,k+1) * w;
%               TempRBuffer(i+1) = msr.mInput(2,k+1) * w;
%           end
          
           idx_i = 0: msr.frameLen-1;
           idx_k = idx_i + msr.mInputPos;
            idx_k_sub = idx_k >=msr.frameLen;
            idx_k(idx_k_sub) = idx_k(idx_k_sub) - msr.frameLen;
            TempLBuffer = msr.mInput(1,idx_k+1) .* msr.analysisWnd.';
            TempRBuffer = msr.mInput(2,idx_k+1) .* msr.analysisWnd.';

           xr_fft = fft(TempRBuffer, msr.fftLen);
           xr_fft = xr_fft(1:freqs);

           xl_fft = fft(TempLBuffer, msr.fftLen);
           xl_fft = xl_fft(1:freqs);

            lI = imag(xl_fft);
            lR = real(xl_fft);
            rI = imag(xr_fft);
            rR = real(xr_fft);
            magnitudeLeft = sqrt(lI.^2 + lR.^2);
            magnitudeRight = sqrt(rI.^2 + rR.^2);
            minMag = min(magnitudeLeft, magnitudeRight);
           
            ratMinLeft = minMag./(magnitudeLeft + FLT_EPSILON);
            idx0 = abs(ratMinLeft - msr.smoothedFunc(9,:)) < msr.smoothing;
            msr.smoothedFunc(9,idx0) = ratMinLeft(idx0); 
            idx1 = ratMinLeft > msr.smoothedFunc(9,:);
            msr.smoothedFunc(9,idx1) = msr.smoothedFunc(9,idx1) + msr.smoothing;
            msr.smoothedFunc(9,~(idx0|idx1)) = msr.smoothedFunc(9,~(idx0|idx1)) - msr.smoothing;
           
            ratMinRight = minMag/(magnitudeRight + FLT_EPSILON);
            idx0 = abs(ratMinRight - msr.smoothedFunc(8,:)) < msr.smoothing;
            msr.smoothedFunc(8,idx0) = ratMinRight; 
            idx1 = ratMinRight > msr.smoothedFunc(8,:);
            msr.smoothedFunc(8,idx1) = msr.smoothedFunc(8,idx1) + msr.smoothing;            
            msr.smoothedFunc(8,~(idx0|idx1)) = msr.smoothedFunc(8,~(idx0|idx1)) - msr.smoothing;
            
            magnitudeSum = magnitudeLeft + magnitudeRight; 
            idx = magnitudeSum < DBL_EPSILON;
            magDiff(idx) = 0;
            magDiff(~idx) = (magnitudeRight(~idx) - magnitudeLeft(~idx))./ magnitudeSum(~idx);
            magDiff = max(-1.0, min(1.0, magDiff));
            phaseL = atan2(lI,lR);
            phaseR = atan2(rI,rR);
            phaseDiff = abs(phaseR - phaseL);
            idx  = phaseDiff > pi; 
            phaseDiff(idx) = 2*pi - phaseDiff(idx);
             
%             x = zeros(1,freqs);
%             y = x;
%             for idx = 1:freqs
%              [x(idx),y(idx)] = cartesianMap(msr.magPhaseDiff2Cartesian,magDiff(idx),phaseDiff(idx));
%             end
            [x,y] = cartesianMap(msr.magPhaseDiff2Cartesian,magDiff,phaseDiff);


             leftSineTerm = sin(phaseL);
             leftCosineTerm = cos(phaseL);

             rightSineTerm = sin(phaseR);
             rightCosineTerm = cos(phaseR);

             centrePhase = atan2((lI + rI),(lR + rR));
             centreCosineTerm = cos(centrePhase); 
             centreSineTerm = sin(centrePhase); 

%              x0 = zeros(1,freqs);
%              y0 = x0;
%              u = x0;
%              v = x0;
%              for idx = 1:freqs
%                  [u(idx),x0(idx)] = map_to_grid(x(idx), 21); 
%                  [v(idx),y0(idx)] = map_to_grid(y(idx), 21);
%              end
             [u,x0] = map_to_grid(x, 21); 
             [v,y0] = map_to_grid(y, 21);
             
            magSqrt = sqrt(magnitudeLeft.^2 + magnitudeRight.^2);

            opLFwd = xl_fft - xr_fft.* msr.smoothedFunc(8,:); 
            opRFwd = xr_fft - xl_fft.* msr.smoothedFunc(9,:); 

            pmix = msr.mix;
            minusPmix = 1.0 - pmix;
            
            magnitude = lerpCompCplx_VecComputVersion(1,magSqrt,msr.chmaskPtr,u,v,x0,y0);
            real_part = magnitude .* leftCosineTerm;
            imag_part = magnitude .* leftSineTerm;
            msr.timeDomainOut(1,1:freqs) = real_part + 1i*imag_part;
% 
           magnitude = lerpCompCplx_VecComputVersion(2,magSqrt,msr.chmaskPtr,u,v,x0,y0);
           real_part = magnitude .* rightCosineTerm;
           imag_part = magnitude .* rightSineTerm;
           msr.timeDomainOut(2,1:freqs) = real_part + 1i*imag_part;

            magnitude = lerpCompCplx_VecComputVersion(3,magSqrt,msr.chmaskPtr,u,v,x0,y0);
            real_part = magnitude .* centreCosineTerm;
            imag_part = magnitude .* centreSineTerm;
            msr.timeDomainOut(3,1:freqs) = real_part + 1i*imag_part;

            magnitude = lerpCompCplx_VecComputVersion(4,magSqrt,msr.chmaskPtr,u,v,x0,y0);
            real_part = magnitude .* leftCosineTerm;
            imag_part = magnitude .* leftSineTerm;
            msr.timeDomainOut(4,1:freqs) = (real_part + 1i*imag_part)*minusPmix + opLFwd * pmix;

            magnitude = lerpCompCplx_VecComputVersion(5,magSqrt,msr.chmaskPtr,u,v,x0,y0);
            real_part = magnitude .* rightCosineTerm;
            imag_part = magnitude .* rightSineTerm;
            msr.timeDomainOut(5,1:freqs) = (real_part + 1i*imag_part)*minusPmix + opRFwd * pmix;   

           for j = 1:mode
              tmp = [msr.timeDomainOut(j,1:freqs) conj(msr.timeDomainOut(j,1:freqs))];
              msr.timeDomainOut(j,1:msr.fftLen) = real(ifft(tmp,msr.fftLen)).';
           end
           msr.mOutputBufferCount = msr.mOutputBufferCount + 1;

         if msr.mOutputBufferCount > MAX_OUTPUT_BUFFERS
              return;
         end

          idx_kk = 0:msr.ovpLen-1;
            idx_mm = idx_kk + msr.mInputPos + msr.smpShift;
            idx_mm_sub = idx_mm >= msr.frameLen;
            idx_mm(idx_mm_sub) = idx_mm(idx_mm_sub) - msr.frameLen;
   
            msr.mOutputBuffer(1,idx_kk+1) = msr.mOverlapStage2dash(1,idx_kk+1) + msr.timeDomainOut(1,1+idx_kk  + msr.smpShift) .* msr.synthesisWnd(idx_kk+1).' + msr.mInputSub(1,idx_mm+1) * msr.minuscrossBass;%20%的低频混到左声道去
            msr.mOverlapStage2dash(1,idx_kk +1) = msr.timeDomainOut(1,1+msr.smpShift + msr.ovpLen + idx_kk ).* msr.synthesisWnd(1+idx_kk  + msr.ovpLen).';
            msr.mOutputBuffer(2,idx_kk +1) = msr.mOverlapStage2dash(2,idx_kk +1) + msr.timeDomainOut(2,1+idx_kk  + msr.smpShift).* msr.synthesisWnd(idx_kk +1).' + msr.mInputSub(1,idx_mm+1) * msr.minuscrossBass;%20%的低频混到右声道去
            msr.mOverlapStage2dash(2,idx_kk +1) = msr.timeDomainOut(2,1+msr.smpShift + msr.ovpLen + idx_kk ) .* msr.synthesisWnd(1+idx_kk  + msr.ovpLen).';

            for ll = 3:5
                msr.mOutputBuffer(ll,idx_kk +1) = msr.mOverlapStage2dash(ll,idx_kk +1) + msr.timeDomainOut(ll,1+idx_kk  + msr.smpShift) .* msr.synthesisWnd(idx_kk +1).';
                msr.mOverlapStage2dash(ll,idx_kk +1) = msr.timeDomainOut(ll,1+msr.smpShift + msr.ovpLen + idx_kk ).* msr.synthesisWnd(1+idx_kk  + msr.ovpLen).';
            end
             msr.mOutputBuffer(6,idx_kk +1) = (msr.mInputSub(1,idx_mm+1) + msr.mInputSub(2,idx_mm+1)) * msr.crossBass;
           msr.mInputSamplesNeeded = msr.ovpLen;
end
