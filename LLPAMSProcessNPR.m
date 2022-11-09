function   LLPAMSProcessNPR
global msr
global MAX_OUTPUT_BUFFERS
global mode
global FLT_EPSILON            
global DBL_EPSILON
  % find dimensions
        freqs = msr.halfLen + 2; % two-sided + nyquist point
        TempLBuffer = zeros(1,msr.frameLen);
        TempRBuffer = zeros(1,msr.frameLen);
          for i = 0: msr.frameLen-1
              k = i + msr.mInputPos;
              if (k >= msr.frameLen)
                    k = k - msr.frameLen;
              end
              w = msr.analysisWnd(i+1);
              TempLBuffer(i+1) = msr.mInput(1,k+1) * w;
              TempRBuffer(i+1) = msr.mInput(2,k+1) * w;
          end

           xr_fft = fft(TempRBuffer, msr.fftLen);
           xr_fft = xr_fft(1:freqs);

           xl_fft = fft(TempLBuffer, msr.fftLen);
           xl_fft = xl_fft(1:freqs);

        for i = 1:freqs
            lI = imag(xl_fft(i));
            lR = real(xl_fft(i));
            rI = imag(xr_fft(i));
            rR = real(xr_fft(i));
            magnitudeLeft = sqrt(lI^2 + lR^2);
            magnitudeRight = sqrt(rI^2 + rR^2);
            minMag = min(magnitudeLeft, magnitudeRight);
            ratMinLeft = minMag/(magnitudeLeft + FLT_EPSILON);
             if (abs(ratMinLeft - msr.smoothedFunc(9,i)) < msr.smoothing)
                 msr.smoothedFunc(9,i) = ratMinLeft; 
             elseif (ratMinLeft > msr.smoothedFunc(9,i))
                 msr.smoothedFunc(9,i) = msr.smoothedFunc(9,i) + msr.smoothing;
             else
                 msr.smoothedFunc(9,i) = msr.smoothedFunc(9,i) - msr.smoothing;
             end
             ratMinRight = minMag/(magnitudeRight + FLT_EPSILON);
             if (abs(ratMinRight - msr.smoothedFunc(8,i)) < msr.smoothing)
                 msr.smoothedFunc(8,i) = ratMinRight; 
             elseif (ratMinRight > msr.smoothedFunc(8,i))
                 msr.smoothedFunc(8,i) = msr.smoothedFunc(8,i) + msr.smoothing;
             else
                 msr.smoothedFunc(8,i) = msr.smoothedFunc(8,i) - msr.smoothing;
             end
             magnitudeSum = magnitudeLeft + magnitudeRight; 
             if(magnitudeSum < DBL_EPSILON)
                 magDiff = 0;
             else
                 magDiff = (magnitudeRight - magnitudeLeft)/ magnitudeSum;
             end
             magDiff = max(-1.0, min(1.0, magDiff));
             phaseL = atan2(lI,lR);
             phaseR = atan2(rI,rR);
             phaseDiff = abs(phaseR - phaseL);
             if (phaseDiff > pi) 
                    phaseDiff = 2*pi - phaseDiff;
             end
             [x,y] = cartesianMap(msr.magPhaseDiff2Cartesian,magDiff,phaseDiff);

             leftSineTerm = sin(phaseL);
             leftCosineTerm = cos(phaseL);

             rightSineTerm = sin(phaseR);
             rightCosineTerm = cos(phaseR);

             centrePhase = atan2((lI + rI),(lR + rR));
             centreCosineTerm = cos(centrePhase); 
             centreSineTerm = sin(centrePhase); 

             [u,x] = map_to_grid(x, 21); 
             [v,y] = map_to_grid(y, 21);

             magSqrt = sqrt(magnitudeLeft^2 + magnitudeRight^2);

            opLFwd = xl_fft(i) - xr_fft(i) * msr.smoothedFunc(8,i); 
            opRFwd = xr_fft(i) - xl_fft(i)  * msr.smoothedFunc(9,i); 

            pmix = msr.mix;
            minusPmix = 1.0 - pmix;

            magnitude = lerpCompCplx(1,i,magSqrt,msr.chmaskPtr,u,v,x,y);
            real_part = magnitude * leftCosineTerm;
            imag_part = magnitude * leftSineTerm;
            msr.timeDomainOut(1,i) = real_part + 1i*imag_part;

           magnitude = lerpCompCplx(2,i,magSqrt,msr.chmaskPtr,u,v,x,y);
           real_part = magnitude * rightCosineTerm;
           imag_part = magnitude * rightSineTerm;
           msr.timeDomainOut(2,i) = real_part + 1i*imag_part;

            magnitude = lerpCompCplx(3,i,magSqrt,msr.chmaskPtr,u,v,x,y);
            real_part = magnitude * centreCosineTerm;
            imag_part = magnitude * centreSineTerm;
            msr.timeDomainOut(3,i) = real_part + 1i*imag_part;

            magnitude = lerpCompCplx(4,i,magSqrt,msr.chmaskPtr,u,v,x,y);
            real_part = magnitude * leftCosineTerm;
            imag_part = magnitude * leftSineTerm;
            msr.timeDomainOut(4,i) = (real_part + 1i*imag_part)*minusPmix + opLFwd * pmix;

            magnitude = lerpCompCplx(5,i,magSqrt,msr.chmaskPtr,u,v,x,y);
            real_part = magnitude * rightCosineTerm;
            imag_part = magnitude * rightSineTerm;
            msr.timeDomainOut(5,i) = (real_part + 1i*imag_part)*minusPmix + opRFwd * pmix;   

        end

         for j = 1:mode
             tmp = [msr.timeDomainOut(j,1:freqs) conj(msr.timeDomainOut(j,1:freqs))];
             msr.timeDomainOut(j,1:msr.fftLen) = real(ifft(tmp,msr.fftLen)).';
         end
            msr.mOutputBufferCount = msr.mOutputBufferCount + 1;

         if msr.mOutputBufferCount > MAX_OUTPUT_BUFFERS
              return;
         end

        for kk = 0:msr.ovpLen-1  
            mm = kk + msr.mInputPos + msr.smpShift;
            if(mm >= msr.frameLen)
                mm = mm - msr.frameLen;
            end

            msr.mOutputBuffer(1,kk+1) = msr.mOverlapStage2dash(1,kk+1) + (msr.timeDomainOut(1,1+kk + msr.smpShift) * msr.synthesisWnd(kk+1)) + msr.mInputSub(1,mm+1) * msr.minuscrossBass;%20%的低频混到左声道去
            msr.mOverlapStage2dash(1,kk+1) = msr.timeDomainOut(1,1+msr.smpShift + msr.ovpLen + kk) * msr.synthesisWnd(1+kk + msr.ovpLen);
            msr.mOutputBuffer(2,kk+1) = msr.mOverlapStage2dash(2,kk+1) + (msr.timeDomainOut(2,1+kk + msr.smpShift) * msr.synthesisWnd(kk+1)) + msr.mInputSub(1,mm+1) * msr.minuscrossBass;%20%的低频混到右声道去
            msr.mOverlapStage2dash(2,kk+1) = msr.timeDomainOut(2,1+msr.smpShift + msr.ovpLen + kk) * msr.synthesisWnd(1+kk + msr.ovpLen);

           for ll = 3:5
                msr.mOutputBuffer(ll,kk+1) = msr.mOverlapStage2dash(ll,kk+1) + (msr.timeDomainOut(ll,1+kk + msr.smpShift) * msr.synthesisWnd(kk+1));
                msr.mOverlapStage2dash(ll,kk+1) = msr.timeDomainOut(ll,1+msr.smpShift + msr.ovpLen + kk) * msr.synthesisWnd(1+kk + msr.ovpLen);
           end
            msr.mOutputBuffer(6,kk+1) = (msr.mInputSub(1,mm+1) + msr.mInputSub(2,mm+1)) * msr.crossBass;
        end
            msr.mInputSamplesNeeded = msr.ovpLen;
end
