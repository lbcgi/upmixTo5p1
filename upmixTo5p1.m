clear;

% load angleMapper.mat
load chmaskPtr5p1.mat
load AmpDiffPhaseDiff2xy.mat

% global angleMapper
global ANALYSIS_OVERLAP
global msr
global lr2
global TBL_SIZE
global FFTSIZE
global MAX_OUTPUT_BUFFERS
global mode
global FLT_EPSILON            
global DBL_EPSILON      

DBL_EPSILON = 2.2204460492503131e-016;
FLT_EPSILON = 1.192092896e-07;

% filename = 'sines_with_delays_lr.wav';%文件名
filename = '彭芳 - 老人与海.wav';
[x, fs] = audioread(strcat(filename));

mode = 5; % 5.1

MAX_OUTPUT_BUFFERS = 2;
TBL_SIZE = 64;
FFTSIZE = 8192;
ANALYSIS_OVERLAP = 4;

msr.chmaskPtr = chmaskPtr;
msr.mBitRev = zeros(1,FFTSIZE);
msr.mOutputReadSampleOffset = 0;
msr.mOutputBufferCount = 0;
msr.mInputSub = zeros(2,FFTSIZE);
msr.mOutputBufferCount = 0;
msr.mInputPos = 0;
msr.magPhaseDiff2Cartesian = AmpDiffPhaseDiff2xy;
msr.spread = pi / 2.0;
msr.separation = 0.1;
msr.spatialEnhancement = 1.0;
msr.XYaxis = linspace(-1,1,TBL_SIZE);

msr.fs = fs;
msr.fc = 70.0;%Subwoofer cross frequency point

msr.mInput = zeros(2,FFTSIZE);
msr.frameLen = 1536;%1536
msr.fftLen = 2048;
msr.timeDomainOut = zeros(7,FFTSIZE);
msr.ovpLen = msr.frameLen / ANALYSIS_OVERLAP;
msr.smpShift = msr.frameLen - 2*msr.ovpLen;
msr.mOutputBuffer = zeros(7,msr.ovpLen);
msr.mOverlapStage2dash = zeros(7,FFTSIZE / ANALYSIS_OVERLAP);
msr.halfLen = msr.fftLen/2 - 1;
msr.procUpTo = 24000/(msr.fs/msr.fftLen);
msr.mInputSamplesNeeded = msr.ovpLen;

[msr.analysisWnd,msr.synthesisWnd] = getAsymmetricWindow(msr.frameLen,msr.ovpLen,msr.smpShift);

% LUTSurroundRefreshParameter;

sm = 15;
msr.smoothing = sm / msr.fs * (msr.frameLen / ANALYSIS_OVERLAP);
msr.mix = 1.0 - tanh(0.3 * sqrt(sm));
msr.smoothedFunc = 0.5*ones(9,msr.fftLen/2+1);

msr.TempLBuffer = zeros(1,FFTSIZE);
msr.TempRBuffer = zeros(1,FFTSIZE);
% LLbitReversalTbl;

lr2 = initLR2(msr.fs,max(40,msr.fc));
lr2.lp1_xm0 = 0;
lr2.lp2_xm0 = 0;
lr2.hp1_xm0 = 0;
lr2.hp2_xm0 = 0;
lr2.lp1_xm1 = 0;
lr2.lp2_xm1 = 0;
lr2.hp1_xm1 = 0;
lr2.hp2_xm1 = 0;

crossRat = 0.8;
msr.crossBass = max(0,min(crossRat,1));
msr.minuscrossBass = 1- msr.crossBass;

% x = x(1:msr.fs*15,:);
xl = x(:, 1); % left
xr = x(:, 2); % right

totalPCMFrameCount = length(xl);
frameCountProvided = 128;
readcount = ceil(totalPCMFrameCount/frameCountProvided);
wave_final = zeros(totalPCMFrameCount,mode+1);
pointerOffset = 0;
for i = 0 : readcount - 2
    
    pointerOffset = frameCountProvided * i;
    inputs0 = xl(pointerOffset+1:pointerOffset+frameCountProvided);
    inputs1 = xr(pointerOffset+1:pointerOffset+frameCountProvided);
    
    offset = 0;     
    while (offset < frameCountProvided)
        processing = min(frameCountProvided - offset, msr.ovpLen);
        tmp = LUTSurroundProcessSamples5_1(inputs0,inputs1,frameCountProvided);
        wave_final(pointerOffset+1:pointerOffset+frameCountProvided,:)= tmp;
        offset = offset + processing;   
    end
end

% audiowrite('sines_with_delays_lr.proced.mb.wav',wave_final,msr.fs);
audiowrite('彭芳 - 老人与海.proced.mb.wav',wave_final,msr.fs);
