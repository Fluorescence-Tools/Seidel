function [result] = IIT_FRCAnalysis_C(im1, im2)

%   Default parameters
pixelSize = 0.025;
meanFilterWidth = 3;
maxIterations = 200;
roiRadius = 5;
theta = 0;

im1 = double(im1);
im2 = double(im2)

N = size(im1,1);                        %   Number of pixels;

if mod(N,2)==0
    Scale = (0:(N/2-1))/(pixelSize*N);
else
    Scale = (0:(N/2))/(pixelSize*N);
end

Scale = Scale';
Hm =IIT_Hamming(N,N);                   %   Hamming function
%Shift = [0 0];
Shift = zeros(1,2,d);

fftA = fftshift(fft2(im1 .* Hm));
fftB = fftshift(fft2(im2 .* Hm));

[FRC, threeSigma, fiveSigma] = IIT_performFRC(fftA, fftB, h, w, theta );


sFRC.smallAngles = IIT_meanFilter(FRC.smallAngles, meanFilterWidth);         % smooth FRC via mean filter
sFRC.largeAngles = IIT_meanFilter(FRC.largeAngles, meanFilterWidth);

fixedThreshold = 1/7 * ones(size(Scale,1),1);

% CutOff frequencies for different threshold criteria - build a
% function
cutOff.fixedThreshold_smallAngles = IIT_findIntercept(sFRC.smallAngles, fixedThreshold);
cutOff.threeSigma_smallAngles = IIT_findIntercept(sFRC.smallAngles, threeSigma);
cutOff.fiveSigma_smallAngles = IIT_findIntercept(sFRC.smallAngles, fiveSigma);
cutOff.fixedThreshold_largeAngles = IIT_findIntercept(sFRC.largeAngles, fixedThreshold);
cutOff.threeSigma_largeAngles = IIT_findIntercept(sFRC.largeAngles, threeSigma);
cutOff.fiveSigma_largeAngles = IIT_findIntercept(sFRC.largeAngles, fiveSigma);

[cutOff, resolution] = measureResolution(Scale, cutOff, params);

result.FRC = FRC;
result.sFRC = sFRC;
result.Resolution = resolution;
result.Im1 = im1;
result.Im2 = im2;
result.Scale = Scale ;
result.N = N;
result.Shift = Shift;
result.PixelSize = pixelSize;
result. MeanFilterWidth = meanFilterWidth;
result.DriftCorrected = driftCorrection;
result.ThreeSigma = threeSigma;
result.FiveSigma = fiveSigma;
result.FixedThreshold = fixedThreshold;
result.Theta = theta;
result.CutOff = cutOff;
result.Title = title;
end