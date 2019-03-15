function [result] = IIT_FRCAnalysis(im1, im2, varargin)
% [result] = IIT_FRCAnalysis (im1, im2, varargin)
%
% Performs FRC analysis on two images. If images are 3D, FRC is computed for each
% plane.
%
% Input arguments:
%   im1:                2D or 3D matrix containing photon counts of the first
%                       measurement
%   im2:                2D or 3D matrix containing photon counts of the second
%                       measurement
%
% Varargin:
%   pixelSize           (default 0.025 um)
%   driftCorrection:    if true, perform drift correction before FRC
%                       analysis (default false);
%   meanFilterWidth:    radius of the mean filter used to perform smooth
%                       the FRC curve (default 3 frequency bins);
%   maxIterations:      maximum iterations for the gaussian fitting
%                       performed when applying drift detection algoritm (default 200)
%   roiRadius:          radius of the crosscorrelation function considered for gaussian fitting
%                       when applying drift detection algoritm (default 5 pixels)
%   theta:              if set, for each ring split the FRC analysis for pixels in the
%                       frequency domain for which -theta < angle < +theta
%                       (0 rad). 0 <= theta <= (pi/2)
%   title:              the title of the displayed image
%
%
% Output parameters:
%   result:             struct containing all parameters and result of the
%                       FRC analysis
%
% History
% Version 1.0
% Date: 11/06/2017
%
% Authors: Giorgio Tortarolo, Marco Castello, Giuseppe Vicidomini
% Molecular Microscopy and Spectroscopy
% Istituto Italiano di Tecnologia, Genoa, Italy


%   Default parameters
params.pixelSize = 0.025;
params.driftCorrection = false;
params.meanFilterWidth = 3;
params.maxIterations = 200;
params.roiRadius = 5;
params.theta = 0;
params.title = '';

if (nargin > 2)
    params = IIT_read_params(params, varargin);
end

if (params.theta > pi/2 || params.theta < 0)
    warning('theta > pi/2 | theta < 0; forcing to 0. ');
    params.theta = 0;
end

try
    %check if the two images have same dimension and are squared
    assert( length(size(im1)) == length(size(im2)),'MATLAB:differentDimensions', 'length(size(im1)) != length(size(im2)');
    assert(nnz( size(im1) == size(im2) ) == length(size(im1)), 'MATLAB:differentSizes', 'size(im1) != size(im2)');
    [h, w, d] = size(im1);
    assert( h == w , 'MATLAB:nonSquared','h != w');
    
    im1 = double(im1);
    im2 = double(im2);
    N = size(im1,1);                        %   Number of pixels
    
    if mod(N,2)==0
        Scale = (0:(N/2-1))/(params.pixelSize*N);
    else
        Scale = (0:(N/2))/(params.pixelSize*N);
    end
    
    Scale = Scale';
    Hm =IIT_Hamming(N,N);                   %   Hamming function
    %Shift = [0 0];
    Shift = zeros(1,2,d);
    
    for cDepth = 1 : d
        tmpIm1 = im1(:,:,cDepth);
        tmpIm2 = im2(:,:,cDepth);
        
        fftA = fftshift(fft2(tmpIm1 .* Hm));
        fftB = fftshift(fft2(tmpIm2 .* Hm));
        
        if params.driftCorrection
            
            [~ , Shift(:,:,cDepth)]= IIT_DriftDetect (tmpIm1, tmpIm2, ...
                'maxIterations', params.maxIterations, 'roiRadius', params.roiRadius);                                    %Find shift
            [xF,yF] = meshgrid(-floor(h/2):-floor(h/2)+h-1,-floor(w/2):-floor(w/2)+w-1);        %Define shift in frequency domain
            fftB=fftB.*exp(-1i*2*pi.*((-xF*Shift(1))/h+(-yF*Shift(2))/w))+eps;                  %Perform the shift
            
            fprintf(strcat('\n', 'SubPixel Shift:', '\n', ...
                num2str(Shift(1)),' pixels in x;' , '\n', ...
                num2str(Shift(2)),' pixels in y.', '\n\n' ));
            
        end
        
        
        [FRC, threeSigma, fiveSigma] = IIT_performFRC(fftA, fftB, h, w,params.theta );
        
        
        sFRC.smallAngles = IIT_meanFilter(FRC.smallAngles, params.meanFilterWidth);         % smooth FRC via mean filter
        sFRC.largeAngles = IIT_meanFilter(FRC.largeAngles, params.meanFilterWidth);
        
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
        result.PixelSize = params.pixelSize;
        result. MeanFilterWidth = params.meanFilterWidth;
        result.DriftCorrected = params.driftCorrection;
        result.ThreeSigma = threeSigma;
        result.FiveSigma = fiveSigma;
        result.FixedThreshold = fixedThreshold;
        result.Theta = params.theta;
        result.CutOff = cutOff;
        result.Title = params.title;
        
        
        %plotting
        IIT_plotFRC(result);
        
        outputString = ['\n', 'Retrieved FRC resolution for ', params.title ,':', '\n'];
        if isempty(params.title)
            outputString = ['\n', 'Retrieved FRC resolution:', '\n'];
        end
        
        fprintf( ...
            [ outputString, ...
            num2str(min(round(resolution.fixedThreshold_smallAngles), ...
            round(resolution.fixedThreshold_largeAngles))),' nm.' , '\n\n'] );
        
    end
catch ME
    rethrow (ME);
    
end
end

function [cutOff, resolution] = measureResolution(Scale, cutOff, params)
fields = fieldnames(cutOff);

for i=1:numel(fields)
    if cutOff.(fields{i}) == -1
        if ( ...
                params.theta ~= 0 && ~isempty(strfind(fields{i}, '_smallAngles')) ...
                || ...
                params.theta ~= pi/2 && ~isempty(strfind(fields{i}, '_largeAngles')) ...
                )
            warning(strcat('FRC curve is not going under threshold: ', fields{i}));
        end
        cutOff.(fields{i}) = NaN;
        resolution.(fields{i}) = NaN;
    else
        resolution.(fields{i}) = 1000/Scale(cutOff.(fields{i}));
    end
end
end
