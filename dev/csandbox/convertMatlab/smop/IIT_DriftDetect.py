# Generated with SMOP  0.41-beta
from libsmop import *
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_DriftDetect.m

    
@function
def IIT_DriftDetect(Im1=None,Im2=None,varargin=None,*args,**kwargs):
    varargin = IIT_DriftDetect.varargin
    nargin = IIT_DriftDetect.nargin

    #   [ShiftInteger, Shift] = IIT_DriftDetect( Im1 , Im2 )
#   Detects the drift between two squared images of the same size
    
    assert_((size(Im1,1) == size(Im2,1)) and (size(Im1,2) == size(Im2,2)),'size(Im1) ~= size(Im2)')
    assert_((size(Im1,3) == 1) and (size(Im2,3) == 1),'( size(Im1,3) || size(Im2,3) ) ~= 1')
    params.maxIterations = copy(200)
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_DriftDetect.m:15
    params.roiRadius = copy(5)
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_DriftDetect.m:16
    if (nargin > 2):
        params=IIT_read_params(params,varargin)
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_DriftDetect.m:19
    
    h,w=size(Im1,nargout=2)
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_DriftDetect.m:22
    # Crosscorrelation in the frequency domain
    Nr=dot(size(Im1,1),2) - 1
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_DriftDetect.m:25
    Nc=dot(size(Im1,2),2) - 1
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_DriftDetect.m:26
    corr=fftshift(ifft2(multiply(fft2(Im1,Nr,Nc),conj(fft2(Im2,Nr,Nc)))))
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_DriftDetect.m:27
    __,SND=max(ravel(corr),nargout=2)
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_DriftDetect.m:29
    IJ,JI=ind2sub(size(corr),SND,nargout=2)
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_DriftDetect.m:30
    ShiftInteger[1]=h - JI
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_DriftDetect.m:31
    ShiftInteger[2]=w - IJ
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_DriftDetect.m:32
    #  gaussianfitting, centered according to integerShift
    corr=real(corr)
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_DriftDetect.m:35
    corr_Size=size(corr)
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_DriftDetect.m:36
    shiftedCorr=IIT_ShiftImage(corr,ShiftInteger)
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_DriftDetect.m:38
    sub_corr=shiftedCorr(arange(floor(corr_Size(2) / 2) - params.roiRadius,floor(corr_Size(2) / 2) + params.roiRadius),arange(floor(corr_Size(1) / 2) - params.roiRadius,floor(corr_Size(1) / 2) + params.roiRadius))
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_DriftDetect.m:39
    #   Gaussian Fitting
    GaussianParameters=IIT_G2DFit_gaussian2DFittingRotational(sub_corr,params.maxIterations,1)
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_DriftDetect.m:44
    Shift[1]=params.roiRadius - GaussianParameters.ux + ShiftInteger(1) + 1
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_DriftDetect.m:46
    Shift[2]=params.roiRadius - GaussianParameters.uy + ShiftInteger(2) + 1
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_DriftDetect.m:47
    return ShiftInteger,Shift
    
if __name__ == '__main__':
    pass
    