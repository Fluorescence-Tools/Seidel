# Generated with SMOP  0.41-beta
from libsmop import *
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_meanFilter.m

    
@function
def IIT_meanFilter(data=None,mean_width=None,*args,**kwargs):
    varargin = IIT_meanFilter.varargin
    nargin = IIT_meanFilter.nargin

    #   [ smoothData ] = IIT_meanFilter (data, mean_width) 	
#   Smooth data using a moving mean filter of (2 * mean_width +1) elements
    
    d=size(data,1)
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_meanFilter.m:5
    smoothData=zeros(d,1)
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_meanFilter.m:6
    for i in arange(1,d).reshape(-1):
        if i <= mean_width:
            smoothData[i]=mean(data(arange(1,i + mean_width)))
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_meanFilter.m:10
        else:
            if i > mean_width and i < d - mean_width:
                smoothData[i]=mean(data(arange(i - mean_width,i + mean_width)))
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_meanFilter.m:13
            else:
                if i >= d - mean_width:
                    smoothData[i]=mean(data(arange(i - mean_width,d)))
# K:/vanderVoortN/python/sandbox/convertMatlab/IIT_meanFilter.m:16
    
    return smoothData
    
if __name__ == '__main__':
    pass
    