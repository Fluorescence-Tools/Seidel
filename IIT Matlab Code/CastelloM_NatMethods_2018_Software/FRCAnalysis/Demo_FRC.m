%%   load Data from example directory

selpath = uigetdir('Select the folder containing FRC Analysis data for test:');
row = 400;
col = 400;

separator = '\';
if (ismac)
    separator = '/';
end

fin=fopen([selpath, separator, 'Measure1_CLSM_0_2_AU.raw'], 'r');
I=fread(fin,row*col,'uint16=>uint16'); 
Z=double(reshape(I,row,col));
Z=Z';
Measure1.CLSM02 = Z';

fin=fopen([selpath, separator, 'Measure1_CLSM_1_4_AU.raw'], 'r');
I=fread(fin,row*col,'uint16=>uint16'); 
Z=double(reshape(I,row,col));
Z=Z';
Measure1.CLSM14 = Z';

fin=fopen([selpath, separator, 'Measure1_ISM.raw'], 'r');
I=fread(fin,row*col,'uint16=>uint16'); 
Z=double(reshape(I,row,col));
Z=Z';
Measure1.ISM = Z';

fin=fopen([selpath, separator, 'Measure2_CLSM_0_2_AU.raw'], 'r');
I=fread(fin,row*col,'uint16=>uint16'); 
Z=double(reshape(I,row,col));
Z=Z';
Measure2.CLSM02 = Z';

fin=fopen([selpath, separator, 'Measure2_CLSM_1_4_AU.raw'], 'r');
I=fread(fin,row*col,'uint16=>uint16'); 
Z=double(reshape(I,row,col));
Z=Z';
Measure2.CLSM14 = Z';

fin=fopen([selpath, separator, 'Measure2_ISM.raw'], 'r');
I=fread(fin,row*col,'uint16=>uint16'); 
Z=double(reshape(I,row,col));
Z=Z';
Measure2.ISM = Z';

cd(selpath);
cd .. 
cd ..
warning('off','MATLAB:rmpath:DirNotFound');
rmpath(genpath(cd));

addpath([cd,separator,'FRCAnalysis/Functions']);

clear Z I row col fin selpath separator
%% perform FRC analysis

FRCresult.CLSM02 = IIT_FRCAnalysis(Measure1.CLSM02, Measure2.CLSM02, 'pixelSize',0.0376, 'title', 'CLSM 0.2 A.U.'); 
FRCresult.CLSM14 = IIT_FRCAnalysis(Measure1.CLSM14, Measure2.CLSM14, 'pixelSize',0.0376, 'title', 'CLSM 1.4 A.U.'); 
FRCresult.ISM = IIT_FRCAnalysis(Measure1.ISM, Measure2.ISM, 'pixelSize',0.0376, 'title', 'ISM'); 

