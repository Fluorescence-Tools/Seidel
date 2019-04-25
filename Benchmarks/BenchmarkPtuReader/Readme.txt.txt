Benchmark compares performance of existing Ani image reading and exporting routines and single-purpose SplitOnTacs routine for FRC image generation.

file
N:\Singlem\singlem19-2\April\18_Calibrations_NV\10micros\10micros_1.ptu
Some 20nmCrimson bead 100% STED data

Ani file generation:
1.File is loaded into AnI
2.image size is 512*512 pixels, lines = 1
3.Tac gates from 26 to 200 is applied
4.Image is exported using the File > Export > All Images as TIFF *.tiff

SplitOnTacs routine file generation:
1.Routine genimABfromptu is executed using uselines = [1] and gate = 3228
2.imA and imB is saved in folder 10micros_1

Comparison:
1.Open images ImA.tiff, imB.tiff, Ani_export
2. Arithmetic imA+imB -> SumAB
3. Arithmatic Ani_export / 8192 -> Ani_export (This accounts for the stretching done by Ani export routine)
4.Arithmetic SumAB - Ani_export = Result
Result should be 0 for identical performance. Two deviations are observed: 
*SumAB has more photons than Ani_export, possibly due to different tac settings handling
*Ani places a photon often a pixel earlier than SplitOnTacs.
Both deviances are minor and test is considered succesfull.