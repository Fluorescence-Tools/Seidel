#this file contains arcane functions that should not exist
#author: Nicolaas van der Voort
#July 17, 2020
#HHU, Dusseldorf


import pandas as pd
import os
import aid_functions as aid
import copy

def renameDataFrameForMargarita(DataFrame):
    """takes a DataFrame and alters the column names to be Margarita readable"""
    renames = {
        'NG' : 'Green Count Rate (KHz)',
        'NR' : 'Red Count Rate (KHz)',
        'NY' : 'Yellow Count Rate (KHz)',
        'posxG' : 'peak_x_cI',
        'posyG' : 'peak_y_cI',
        'posxY' : 'peak_x_cII',
        'posyY' : 'peak_y_cII'
    }
    for oldname, newname in renames.items():
        try:
            DataFrame.rename(columns = {oldname : newname}, inplace = True)
        except:
            print('name %s not found' %oldname)
    return DataFrame

def genROIFilesFromDataFrame(dfrm, outpath, swapxy = True, 
                             appendzeros = True):
    """takes a .pg4 file and extracts data to match AnI ROI format
    AnI ROI format:
        xl # x-coord top left corner of ROI
        Lx # length of box in x
        yl # y-coord corner
        Ly # length of box in y
        x0 # x center of box
        y0 # y center of box
        BG # number of background photons integrated over the ROI
        
    Each file gets a seperate *.roi file, the start roi names should match
        filename_Green Photons.roi # for Green
        filename_Yellow Photons.roi # for Yellow"""
    length = len(dfrm)
    #do some null initialisations to get code running
    filepath = ''
    Nfiles = 0
    for i in range(length):
        row = dfrm.loc[i]
        if filepath != row['filepath']: # signals new file
            #save previous file
            file = os.path.splitext(os.path.split(filepath)[1])[0]
            if i % 100 == 0: print ('done row %i ' % i + file) # some output
            outG = os.path.join(outpath, file + '_Green Photons.roi')
            outY = os.path.join(outpath, file + '_Yellow Photons.roi')
            try: 
                for roi, out in zip([roiG, roiY], [outG, outY]):
                    #add 0 row (surens convention)
                    if appendzeros:
                        for column in ['xl', 'Lx', 'yl', 'Ly', 'x0', 'y0', 'BG']:
                            roi.at[j, column] = 0
                    savekws = {'index': False, 'header' : False, 'sep': '\t',
                               'float_format' : '%i'}
                    #there's a bug in pandas that print floats and ints the same way
                    #now the BG is rounded down.
                    roi.to_csv(out, **savekws)
            #the first time the loop is run the roiG and roiY are not defined
            except NameError: pass
            #start empty dataframe
            roiG = pd.DataFrame()
            roiY = pd.DataFrame()
            filepath = row['filepath']
            #pandas need indexing of the rows
            Nfiles +=1
            j = 0
        
        #go to next row if it is not a pair
        if row['xstartG'] == 0 or row['xstartY'] == 0: continue
    
        #build a roi row for each color
        for roi, colorid, colorname in \
            zip([roiG, roiY], ['cI', 'cII'], ['G', 'Y']):
            xcenter = round(row['peak_x_' + colorid]) #transfer to pixel number
            ycenter = round(row['peak_y_' + colorid])
            #AnI has swapped x and y convention (I guess AnI is correct)
            if swapxy:
                [ycenter, xcenter] = [xcenter, ycenter]
            bg = row[colorname + 'bg']
            #AnI has reversed
            roi.at[j, 'xl'] = xcenter - 3
            roi.at[j, 'Lx'] = 7
            roi.at[j, 'yl'] = ycenter - 3
            roi.at[j, 'Ly'] = 7
            roi.at[j, 'x0'] = xcenter
            roi.at[j, 'y0'] = ycenter
            roi.at[j, 'BG'] = 49 * bg
        j += 1
    print('Nfiles is ' + str(Nfiles))
    return 0

def movePtuWithoutSpots(filedir, dfrm, targetdir = ''):
    """some files have no spots, hence, no roifile for them is generated.
    This is a problem for AnI analysis.
    This routine moves all files in filedir which are not mentioned in dfrm,
    to a subfolder called 'nospots' such that no problem in Ani occurs
    
    targetdir kw allows for setting the relative path where the files are moved 
        to. If it is not given, '$filedir + withoutSpots' is chosen
        
    Three dirs are defined:
        filedir # the dir from where the files are being moved
        origdir"""
    filled_files = [os.path.split(fpath)[1] for fpath in dfrm['filepath']]
    if not targetdir:
        targetdir = os.path.join(filedir, 'withoutSpots')
    aid.trymkdir(targetdir)
    filenames = os.listdir(filedir)
    Nfull = 0
    for filename in filenames:
        if not filename.endswith('.ptu'): continue
        if filename in filled_files:
            Nfull += 1
        else: 
            sourcepath = os.path.join(filedir, filename)
            targetpath = os.path.join(targetdir, filename)
            os.rename(sourcepath, targetpath)
            print ('moving empty file ´%s' %filename)
    print('%i unique files have spots' % Nfull)
    
#below functions were once started, but unfinished as they seemed no longer needed
#left here for now to see if it is needed in the future still.
# def mergeAnIROIFittingWithdfrm(AnIROIFitDir, dfrm, color = 'cII'):
    # """
    # after the ROI fitting, three fit files for each file are generated.
    # One for each Green, Red and Yellow. 
    # There is a lot of redundancy in the fit information provided, it is kept
    # here for general applicability, but really creates a lot of clutter.
    # Both fileformats keep uniquely: the name of the original .ptu file
                                    # the center of the ROI
    # These identifiers are used to link the two together. One is not enough,
    
    # AnIROIFitDir:   the directory where AnI outputs its batch fitting
                    # It contains a Green, Red and Yellow file for each ptu file.
    # dfrm:   DataFrame where you want to be appended to, should contain filename
            # and roi center information. Analysed Seidel data usually does.
    # color:  The color used to generate the ROI center. It is advised to choose 
            # yellow because it is more precise. Seidel uses different ROIs for
            # each color channel, but that is thrown away here for simplicity.
    # """
    # #don't modify what you're iterating over
    # dfrm_cpy = copy.deepcopy(dfrm)
    # #get list of file
    # ARFdir = AnIROIFitDir
    # files = [file for file in os.listdir(ARFdir) if ispx4(file)][:10]
    # #load file as Dataframe
    # for file in files:
        # ARDdfrm = pd.read_csv(os.path.join(ARFdir, file), sep = '\t')
        # #loop over rows in ARDdfrm
        # for ARDindex, ARDrow in ARDdfrm.iterrows():
            # #skip emty rows (AnI includes them, but they are not needed):
            # if ARDrow['Duration (ms)'] == 0: continue
            # #find corresponding row in original dfrm
            # charid = file.find('_Frames') #some arcane shit that allows me to strip the file to its original length
            # basefile = file[:charid]
            # rowids = dfrm['filepath'].str.contains(basefile)
            # #loop over the rows untill one row with matching center is found
            # for index, row in dfrm[rowids].iterrows():
                # # again swapping of the axes
                # if round(row['peak_x_' + color]) == ARDrow['Y pixel'] and \
                    # round(row['peak_y_' + color]) == ARDrow['X pixel']:
                    # # append to row to dfrm_cpy
                    # print('found matching rowid %i' % index)
            # #Q: how to deal with rows that share the same name?
    # #return drfm_copy
# def ispx4(filename):
    # return '.pg4' in filename or \
            # 'pr4' in filename or \
            # 'py4' in filename
    
