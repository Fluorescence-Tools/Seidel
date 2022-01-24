""" CD95 Project, AG Seidel, AG Monzel,
Nicolaas van der Voort, July 7, 2021"""
import numpy as np
import quickpbsa as pbsa
import cpp_wrappers
import os
import matplotlib.pyplot as plt
import copy
import pickle
import aid_functions as aid

#for debugging
import traceback

#%%
def collect_ptu(fromdirs, outdir):
    aid.trymkdir(outdir)
    for fromdir in fromdirs:
        files = os.listdir(fromdir)
        for file in files:
            if file.endswith('.ptu'):
                ffile = os.path.join(fromdir, file)
                outfile = os.path.join(outdir, file)
                os.rename(ffile, outfile)
                
def getTraces(fname, Chnumbers, counttime = 25e-9):
    NumRecords = cpp_wrappers.ptuHeader_wrap (fname)
    eventN, tac, t, can = cpp_wrappers.ptu_wrap(fname, NumRecords)
    channels = {}
    for Chnumber in Chnumbers:
        channels['ch'+str(Chnumber)] = t[can==Chnumber].astype(np.float64)
    startt = t[0]
    for ch in channels.values():
        ch -= startt
        ch *= counttime
    return channels
def binTrace(channels, step = 5e-3):
    binneddata = []
    maxt = 0
    for ch in channels.values():
        maxt = max(maxt, max(ch))
    xlim = (0, maxt)
    bin_edges = np.arange(*xlim,step)
    bincenters = bin_edges + step / 2.
    binneddata.append(bincenters[:-1])
    header = 'time[s]\t'
    for name, ch in channels.items():
        y, xborders = np.histogram(ch, bin_edges)
        binneddata.append(y.astype(np.float))
        header += name + '\t'
        if 'sumy' not in locals(): 
            sumy = y
        else: 
            sumy +=y
    return sumy

def plotFittedTrace(sumtrace, meantrace, fluortrace, fluortrace_final, means, step, \
                    plotout = None, datarange = None):
    #zoom in on the start of the trace to visualize quick bleacing events
    if datarange:
        sumtrace = sumtrace[-datarange:]
        meantrace = meantrace[-datarange:]
        fluortrace = fluortrace[-datarange:]
        fluortrace_final = fluortrace_final[-datarange:]
    time = np.arange(len(sumtrace)) * step
    #fany rescaling of axes such that curves overlap
    bg = min(means); 
    try: 
        mf = abs(means[1] - means[0])
    except IndexError:
        mf = 1
    ax1_max = max(sumtrace) * 1.05
    ax2_max = (ax1_max - bg) / mf
    ax2_min = - bg / mf
    fig = plt.figure()
    plt.plot(time, sumtrace, label='Trace')
    plt.plot(time, meantrace, label='Means')
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ax1.set_xlabel('time(s)')
    ax1.set_ylim(0, ax1_max)
    ax1.set_ylabel('photon counts / %.0f ms' % (step * 1e3))
    ax2.plot(time, fluortrace, 'r-.', label='prelim Fluorophores')
    ax2.plot(time, fluortrace_final, 'k--', label='final Fluorophores')
    ax2.set_ylabel('Fluorophores', color='r')
    ax2.set_ylim(ax2_min, ax2_max)
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc=2)
    if plotout: plt.savefig(plotout, bbox_inches = 'tight', dpi = 600)
    plt.show()
    
def tracesFromPtu(ffiles, Chnumbers, timestep):
    sumtraces = []
    for ffile in ffiles:
        #get the channels
        channels = getTraces(ffile, Chnumbers)
        #bin the channels
        sumtraces.append(np.flipud(binTrace(channels, step = timestep)))
    return sumtraces
    
def confocalPbsa(sumtraces, threshold, ffiles, timestep, Chnumbers, \
        lstout = None, verbose = True):
    nbad = 0
    traceLst = []
    for sumtrace, ffile in zip(sumtraces, ffiles):
        print(ffile)
        # preliminary step detection
        steppos, means, variances, posbysic, niter = \
            pbsa.steps_preliminary.kv_single_fast(sumtrace, threshold,100)

        #calculate all the to be plotted variables
        # diffs are the numbers of frames between steps
        diffs = np.diff(np.hstack([0, steppos, len(sumtrace)]))
        # calculate mean trace 
        meantrace = np.repeat(means, diffs)
        
        # calculate fluorophores when counting every step with single occupancy
        fluors = np.cumsum(np.hstack((0, np.sign(np.diff(means)))))
        #NV added
        #adress a bug where at the end of a trace there are still
        #fluorophores. This messes with the zero level.
        fluors -= min(fluors) 
        #end added
        fluortrace_prelim = np.repeat(fluors, diffs)
        
        #build a file name
        filedir, file = os.path.split(ffile.decode())
        fileroot, ext = os.path.splitext(file)
        plotoutdir = os.path.join(filedir, 'plotout')
        plotout = os.path.join(plotoutdir, fileroot +'thr%.0f' % threshold + '.svg')
        aid.trymkdir(plotoutdir)
        # call to improve_steps_single
        #sometimes the improve_steps_single finds no steps and runs into an error
        
        try:
            success, sicmin, steppos_out, step_out = \
                pbsa.steps_refinement.improve_steps_single(sumtrace, steppos, means, variances, posbysic)
            
        except IndexError as e:
            # if verbose
            # # plot the curve that was unsuccessfull. try to troubleshoot the problem
            # time = np.arange(len(sumtrace)) * timestep
            # plt.plot(time, sumtrace)
            # plt.plot(time, meantrace)
            # plt.title('refinement step unsuccesfull for %s' % file)
            # twostepoutdir = os.path.join(filedir, 'twostepplots')
            # twostepplotout = os.path.join(twostepoutdir, fileroot + 'png')
            # aid.trymkdir(twostepoutdir)
            # plt.savefig(twostepplotout, dpi = 600, bbox_inches = 'tight')
            # plt.show()
            nbad += 1
            success = False
            sicmin = -1
            steppos_out = steppos # take the pos from the preliminary step
            step_out = np.cumsum(np.sign(np.diff(means))) # take the means as step values
        # calculate fluorophore trace
        diffs_final = np.diff(np.hstack([0, steppos_out, len(sumtrace)]))
        fluors_final = np.cumsum(np.hstack([0, step_out]))
        #address a bug where at the end of a trace there are still fluorophores
        fluors_final -= min(fluors_final) # set minimal step to 0
        fluortrace_final = np.repeat(fluors_final, diffs_final)
        if verbose:
            #plot the stuff
            plotFittedTrace(sumtrace, meantrace, fluortrace_prelim, fluortrace_final, \
                            means,\
                            timestep, \
                            plotout = plotout)
        #save all the crap, what should be saved exactly?
        #all these if statements make this such a mess, is there a better approach?
        #here al the single parameter properties come
        onep = {}
        if len(means) == 1:
            onep['Nfluors_kv'] = 0
        else:
            onep['Nfluors_kv'] = max(fluors)
        
        onep['Nfluors_final'] = max(fluortrace_final)
        onep['Nsteps_final'] = len(fluors_final) -1
        #how many steps were larger than 1 fluo and by how much
        onep['Nmultisteps'] = np.sum(np.abs(step_out) -1)
        #workaround needed in case no step is found
        if len(means) < 2:
            stepsize = -1
            stepvariance = -1
        else:
            stepsize = means[1]-means[0]
            stepvariance = variances[1] - variances[0]
        onep['stepsize'] = stepsize
        onep['stepvariance'] = stepvariance
        onep['bgsize'] = means[0]
        onep['bgvariance'] = variances[0]
        onep['niter'] = niter
        onep['timestep'] = timestep
        onep['fitsucces'] = success
        onep['sicmin'] = sicmin
        onep['threshold'] = threshold
        onep['Nfluors_end_kv'] = fluors[0]
        onep['Nfluors_end_refined'] = fluors_final[0]
        #here all the multiparemeter properties come
        multip = {}
        multip['steppos_kv'] = steppos
        multip['means'] = means
        multip['diffs_kv'] = np.diff(means)
        multip['meantrace'] = meantrace
        multip['variances'] = variances
        multip['posbysic'] = posbysic
        multip['trace'] = sumtrace
        multip['Chnumbers'] = Chnumbers
        multip['steppos_out'] = steppos_out
        multip['step_out'] = step_out
        multip['tracefit'] = fluortrace_final * stepsize + means[0]
        multip['fluortrace_kv'] = fluortrace_prelim
        multip['fluortrace_final'] = fluortrace_final
        multip['file'] = file
        traceinfo = {'onep': onep, 'multip':multip}
        traceLst.append(traceinfo)
    print('number of traces with one or less steps: %i' % nbad)
    #save if identifier is given
    if lstout:
        aid.trymkdir(os.path.split(lstout)[0])
        aid.savepickle(traceLst, lstout)
        #save csv of 1 parameter properties
        csvout = os.path.splitext(lstout)[0] + '.csv'
        printOnepToCsv(traceLst, csvout)
    return traceLst
    
def filterTraceProperty(traceLst, property, minval = None, maxval = None):
    """filter a trace List based on some property.
    property is a string and should represent a single value property"""
    #we don't want to edit the original
    traceLst = copy.deepcopy(traceLst)
    badLst = []
    for i, trace in enumerate(traceLst):
        val = trace['onep'][property]
        if minval is not None:
            if val <= minval:
                badLst.append(i)
                continue
        if maxval is not None:
            if val >= maxval:
                badLst.append(i)
                continue
    #pop last to first
    for i in badLst[::-1]:
        traceLst.pop(i)
    return traceLst
def getTraceProperty(traceLst, property, oneormulti = 'onep'):
    out = []
    for trace in traceLst:
        out.append(trace[oneormulti][property])
    return out
def plotTraceLst(traceLst, outfolder = None, **kwargs):
    for trace in traceLst:
        fluortrace_prelim = trace['multip']['fluortrace_kv']
        sumtrace = trace['multip']['trace']
        meantrace = trace['multip']['meantrace']
        fluortrace_final = trace['multip']['fluortrace_final']
        means = trace['multip']['means']
        timestep = trace['onep']['timestep']
        if outfolder: plotout = os.path.join(outfolder, \
            trace['multip']['file'][:-4] + '.svg')
        else: plotout = None
        print(trace["multip"]['file'][-20:])
        plotFittedTrace(sumtrace, meantrace, fluortrace_prelim, fluortrace_final, \
                            means,\
                            timestep, \
                            plotout = plotout, \
                            **kwargs)
                            
def printOnepToCsv(traceLst, csvout = None):
    import pandas as pd
    oneponly = [el['onep'] for el in traceLst]
    outdfrm = pd.DataFrame.from_dict(oneponly)
    outdfrm.to_csv(csvout)
    return outdfrm
    
def calcVarAlpha(alpha, n):
    """
    calculate the variance on alpha assuming
    alpha is the average amount of time in the on statements
    n is the average amount of transitions in selected time period
    Formula originates from FRET lines paper from Oleg"""
    #alpha = <kon / (kon + koff)>
    #n = <kon + koff> * tbin
    return alpha * (1-alpha) * 2 / n * (1 + (np.exp(-n)-1)/n)
def calcVarSOleg(Nfl, alpha, n):
    """using standard propagation of the variance for formula
    S = Nfl * alpha
    we obtain this formula, however it does not match MC simulations"""
    return Nfl**2 * calcVarAlpha(alpha, n) + alpha**2*Nfl

def calcVarSOleg2(Nfl, alpha, n):
    """simular to calcVarSOleg
    the square in the second therm is removed, then this formula matches MC
    simulations"""
    return Nfl**2 * calcVarAlpha(alpha, n) + alpha*Nfl