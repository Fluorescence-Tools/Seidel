""" CD95 Project, AG Seidel, AG Monzel,
Nicolaas van der Voort, July 7, 2021"""
import numpy as np
import quickpbsa as pbsa
import cpp_wrappers
import os
import matplotlib.pyplot as plt

#%%
def collect_ptu(fromdirs, outdir):
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
                    plotout = None):
    time = np.arange(len(sumtrace)) * step
    #fany rescaling of axes such that curves overlap
    bg = means[0]; mf = means[1] - means[0]
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
    ax2.plot(time, fluortrace, 'r--', label='prelim Fluorophores')
    ax2.plot(time, fluortrace_final, 'k--', label='final Fluorophores')
    ax2.set_ylabel('Fluorophores', color='r')
    ax2.set_ylim(ax2_min, ax2_max)
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc=2)
    if plotout: plt.savefig(plotout, bbox_inches = 'tight', dpi = 600)
    plt.show()
    
def confocalPbsa(ffiles, Chnumbers, timestep, threshold):
    nbad = 0
    for ffile in ffiles:
        #get the channels
        channels = getTraces(ffile, Chnumbers)
        #bin the channels
        sumtrace = np.flipud(binTrace(channels, step = timestep))
        # preliminary step detection
        steppos,means,variances,posbysic,niter= pbsa.steps_preliminary.kv_single_fast(sumtrace, threshold,100)
        # call to improve_steps_single
        #sometimes the improve_steps_single finds no steps and runs into an error
        try:
            success, sicmin, steppos_out, step_out = \
                pbsa.steps_refinement.improve_steps_single(sumtrace, steppos, means, variances, posbysic)
        except IndexError:
            nbad += 1
            continue
        #calculate all the to be plotted variables
        # diffs are the numbers of frames between steps
        diffs = np.diff(np.hstack([0, steppos, len(sumtrace)]))
        # calculate mean trace 
        meantrace = np.repeat(means, diffs)
        
        # calculate fluorophores when counting every step with single occupancy
        fluors = np.cumsum(np.hstack((0, np.sign(np.diff(means)))))
        fluortrace = np.repeat(fluors, diffs)
        
        # calculate fluorophore trace
        diffs_final = np.diff(np.hstack([0, steppos_out, len(sumtrace)]))
        fluors_final = np.cumsum(np.hstack([0, step_out]))
        fluortrace_final = np.repeat(fluors_final, diffs_final)
        
        #build a file name
        filedir, file = os.path.split(ffile.decode())
        fileroot, ext = os.path.splitext(file)
        plotout = os.path.join(filedir, 'plotout', fileroot + '.png')
        
        #plot the stuff
        plotFittedTrace(sumtrace, meantrace, fluortrace, fluortrace_final, \
                        means,\
                        timestep, \
                        plotout = plotout)
        #save all the crap, what should be saved exactly?
    print('number of skipped traces: %i' % nbad)