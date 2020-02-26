import os,sys,json,array
import numpy as np
import ROOT as rt
from math import sqrt

from oscillate_flux_sterile import read_flux, oscillate_flux_sterile

def define_axes( nbinsx, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax ):
    xbins = np.logspace(np.log10(xmin),np.log10(xmax), nbinsx)
    ybins = np.logspace(np.log10(ymin),np.log10(ymax), nbinsy)
    zbins = np.logspace(np.log10(zmin),np.log10(zmax), zbinsy)    
    return xbins,ybins

def make_histogram( histname, nbinsx, xmin, xmax, nbinsy, ymin, ymax ):
    xarr = array.array( 'f', [0.0]*(nbinsx+1) )
    yarr = array.array( 'f', [0.0]*(nbinsy+1) )

    xaxis,yaxis = define_axes( nbinsx, xmin, xmax, nbinsy, ymin, ymax )
    for xb in xrange(nbinsx):
        xarr[xb] = xaxis[xb]
    xarr[nbinsx] = xmax
    print xaxis
    print yaxis

    for yb in xrange(nbinsy):
        yarr[yb] = yaxis[yb]
    yarr[nbinsy] = ymax

    hist = rt.TH2F( histname, "", nbinsx, xarr, nbinsy, yarr )
    return hist

def make_param_array( nbinsU, minU, maxU, nbinsdm2, mindm2, maxdm2 ):

    uaxis  = np.logspace(np.log10(minU),  np.log10(maxU),   nbinsU)
    dmaxis = np.logspace(np.log10(mindm2),np.log10(maxdm2), nbinsdm2)
    return (dmaxis,uaxis,uaxis,uaxis)

def make_flux_file( filename, oscflux ):
    fluxout = open(filename,'w')
    for ibin in xrange(oscflux.shape[0]):
        print>>fluxout,oscflux[ibin,0]," ",oscflux[ibin,1]," ",oscflux[ibin,2]," ",oscflux[ibin,3]," ",oscflux[ibin,4]," ",oscflux[ibin,5]," ",oscflux[ibin,6]
    fluxout.close()

def get_channel_data( chanfile ):
    """ collect channel information from files channels folder """
    
    channeldata = { "channame":[],
                    "channum":[],
                    "cpstate":[],
                    "flavor":[],
                    "num_target_factor":[],
                    "numchans":0}

    chanfilename = "../../channels/channels_"+chanfile+".dat";
    if not os.path.exists(chanfilename):
        raise ValueError( "could not find channel file: ",chanfilename)
    
    chans = open(chanfilename,'r')
    for l in chans.readlines():
        l = l.strip()
        info = l.split()
        channeldata["channame"].append( info[0] )
        channeldata["channum"].append( int(info[1]) )
        channeldata["cpstate"].append( info[2] )
        channeldata["flavor"].append(  info[3] )
        channeldata["num_target_factor"].append( float(info[4]) )
        channeldata["numchans"] += 1
    chans.close();

    print("Number of channels: ",channeldata["numchans"])
    return channeldata

    
def read_output_data( fluxname, chanfile, exptconfig ):

    channeldata = get_channel_data( chanfile )

    # Now get the info from the files.  Assumption is that all have the same binning.
    nchannels = channeldata["numchans"]
    maxpoints = 1000

    # setup data arrays
    import numpy as np
    data = { "total_events":np.zeros( maxpoints ),
             "total_events_smeared":np.zeros( maxpoints ),
             "en":np.zeros( maxpoints ),
             "eng":np.zeros( maxpoints ),
             "events":np.zeros( (nchannels,maxpoints) ),
             "events_smeared":np.zeros( (nchannels,maxpoints) ),
             "nue_es_smeared":np.zeros( maxpoints ),
             "nuebar_es_smeared":np.zeros( maxpoints ),
             "nux_es_smeared":np.zeros( maxpoints ),
             "channeldata":channeldata }

    # read in data
    for ifile in xrange(nchannels):

        # read unsmeared data        
        filename = "out/"+fluxname+"_"+channeldata["channame"][ifile]+"_"+exptconfig+"_events.dat"
        f = open(filename,'r')
        ll = f.readlines()
        npoints = 0
        #print("parsing ",filename)
        for i,l in enumerate(ll):
            l = l.strip()

            if "---" in l or "Total" in l:
                continue
            
            info = l.split()            
            data["en"][i]  = float(info[0])*1000.0
            data["eng"][i] = float(info[0])*1000.0
            nevents = float(info[1])
            if abs(nevents)>1e9 or info[1].strip()=="NaN":
                nevents = 0.0
            data["events"][ ifile, i ] = nevents

            # Account for the number of targets relative to reference target-- for unweighted file only
            #data["events"][ ifile, i ] *= channeldata["num_target_factor"][ifile]

            data["total_events"][i] += data["events"][ifile,i]
            npoints += 1
        #print(" found %d points for channel %d. integral=%f"%(npoints,ifile,np.sum( data["events"][ifile,:] ) ))
        f.close()

        # read smeared data
        filename_smeared = "out/"+fluxname+"_"+channeldata["channame"][ifile]+"_"+exptconfig+"_events_smeared.dat";
        f = open(filename_smeared,'r')
        ll = f.readlines()
        npoints_smeared = 0        
        for i,l in enumerate(ll):
            l = l.strip()

            if "---" in l or "Total" in l:
                continue
            
            info = l.split()
            nevents = float(info[1])
            if abs(nevents)>1e9 or info[1].strip()=="NaN":
                data["events_smeared"][ ifile, i ] = 0.0
            else:
                data["events_smeared"][ ifile, i ] = nevents

            # Account for the number of targets relative to reference target-- for unweighted file only
            #data["events_smeared"][ ifile, i ] *= channeldata["num_target_factor"][ifile]

            data["total_events_smeared"][i] += data["events_smeared"][ifile,i]
            npoints_smeared += 1

            for es_chans in ["nue_e","nuebar_e","numu_e","numubar_e","nutau_e","nutaubar_e"]:
                if channeldata["channame"][ifile] == es_chans:
                    if es_chans in ["nue_e","nuebar_e"]:
                        data["%ss_smeared"%(es_chans)][i] += data["events_smeared"][ifile,i]
                    else:
                        data["nux_es_smeared"][i] += data["events_smeared"][ifile,i]
            
        #print(" found %d points for smeared channel %d; integral=%f"%(npoints_smeared,ifile, np.sum( data["events_smeared"][ifile,:] ) ))

    return data

def ana_data( oscdata, nulldata ):
    """
[ 6   nue_Ar40_marley1 ]  327.01578253498656
[ 7   nuebar_Ar40 ]  0.0
[ 8   nc_nue_Ar40 ]  26.45700936979586
[ 9   nc_numu_Ar40 ]  10.90995486894403
[ 10   nc_nutau_Ar40 ]  0.006368001369250145
[ 11   nc_nuebar_Ar40 ]  0.009812826296037637
[ 12   nc_numubar_Ar40 ]  44.057154025162575
[ 13   nc_nutaubar_Ar40 ]  0.009812826296037637
    """
    totnc_null = 0.0
    totnc_osc  = 0.0
    for fid in xrange(8,14):
        totnc_null += np.sum(nulldata["events"][fid,:])
        totnc_osc  += np.sum(oscdata["events"][fid,:])
    diff_nc = totnc_osc-totnc_null
    onebin_nc_chi2 = diff_nc*diff_nc/totnc_null

    cc_null = nulldata["events"][6,:] + nulldata["events"][7,:]
    cc_osc  = oscdata["events"][6,:]  + oscdata["events"][7,:]
    diff    = cc_osc-cc_null
    diff2   = np.power(diff,2)
    chi2    = diff2[cc_null>0]/cc_null[cc_null>0]
    ndf     = len(cc_null[cc_null>0])

    print "cc_osc: ",cc_osc[cc_null>0]
    print "cc_null: ",cc_null[cc_null>0]
    print "diff: ",diff[cc_null>0]
    print "chi2: ",chi2,ndf
    schi2 = (np.sum(chi2) + (totnc_osc-totnc_null)*(totnc_osc-totnc_null)/totnc_null)/float(ndf+1-1)
    print "sum(diff): ",np.abs(diff[cc_null>0]).sum()
    print "sum(diff2):",np.sum(diff2[cc_null>0])
    print "sum(chi2):",np.sum(chi2)
    print "chi2/ndf: ",schi2

    tot_cc_osc = np.sum(cc_osc[cc_null>0])
    tot_cc_nul = np.sum(cc_null[cc_null>0])
    tot_cc_diff = tot_cc_osc-tot_cc_nul
    print "tot_cc(null)=",tot_cc_nul
    print "tot_cc(osc)=",tot_cc_osc
    print "onebin_diff = ",tot_cc_diff
    onebin_chi2 = tot_cc_diff*tot_cc_diff/tot_cc_nul
    print "chi2(one-bin) = ",onebin_chi2," + ",onebin_nc_chi2
        
    return 0
    

if __name__ == "__main__":
    # we want to generate a sensitivity plot
    # it's a log-log plot, so we need to define the bins ourselves

    (dm2_axis, ue4_axis, um4_axis, ut4_axis) = make_param_array( 10, 0.1, 1.0, 1, 7.0, 7.0 )
    flux_unosc = read_flux( "fluxes/stpi.dat" )

    channame    = "argon_marley1"
    expt_config = "ar40kt"

    pnull = os.popen("./supernova.pl stpi %s %s 0"%(channame,expt_config))
    lines = pnull.readlines()
    nulldata = read_output_data("stpi","argon_marley1","ar40kt")
    print "[ NULL HYPO ]"
    for ic,chan in enumerate(nulldata["channeldata"]["channame"]):
        nulldata["events"][ic,:] *= 3.0*3.14e7*0.000612/40.0*(20.0*20.0)/(29.0*29.0)*(5000.0/8766.0)
        print "[",ic," ",chan,"] ",np.sum(nulldata["events"][ic,:])

    print "[Start osc calcs]"
    for dm2 in dm2_axis:
        for ue4 in ue4_axis:
            for um4 in um4_axis:
                for ut4 in ut4_axis:
                    oscpars = {"dm2":dm2,
                               "Ue4":sqrt(ue4),
                               "Um4":sqrt(um4),
                               "Ut4":sqrt(ut4),
                               "L_m":29.0}
                    sin2_mue = 4.0*ue4*um4
                    print "[run osc: sin2_{mue}]",sin2_mue
                    oscflux = oscillate_flux_sterile( flux_unosc, oscpars )
                    make_flux_file( "fluxes/osc_stpi.dat", oscflux )
                    #pline = os.popen("./supernova.pl osc_stpi %s %s 0"%(channame,expt_config))
                    #lines = pline.readlines()
                    os.system("./supernova.pl osc_stpi %s %s 0"%(channame,expt_config))
                    #for l in lines:
                    #    print l.strip()
                    data = read_output_data("osc_stpi","argon_marley1","ar40kt")
                    for ic,chan in enumerate(data["channeldata"]["channame"]):
                        data["events"][ic,:] *= 3.0*3.14e7*0.000612/40.0*(20.0*20.0)/(29.0*29.0)*(5000.0/8766.0)
                        print "[",ic," ",chan,"] ",np.sum(data["events"][ic,:])
                    ana_data( data, nulldata )
                    break
                break
            break
        break
                    
    print "dm2 axis: ",dm2_axis
    print "ue4 axis: ",ue4_axis

    
    pass
