import os,sys,time
import numpy as np
import ROOT as rt

def read_flux( fluxfile, livetime_hours=5000.0, distance=29.0, normalize=False ):
    """
    import flux file into numpy array
    """
    s = time.time()
    f = open(fluxfile,'r')
    fl = f.readlines()
    flux = np.zeros((len(fl), 7 ))
    for i,l in enumerate(fl):
        l = l.strip()
        info = l.split()
        for j in xrange(7):
            flux[i,j] = float(info[j])
    if normalize:
        flux[:,1:] *= 3.14e7*(20.0*20.0)/(distance*distance)*(livetime_hours/8766.0)
    print "time to load flux from file: %.3f secs"%(time.time()-s)        
    return flux
        
    
def oscillate_flux_sterile( flux, params ):
    """
    flux[ energy, nue, numu, nutau, nuebar numubar, nutaubar ]
    """
    s = time.time()
    dm2  = params["dm2"]
    L    = params["L_m"]*0.001
    sin2dm_arg = 1.27*dm2*L/flux[:,0]
    sin2dm = np.sin( sin2dm_arg )
    sin2dm = np.power( sin2dm, 2 )
    
    if "Ue4" in params and "Um4" in params and "Ut4" in params:
        # use abs value of matrix values
        Ue4 = params["Ue4"]
        Um4 = params["Um4"]
        Ut4 = params["Ut4"]
        prob_ee = 1.0 - 4.0*(1-Ue4*Ue4)*(Ue4*Ue4)*sin2dm
        prob_mm = 1.0 - 4.0*(1-Um4*Um4)*(Um4*Um4)*sin2dm
        prob_tt = 1.0 - 4.0*(1-Ut4*Ut4)*(Ut4*Ut4)*sin2dm
        prob_em = 4*(Ue4*Ue4)*(Um4*Um4)*sin2dm
        prob_et = 4*(Ue4*Ue4)*(Ut4*Ut4)*sin2dm
        prob_me = 4*(Um4*Um4)*(Ue4*Ue4)*sin2dm
        prob_mt = 4*(Um4*Um4)*(Ut4*Ut4)*sin2dm
        prob_te = 4*(Ut4*Ut4)*(Ue4*Ue4)*sin2dm
        prob_tm = 4*(Ut4*Ut4)*(Um4*Um4)*sin2dm

        
    elif "sin2_{ee}" in params and "sin2_{mm}" in params and "sin2_{me}" in params:
        s2ee = params["sin2_{ee}"]
        s2mm = params["sin2_{mm}"]
        s2me = params["sin2_{me}"]
        prob_ee = 1.0-s2ee*sin2dm
        prob_mm = 1.0-s2mm*sin2dm
        prob_me = s2me*sin2dm
        prob_em = prob_me
        

    oscflux = np.zeros( flux.shape )

    #print "arg: ",sin2dm_arg[:10]
    #print "E: ",flux[:10,0]
    #print "sin2dm: ",sin2dm[:10]

    # copy energies
    oscflux[:,0] = flux[:,0]
    
    # nue
    oscflux[:,1] = flux[:,1]*prob_ee[:] + flux[:,2]*prob_me
    # numu
    oscflux[:,2] = flux[:,2]*prob_mm[:] + flux[:,1]*prob_em
    # tau
    oscflux[:,3] = flux[:,1]*prob_et + flux[:,2]*prob_mt
    # nue-bar
    oscflux[:,4] = flux[:,4]*prob_ee[:] + flux[:,5]*prob_me
    # numu-bar
    oscflux[:,5] = flux[:,5]*prob_mm[:] + flux[:,4]*prob_em
    # nutar-bar
    oscflux[:,6] = flux[:,4]*prob_et + flux[:,5]*prob_mt

    print "osc applied for {} in {} secs".format( params, time.time()-s )
        
    return oscflux

def make_flux_tgraph( flux ):
    npts = flux.shape[0]
    tgraph_v = {"nue":rt.TGraph(npts),
                "numu":rt.TGraph(npts),
                "nutau":rt.TGraph(npts),
                "nuebar":rt.TGraph(npts),
                "numubar":rt.TGraph(npts),
                "nutaubar":rt.TGraph(npts)}

    tots = {}    
    for k in tgraph_v.keys():
        tots[k] = 0.0
        
    for i in xrange(npts):
        tgraph_v["nue"].SetPoint(i, flux[i,0],flux[i,1])
        tgraph_v["numu"].SetPoint(i,flux[i,0],flux[i,2])
        tgraph_v["nutau"].SetPoint(i,flux[i,0],flux[i,3])        
        tgraph_v["nuebar"].SetPoint(i, flux[i,0],flux[i,4])
        tgraph_v["numubar"].SetPoint(i,flux[i,0],flux[i,5])
        tgraph_v["nutaubar"].SetPoint(i,flux[i,0],flux[i,6])
        tots["nue"] += flux[i,1]
        tots["numu"] += flux[i,2]
        tots["nutau"] += flux[i,3]
        tots["nuebar"] += flux[i,4]
        tots["numubar"] += flux[i,5]
        tots["nutaubar"] += flux[i,6]        

    return tgraph_v,tots
        

if __name__ == "__main__":

    flux = read_flux( "stpi.dat" )

    #params = { "sin2_{ee}":0.01,
    #           "sin2_{mm}":0.01,
    #           "sin2_{me}":0.01,
    #          "dm2":7.0,
    #          "L_m":29.0 }
    params = { "Ue4":0.1,
               "Um4":0.1,
               "Ut4":0.1,
               "dm2":7.0,
               "L_m":29.0 }

    print "given pars: sin2_{ee}: ",4.0*(1-params["Ue4"]*params["Ue4"])*params["Ue4"]*params["Ue4"]
    print "given pars: sin2_{mm}: ",4.0*(1-params["Um4"]*params["Um4"])*params["Um4"]*params["Um4"]
    print "given pars: sin2_{tt}: ",4.0*(1-params["Ut4"]*params["Ut4"])*params["Ut4"]*params["Ut4"]    
    print "given pars: sin2_{me}: ",4.0*(params["Ue4"]*params["Ue4"])*(params["Um4"]*params["Um4"])

    oscflux = oscillate_flux_sterile( flux, params )

    gflux,totflux = make_flux_tgraph(flux)
    oflux,totosc  = make_flux_tgraph(oscflux)

    print oscflux[0:10,:]

    c = rt.TCanvas("c","c",1200,600)
    c.Draw()

    for g in [gflux,oflux]:
        g["numu"].SetLineColor(rt.kRed)
        g["nue"].SetLineColor(rt.kBlue)
        g["nutau"].SetLineColor(rt.kGreen+3)        
        g["numubar"].SetLineColor(rt.kMagenta)
        g["nuebar"].SetLineColor(rt.kCyan)
        g["nutaubar"].SetLineColor(rt.kTeal+10)        

    for n,g in oflux.items():
        g.SetLineStyle(2)

    gflux["numu"].Draw("AL")
    for n,g in gflux.items():
        g.Draw("Lsame")
    for n,g in oflux.items():
        g.Draw("Lsame")
    #gflux["numu"].GetXaxis().SetRangeUser(0,0.060)

    for fl in ["nue","numu","nutau","nuebar","numubar","nutaubar"]:
        print "total {}: unosc={} osc={}".format(fl,totflux[fl],totosc[fl])

    c.Update()
    raw_input()
    
    
