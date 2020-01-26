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
    s2ee = params["sin2_{ee}"]
    s2mm = params["sin2_{mm}"]
    s2me = params["sin2_{me}"]
    dm2  = params["dm2"]
    L    = params["L_m"]*0.001

    oscflux = np.zeros( flux.shape )
    sin2dm_arg = 1.27*dm2*L/flux[:,0]
    sin2dm = np.sin( sin2dm_arg )
    sin2dm = np.power( sin2dm, 2 )
    prob_ee = 1.0-s2ee*sin2dm
    prob_mm = 1.0-s2mm*sin2dm
    prob_me = s2me*sin2dm

    #print "arg: ",sin2dm_arg[:10]
    #print "E: ",flux[:10,0]
    #print "sin2dm: ",sin2dm[:10]

    # copy energies
    oscflux[:,0] = flux[:,0]
    
    # nue
    oscflux[:,1] = flux[:,1]*prob_ee[:] + flux[:,2]*prob_me
    oscflux[:,2] = flux[:,2]*prob_mm[:]
    oscflux[:,4] = flux[:,4]*prob_ee[:] + flux[:,4]*prob_me
    oscflux[:,5] = flux[:,5]*prob_mm[:]

    print "osc applied for {} in {} secs".format( params, time.time()-s )
        
    return oscflux

def make_flux_tgraph( flux ):
    npts = flux.shape[0]
    tgraph_v = {"nue":rt.TGraph(npts),
                "numu":rt.TGraph(npts),
                "nuebar":rt.TGraph(npts),
                "numubar":rt.TGraph(npts) }
    for i in xrange(npts):
        tgraph_v["nue"].SetPoint(i, flux[i,0],flux[i,1])
        tgraph_v["numu"].SetPoint(i,flux[i,0],flux[i,2])
        tgraph_v["nuebar"].SetPoint(i, flux[i,0],flux[i,4])
        tgraph_v["numubar"].SetPoint(i,flux[i,0],flux[i,5])

    return tgraph_v
        

if __name__ == "__main__":

    flux = read_flux( "stpi.dat" )

    params = { "sin2_{ee}":0.01,
               "sin2_{mm}":0.01,
               "sin2_{me}":0.01,
               "dm2":7.0,
               "L_m":29.0 }

    oscflux = oscillate_flux_sterile( flux, params )

    gflux = make_flux_tgraph(flux)
    oflux = make_flux_tgraph(oscflux)

    print oscflux[0:10,:]

    c = rt.TCanvas("c","c",1200,600)
    c.Draw()

    for g in [gflux,oflux]:
        g["numu"].SetLineColor(rt.kRed)
        g["nue"].SetLineColor(rt.kBlue)
        g["numubar"].SetLineColor(rt.kMagenta)
        g["nuebar"].SetLineColor(rt.kCyan)

    for n,g in oflux.items():
        g.SetLineStyle(2)

    gflux["numu"].Draw("AL")
    for n,g in gflux.items():
        g.Draw("Lsame")
    for n,g in oflux.items():
        g.Draw("Lsame")
    gflux["numu"].GetXaxis().SetRangeUser(0,0.060)

    c.Update()
    raw_input()
    
    
