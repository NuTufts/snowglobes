import os,sys
from math import sqrt

from oscillate_flux_sterile import read_flux, oscillate_flux_sterile, make_flux_file

if __name__ == "__main__":

    paramfile = sys.argv[1]    
    jobid     = int(sys.argv[2])
    input_fluxfile  = sys.argv[3]
    data_dir        = sys.argv[4]
    outdir          = sys.argv[5]

    channame = "argon_marley1"
    expt_config = "ar40kt"
    L = 29.0

    # get parameters
    fpar = open(paramfile,'r')
    lpar = fpar.readlines()[jobid+1]
    lpar = lpar.strip().split()
    print lpar
    dm2   = float( lpar[0] )
    ue4sq = float( lpar[1] )
    um4sq = float( lpar[2] )
    ut4sq = float( lpar[3] )

    oscpars = {"dm2":dm2,
               "Ue4":sqrt(ue4sq),
               "Um4":sqrt(um4sq),
               "Ut4":sqrt(ut4sq),
               "L_m":L}

    sin2_mue = 4.0*ue4sq*um4sq
    sin2_ee  = 4.0*(1-ue4sq)*ue4sq
    sin2_mm  = 4.0*(1-um4sq)*um4sq
    print "[run osc: sin2_{mue}]",sin2_mue
    print "[run osc: sin2_{ee}]",sin2_ee
    print "[run osc: sin2_{mm}]",sin2_mm

    flux_unosc = read_flux( input_fluxfile )    
    oscflux = oscillate_flux_sterile( flux_unosc, oscpars )
    make_flux_file( data_dir+"/fluxes/osc_stpi.dat", oscflux )

    oscflux_noappear = oscillate_flux_sterile( flux_unosc, oscpars, no_appear=True )
    make_flux_file( data_dir+"/fluxes/osc_stpi_noappear.dat", oscflux_noappear )

    oscflux_nodisappear = oscillate_flux_sterile( flux_unosc, oscpars, no_disappear=True )
    make_flux_file( data_dir+"/fluxes/osc_stpi_nodisappear.dat", oscflux_nodisappear )
    
    os.system("./supernova.pl osc_stpi %s %s 0 %s"%(channame,expt_config,data_dir))
    os.system("./supernova.pl osc_stpi_noappear %s %s 0 %s"%(channame,expt_config,data_dir))
    os.system("./supernova.pl osc_stpi_nodisappear %s %s 0 %s"%(channame,expt_config,data_dir))
    
    
