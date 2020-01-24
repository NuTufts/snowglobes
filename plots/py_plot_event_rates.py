from __future__ import print_function
import os,sys,argparse


def get_channel_data( chanfile ):
    """ collect channel information from files channels folder """
    
    channeldata = { "channame":[],
                    "channum":[],
                    "cpstate":[],
                    "flavor":[],
                    "num_target_factor":[],
                    "numchans":0}

    chanfilename = "../channels/channels_"+chanfile+".dat";
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
             "nux_es_smeared":np.zeros( maxpoints ) }

    # read in data
    for ifile in xrange(nchannels):

        # read unsmeared data        
        filename = "../out/"+fluxname+"_"+channeldata["channame"][ifile]+"_"+exptconfig+"_events.dat"
        f = open(filename,'r')
        ll = f.readlines()
        npoints = 0
        print("parsing ",filename)
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
        print(" found %d points for channel %d. integral=%f"%(npoints,ifile,np.sum( data["events"][ifile,:] ) ))
        f.close()

        # read smeared data
        filename_smeared = "../out/"+fluxname+"_"+channeldata["channame"][ifile]+"_"+exptconfig+"_events_smeared.dat";
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
            
        print(" found %d points for smeared channel %d; integral=%f"%(npoints_smeared,ifile, np.sum( data["events_smeared"][ifile,:] ) ))

    return data

if __name__ == "__main__":

    import numpy as np
    
    data = read_output_data( "stpi", "argon", "ar40kt" )
    print("total number of events: ",np.sum(data["total_events"])*3.14e7*0.000612/40.0*(20.0*20.0)/(29.0*29.0)*(5000.0/8766.0) )

#   // Plots for total events

#   TObjArray changraphlist(0);
#   TObjArray changraph_smearedlist(0);
#   TGraph* gr;
#   TGraph* gr_smeared;

#   gr= new TGraph(numpoints,en,total_events);
#   gr->SetLineWidth(3);
#   gr->SetLineColor(6);
#   TString graphname=fluxname+"_tot_"+exptconfig;
#   gr->SetName(graphname);

#   changraphlist.Add(gr);

#   gr_smeared= new TGraph(numpoints,en,total_events_smeared);
#   gr_smeared->SetLineWidth(3);
#   gr_smeared->SetLineColor(6);
#   gr_smeared->SetLineStyle(2);
#   graphname=fluxname+"_tot_"+exptconfig+"_smeared";
#   gr_smeared->SetName(graphname);

#   changraph_smearedlist.Add(gr_smeared);

#   gr_smeared= new TGraph(numpoints,en,nue_es_smeared);
#   gr_smeared->SetLineWidth(3);
#   gr_smeared->SetLineColor(6);
#   gr_smeared->SetLineStyle(2);
#   graphname=fluxname+"_nue_es_"+exptconfig+"_smeared";
#   gr_smeared->SetName(graphname);

#   changraph_smearedlist.Add(gr_smeared);

#   gr_smeared= new TGraph(numpoints,en,nuebar_es_smeared);
#   gr_smeared->SetLineWidth(3);
#   gr_smeared->SetLineColor(6);
#   gr_smeared->SetLineStyle(2);
#   graphname=fluxname+"_nuebar_es_"+exptconfig+"_smeared";
#   gr_smeared->SetName(graphname);

#   changraph_smearedlist.Add(gr_smeared);

#   gr_smeared= new TGraph(numpoints,en,nux_es_smeared);
#   gr_smeared->SetLineWidth(3);
#   gr_smeared->SetLineColor(6);
#   gr_smeared->SetLineStyle(2);
#   graphname=fluxname+"_nux_es_"+exptconfig+"_smeared";
#   gr_smeared->SetName(graphname);

#   changraph_smearedlist.Add(gr_smeared);

#   // Make the individual channel graphs

#   Int_t water_colors[maxchan]={1,2,2,2,2,2,2,3,4,7,7,7,7,7,7,0,
#                            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
#   Int_t argon_colors[maxchan]={2,2,2,2,2,2,3,4,7,7,7,7,7,7,0,
# 			       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
#   Int_t colors[maxchan];
#   Int_t ic=0;

#   if( exptconfig.CompareTo("ar17kt") == 0){
#     for (ic=0;ic<numchans;ic++) {
#       colors[ic] = argon_colors[ic];
#     } 
#   }
#   else {
#     // works for scint too
#     for (ic=0;ic<numchans;ic++) {
#       colors[ic] = water_colors[ic];
#     } 
#   }
 

#   Int_t m,n;

#   for (m=0;m<numchans;m++) {

#     gr= new TGraph(numpoints,en,events[m]);
#     gr->SetLineWidth(3); 
#     gr->SetLineColor(colors[m]);
#     graphname=fluxname+"_"+channame[channum[m]]+"_"+exptconfig;
#     gr->SetName(graphname);

#     changraphlist.Add(gr);

#     gr_smeared = new TGraph(numpoints,en,events_smeared[m]);
#     gr_smeared->SetLineWidth(3);
#     gr_smeared->SetLineColor(colors[m]);
#     gr_smeared->SetLineStyle(2);
#     graphname=fluxname+"_"+channame[channum[m]]+"_"+exptconfig+"_smeared";
#     gr_smeared->SetName(graphname);

#     changraph_smearedlist.Add(gr_smeared);

#   }

#   // Sum the NC and ES channels

#   Int_t ichan;
#   Double_t es_events[maxpoints],es_events_smeared[maxpoints];
#   Double_t ncO16_events[maxpoints],ncO16_events_smeared[maxpoints],ncC12_events[maxpoints],ncC12_events_smeared[maxpoints];

#   Double_t nc_nu_Pb_1n_events[maxpoints],nc_nu_Pb_2n_events[maxpoints],nc_nubar_Pb_1n_events[maxpoints],nc_nubar_Pb_2n_events[maxpoints];
#   Double_t tot_n_events[maxpoints];

#   for (j=0;j<numpoints;j++) {
#     es_events[j]=0;
#     es_events_smeared[j]=0;
#     ncO16_events[j]=0;
#     ncO16_events_smeared[j]=0;
#     ncC12_events[j]=0;
#     ncC12_events_smeared[j]=0;

#     tot_n_events[j]=0;
#     nc_nu_Pb_1n_events[j]=0;
#     nc_nu_Pb_2n_events[j]=0;
#     nc_nubar_Pb_1n_events[j]=0;
#     nc_nubar_Pb_2n_events[j]=0;

#     Int_t firsteschan, lasteschan;
#     if( exptconfig.CompareTo("ar17kt") == 0 || exptconfig.CompareTo("halo1")==0 || exptconfig.CompareTo("halo2") == 0){
#       firsteschan = 0;
#       lasteschan = 5;
#     } 
#     else {
#       firsteschan = 1;
#       lasteschan = 6;
#     }


#     for (ichan=firsteschan;ichan<=lasteschan;ichan++) {

#       es_events[j] += events[ichan][j];
#       es_events_smeared[j] += events_smeared[ichan][j];

#     }
#     if( exptconfig.CompareTo("wc100kt30prct") == 0 ||exptconfig.CompareTo("wc100kt15prct") == 0 ) {
#       for (ichan=9;ichan<15;ichan++) {
# 	ncO16_events[j] += events[ichan][j];
# 	ncO16_events_smeared[j] += events_smeared[ichan][j];
#       }
#       //	cout << "O16  " << ncO16_events_smeared[j] << endl;
#     }

#     if( exptconfig.CompareTo("scint50kt") == 0) {
#       for (ichan=9;ichan<15;ichan++) {
# 	ncC12_events[j] += events[ichan][j];
# 	ncC12_events_smeared[j] += events_smeared[ichan][j];
#       }

#     }

#     if( exptconfig.CompareTo("halo1") == 0 || exptconfig.CompareTo("halo2") == 0) {
#       nc_nu_Pb_1n_events[j] = events[8][j]+events[10][j]+events[12][j];
#       nc_nubar_Pb_1n_events[j] = events[9][j]+events[11][j]+events[13][j];
#       nc_nu_Pb_2n_events[j] = events[14][j]+events[16][j]+events[18][j];
#       nc_nubar_Pb_2n_events[j] = events[15][j]+events[17][j]+events[19][j];
#       tot_n_events[j] = nc_nu_Pb_1n_events[j]+nc_nubar_Pb_1n_events[j]+
# 	nc_nu_Pb_2n_events[j]+nc_nubar_Pb_2n_events[j]+
# 	events[6][j]+events[7][j];
#       //	cout << "C12  " << ncC12_events_smeared[j] << endl;
#     }
  

    
#   }


#   gr = new TGraph(numpoints,en,es_events);
#   gr->SetLineWidth(3); 
#   gr->SetLineColor(colors[1]);
#   graphname=fluxname+"_es_"+exptconfig;
#   gr->SetName(graphname);

#   changraphlist.Add(gr);
  
#   gr_smeared = new TGraph(numpoints,en,es_events_smeared);
#   gr_smeared->SetLineWidth(3); 
#   gr_smeared->SetLineColor(colors[1]);
#   gr_smeared->SetLineStyle(2); 
#   graphname=fluxname+"_es_"+exptconfig+"_smeared";
#   gr_smeared->SetName(graphname);

#   changraph_smearedlist.Add(gr_smeared);


#     // Oxygen

#   if( exptconfig.CompareTo("wc100kt30prct") == 0 ||exptconfig.CompareTo("wc100kt15prct") == 0 ) {
#     gr = new TGraph(numpoints,en,ncO16_events);
#     gr->SetLineWidth(3); 
#     gr->SetLineColor(colors[9]);
#     graphname=fluxname+"_ncO16_"+exptconfig;
#     gr->SetName(graphname);

#     changraphlist.Add(gr);

#     gr_smeared = new TGraph(numpoints,en,ncO16_events_smeared);
#     gr_smeared->SetLineWidth(3); 
#     gr_smeared->SetLineStyle(2); 
#     gr_smeared->SetLineColor(colors[9]);
#     graphname=fluxname+"_ncO16_"+exptconfig+"_smeared";
#     gr_smeared->SetName(graphname);

#     changraph_smearedlist.Add(gr_smeared);

#    }

#      // C12
#   if( exptconfig.CompareTo("scint50kt") == 0) {
#     gr = new TGraph(numpoints,en,ncC12_events);
#     gr->SetLineWidth(3); 
#     gr->SetLineColor(colors[9]);
#     graphname=fluxname+"_ncC12_"+exptconfig;
#     gr->SetName(graphname);

#     changraphlist.Add(gr);

#     gr_smeared = new TGraph(numpoints,en,ncC12_events_smeared);
#     gr_smeared->SetLineWidth(3); 
#     gr_smeared->SetLineStyle(2); 
#     gr_smeared->SetLineColor(colors[9]);
#     graphname=fluxname+"_ncC12_"+exptconfig+"_smeared";
#     gr_smeared->SetName(graphname);

#     changraph_smearedlist.Add(gr_smeared);
#   }

#   if( exptconfig.CompareTo("halo1") == 0 || exptconfig.CompareTo("halo2") == 0) {

#     gr = new TGraph(numpoints,en,tot_n_events);
#     gr->SetLineWidth(3); 
#     gr->SetLineColor(6);
#     graphname=fluxname+"_tot_n_"+exptconfig;
#     gr->SetName(graphname);

#     changraphlist.Add(gr);

#     gr = new TGraph(numpoints,en,nc_nu_Pb_1n_events);
#     gr->SetLineWidth(3); 
#     gr->SetLineColor(colors[9]);
#     graphname=fluxname+"_nc_nu_Pb_1n_"+exptconfig;
#     gr->SetName(graphname);

#     changraphlist.Add(gr);

#     gr = new TGraph(numpoints,en,nc_nubar_Pb_1n_events);
#     gr->SetLineWidth(3); 
#     gr->SetLineColor(colors[9]);
#     graphname=fluxname+"_nc_nubar_Pb_1n_"+exptconfig;
#     gr->SetName(graphname);

#     changraphlist.Add(gr);

#     gr = new TGraph(numpoints,en,nc_nu_Pb_2n_events);
#     gr->SetLineWidth(3); 
#     gr->SetLineColor(colors[9]);
#     graphname=fluxname+"_nc_nu_Pb_2n_"+exptconfig;
#     gr->SetName(graphname);

#     changraphlist.Add(gr);

#     gr = new TGraph(numpoints,en,nc_nubar_Pb_2n_events);
#     gr->SetLineWidth(3); 
#     gr->SetLineColor(colors[9]);
#     graphname=fluxname+"_nc_nubar_Pb_2n_"+exptconfig;
#     gr->SetName(graphname);

#     changraphlist.Add(gr);


#   }
#   // Output the graphs to a file

#   TString graphfilename = "graphs_"+fluxname+"_"+chanfile+"_"+exptconfig+".root";
#   TFile f(graphfilename,"recreate");
#   changraphlist.Write();
#   changraph_smearedlist.Write();
#   f.Close();

# //---------------------------------------------------------------------------------------------------------------------------------------------------------


# }
        

# def plot_event_rates( fluxname, chanfile, exptconfig, logopt=0, titleopt=0):
#     import ROOT as rt
#     rt.gStyle.SetOptStat(0)

#     # First read in the channels to plot

#     const Int_t maxchan = 32;
#     TString channame[maxchan];
#     Int_t channum[maxchan];
#     Int_t numchans;
#     TString cpstate[maxchan];
#     TString flav[maxchan];
#     Double_t num_target_factor[maxchan];


#     graphs = []
    

#     graphfilename = "graphs_"+fluxname+"_"+chanfile+"_"+exptconfig+".root";
#   TFile f(graphfilename);

#   TGraph* gr;
#   TGraph* gr_smeared;

#   // Get the maximum  
#   TString graphname=fluxname+"_tot_"+exptconfig+"_smeared";
#   gr = (TGraph*)f.Get(graphname);

#   Double_t x1,y1,x2,y2;
#   gr->ComputeRange(x1,y1,x2,y2);
#   //  maxy *=1.6;

#   Double_t xmin;
#   if (exptconfig.CompareTo("ar17kt") == 0 || exptconfig.CompareTo("ar17ktres1")==0|| exptconfig.CompareTo("ar17ktres2") == 0) {
#     xmin = 5.;
#   } else {
#     xmin = 5.;
#   }

#   Double_t xmax = 100.;
#   Double_t ymin = 0.1;
#   Double_t ymax;

#   //  cout << "Max y "<<y2<<endl;
#   if (logopt != 0) {
#     ymax = y2*1.6;
#   } else {
#     ymax = y2*1.1;
#   }
#   //  Double_t ymax = 1100.;

#   TCanvas* canv = new TCanvas("c1"," ",800,700);

#   TString plotexpt;
#   if (exptconfig == "100kt30pc") {
#     plotexpt = "100kt 30%";
#   } else {
#     if (exptconfig == "100kt15pc") {
#       plotexpt = "100kt 15%";
#     } else {
#       plotexpt = exptconfig;
#     }
#   }

#   TString plot_title;
#   if (titleopt == 0) {
#     plot_title = " "; 
#   } else {
#     plot_title = "Flux: "+fluxname+"      Detector: "+plotexpt;
#   }
#   TH2F *hr = new TH2F("hr",plot_title,2,xmin,xmax,2,ymin,ymax);
#   gStyle->SetTitleFontSize(.04);
#   hr->SetXTitle("Energy (MeV) ");
#   hr->SetYTitle("Events per 0.5 MeV");
#   hr->GetXaxis()->SetTitleSize(0.04);
#   hr->GetXaxis()->SetTitleOffset(1.3);
#   hr->GetXaxis()->SetLabelSize(0.04);
#   hr->GetYaxis()->SetTitleSize(0.04);
#   hr->GetYaxis()->SetTitleOffset(1.);
#   hr->GetYaxis()->SetLabelSize(0.04);

#   canv->SetLogy(logopt);

#   canv->cd();
#   hr->Draw();
#   leg = new TLegend(0.6,0.7,0.8,0.75);
#   leg->SetFillColor(10);
#   leg->SetTextFont((Font_t) 62);
#   leg->SetTextSize(0.04);
#   leg->SetBorderSize(0);

#   graphname=fluxname+"_tot_"+exptconfig;
#   gr = (TGraph*)f.Get(graphname);
#   gr->SetLineWidth(4);
#   gr->SetLineColor(6);

#   gr->Print("all");

#   graphname=fluxname+"_tot_"+exptconfig+"_smeared";
#   gr_smeared = (TGraph*)f.Get(graphname);
#   gr_smeared->SetLineWidth(4);
#   gr_smeared->SetLineColor(6);
#   gr_smeared->SetLineStyle(2);

#   gr_smeared->Print("all");

#   // Graph with error bars from the smeared graph-- doesn't seem to be a builtin method

#   const Int_t maxpoints=1000;
#   Int_t numgrpts = gr_smeared->GetN();
#   Double_t x[maxpoints],y[maxpoints];
#   Double_t ex[maxpoints],ey[maxpoints];
 
#   Int_t k;
#   for (k=0;k<numgrpts;k++) {
#     gr_smeared->GetPoint(k,x[k],y[k]);

#     if (y[k]>=0) {
#       ey[k] = sqrt(y[k]);
#     } else {
#       ey[k]=0.;
#     }

#     ex[k]=0.;

#   }

#  TGraphErrors* gr_smeared_errors = new TGraphErrors(numgrpts,x,y,ex,ey);

#   gr->Draw();
#   leg->AddEntry(gr,"Total","l");

#   gr_smeared->Draw("same");
#   //gr_smeared_errors->Draw();
#   leg->AddEntry(gr_smeared,"Total, smeared","l");

#   leg->Draw();

#   TString printfilename = "tot_rates_"+fluxname+"_"+exptconfig+".pdf";

#   canv->Print(printfilename);


# //-------------------------------------------------------------------- interaction rates -------------------------------------------------------------------

#   TCanvas* canv2 = new TCanvas("c2"," ",800,700);


#    canv2->cd();
#    hr->Draw();
#    canv2->SetLogy(logopt);
#    leg2 = new TLegend(0.6,0.6,0.88,0.85);
#    leg2->SetFillColor(10);
#    leg2->SetTextFont((Font_t) 62);
#    leg2->SetTextSize(0.04);
#    leg2->SetBorderSize(0);

   
#    graphname=fluxname+"_tot_"+exptconfig;
#    gr = (TGraph*)f.Get(graphname);
#    gr->SetLineWidth(3);
#    gr->SetLineColor(6);
#    gr->Draw("same");
#    leg2->AddEntry(gr,"Total","l");
   
#    graphname=fluxname+"_es_"+exptconfig;
#    gr = (TGraph*)f.Get(graphname);
#    gr->SetLineWidth(3);
#    gr->SetLineColor(2);
#    gr->Draw("same");
#    leg2->AddEntry(gr,"ES","l");
     
#       if( exptconfig.CompareTo("ar17kt") == 0 || exptconfig.CompareTo("ar17ktres1") == 0 || exptconfig.CompareTo("ar17ktres2")==0){

#      graphname=fluxname+"_nue_Ar40_"+exptconfig;
#      gr = (TGraph*)f.Get(graphname);
#      gr->SetLineWidth(3);
#      gr->SetLineColor(3);
#      gr->Draw("same");
#      leg2->AddEntry(gr,"#nu_{e} ^{40}Ar","l");

#      graphname=fluxname+"_nuebar_Ar40_"+exptconfig;
#      gr = (TGraph*)f.Get(graphname);
#      gr->SetLineWidth(3);
#      gr->SetLineColor(4);
#      gr->Draw("same");
#      leg2->AddEntry(gr,"#bar{#nu}_{e} ^{40}Ar","l");


#    } 
#    else  {

#      if( exptconfig.CompareTo("scint50kt") != 0) {

#      graphname=fluxname+"_ibd_"+exptconfig;
#      gr = (TGraph*)f.Get(graphname);
#      gr->SetLineWidth(3);
#      gr->SetLineColor(1);
#      gr->Draw("same");
#      leg2->AddEntry(gr,"IBD","l");

#      graphname=fluxname+"_nue_O16_"+exptconfig;
#      gr = (TGraph*)f.Get(graphname);
#      gr->SetLineWidth(3);
#      gr->SetLineColor(3);
#      gr->Draw("same");
#      leg2->AddEntry(gr,"#nu_{e}-^{16}O","l");

#      graphname=fluxname+"_nuebar_O16_"+exptconfig;
#      gr = (TGraph*)f.Get(graphname);
#      gr->SetLineWidth(3);
#      gr->SetLineColor(4);
#      gr->Draw("same");
#      leg2->AddEntry(gr,"#bar{#nu}_{e}-^{16}O","l");

#      graphname=fluxname+"_ncO16_"+exptconfig;
#      gr = (TGraph*)f.Get(graphname);
#           gr->SetLineWidth(3);
#      gr->SetLineColor(7);
#      gr->Draw("same");
#      leg2->AddEntry(gr,"NC ^{16}O","l");

#      }

#      else {

#      graphname=fluxname+"_ibd_"+exptconfig;
#      gr = (TGraph*)f.Get(graphname);
#      gr->SetLineWidth(3);
#      gr->SetLineColor(1);
#      gr->Draw("same");
#      leg2->AddEntry(gr,"IBD","l");

#      graphname=fluxname+"_nue_C12_"+exptconfig;
#      gr = (TGraph*)f.Get(graphname);
#      gr->SetLineWidth(3);
#      gr->SetLineColor(3);
#      gr->Draw("same");
#      leg2->AddEntry(gr,"#nu_{e}-^{12}C","l");

#      graphname=fluxname+"_nuebar_C12_"+exptconfig;
#      gr = (TGraph*)f.Get(graphname);
#      gr->SetLineWidth(3);
#      gr->SetLineColor(4);
#      gr->Draw("same");
#      leg2->AddEntry(gr,"#bar{#nu}_{e}-^{12}C","l");

#      graphname=fluxname+"_ncC12_"+exptconfig;
#      gr = (TGraph*)f.Get(graphname);
#      gr->SetLineWidth(3);
#      gr->SetLineColor(7);
#      gr->Draw("same");
#      leg2->AddEntry(gr,"NC ^{12}C","l");

#      }     

#    }

#    leg2->Draw();

#    printfilename = "interaction_rates_"+fluxname+"_"+exptconfig+".pdf";

#    canv2->Print(printfilename);


# //--------------------------------------------------------------- smeared rates channel contribution -------------------------------------------------------

#    TCanvas* canv3 = new TCanvas("c3"," ",800,700);

#    canv3->cd();
#    hr->Draw();
#    canv3->SetLogy(logopt);

#    TLegend* leg3 = new TLegend(0.6,0.6,0.88,0.85);
#    leg3->SetFillColor(10);
#    leg3->SetTextFont((Font_t) 62);
#    leg3->SetTextSize(0.04);
#    leg3->SetBorderSize(0);

   
#    graphname=fluxname+"_tot_"+exptconfig+"_smeared";
#    gr = (TGraph*)f.Get(graphname);
#    gr->SetLineWidth(3);
#    gr->SetLineColor(6);
#    gr->SetLineStyle(2);
#    gr->Draw("same");
#    leg3->AddEntry(gr,"Total","l");
   
#    graphname=fluxname+"_es_"+exptconfig+"_smeared";
#    gr = (TGraph*)f.Get(graphname);
#    gr->SetLineWidth(3);
#    gr->SetLineColor(2);
#    gr->SetLineStyle(2);
#    gr->Draw("same");
#    leg3->AddEntry(gr,"ES","l");
     
#       if( exptconfig.CompareTo("ar17kt") == 0 || exptconfig.CompareTo("ar17ktres1") == 0 || exptconfig.CompareTo("ar17ktres2")==0){

#      graphname=fluxname+"_nue_Ar40_"+exptconfig+"_smeared";
#      gr = (TGraph*)f.Get(graphname);
#      gr->SetLineStyle(2);
#      gr->SetLineWidth(3);
#      gr->SetLineColor(3);
#      gr->Draw("same");
#      leg3->AddEntry(gr,"#nu_{e} ^{40}Ar","l");

#      graphname=fluxname+"_nuebar_Ar40_"+exptconfig+"_smeared";
#      gr = (TGraph*)f.Get(graphname);
#      gr->SetLineStyle(2);
#      gr->SetLineWidth(3);
#      gr->SetLineColor(4);
#      gr->Draw("same");
#      leg3->AddEntry(gr,"#bar{#nu}_{e} ^{40}Ar","l");


#    }

#    if( exptconfig.CompareTo("wc100kt30prct") == 0 || exptconfig.CompareTo("wc100kt15prct") == 0) {

#      graphname=fluxname+"_ibd_"+exptconfig+"_smeared";
#      gr = (TGraph*)f.Get(graphname);
#      gr->SetLineStyle(2);
#      gr->SetLineWidth(3);
#      gr->SetLineColor(1);
#      gr->Draw("same");
#      leg3->AddEntry(gr,"IBD","l");

#      graphname=fluxname+"_nue_O16_"+exptconfig+"_smeared";
#      gr = (TGraph*)f.Get(graphname);
#      gr->SetLineStyle(2);
#      gr->SetLineWidth(3);
#      gr->SetLineColor(3);
#      gr->Draw("same");
#      leg3->AddEntry(gr,"#nu_{e}-^{16}O","l");

#      graphname=fluxname+"_nuebar_O16_"+exptconfig+"_smeared";
#      gr = (TGraph*)f.Get(graphname);
#      gr->SetLineStyle(2);
#      gr->SetLineWidth(3);
#      gr->SetLineColor(4);
#      gr->Draw("same");
#      leg3->AddEntry(gr,"#bar{#nu}_{e}-^{16}O","l");

#      graphname=fluxname+"_ncO16_"+exptconfig+"_smeared";
#      gr = (TGraph*)f.Get(graphname);
#      gr->SetLineStyle(2);
#      gr->SetLineWidth(3);
#      gr->SetLineColor(7);
#      gr->Draw("same");
#      leg3->AddEntry(gr,"NC ^{16}O","l");
#    }

#    if( exptconfig.CompareTo("scint50kt") == 0) {

#      graphname=fluxname+"_ibd_"+exptconfig+"_smeared";
#      gr = (TGraph*)f.Get(graphname);
#      gr->SetLineStyle(2);
#      gr->SetLineWidth(3);
#      gr->SetLineColor(1);
#      gr->Draw("same");
#      leg3->AddEntry(gr,"IBD","l");

#      graphname=fluxname+"_nue_C12_"+exptconfig+"_smeared";
#      gr = (TGraph*)f.Get(graphname);
#      gr->SetLineStyle(2);
#      gr->SetLineWidth(3);
#      gr->SetLineColor(3);
#      gr->Draw("same");
#      leg3->AddEntry(gr,"#nu_{e}-^{12}C","l");

#      graphname=fluxname+"_nuebar_C12_"+exptconfig+"_smeared";
#      gr = (TGraph*)f.Get(graphname);  
#      gr->SetLineStyle(2);
#      gr->SetLineWidth(3);
#      gr->SetLineColor(4);
#      gr->Draw("same");
#      leg3->AddEntry(gr,"#bar{#nu}_{e}-^{12}C","l");

#      graphname=fluxname+"_ncC12_"+exptconfig+"_smeared";
#      gr = (TGraph*)f.Get(graphname);
#      gr->SetLineStyle(2);
#      gr->SetLineWidth(3);
#      gr->SetLineColor(7);
#      gr->Draw("same");
#      leg3->AddEntry(gr,"NC ^{12}C","l");
#    }


#    leg3->Draw();

#    printfilename = "smeared_rates_"+fluxname+"_"+exptconfig+".pdf";

#    canv3->Print(printfilename);

# //-------------------------------------------------------------------true flavor content ES-----------------------------------------------------------------

#   TCanvas* canv4 = new TCanvas("c4"," ",800,700);

#    canv4->cd();
#    hr->Draw();
#    canv4->SetLogy(logopt);

#    leg4 = new TLegend(0.6,0.6,0.88,0.85);
#    leg4->SetFillColor(10);
#    leg4->SetTextFont((Font_t) 62);
#    leg4->SetTextSize(0.04);
#    leg4->SetBorderSize(0);

   
#    graphname=fluxname+"_es_"+exptconfig+"_smeared";
#    gr = (TGraph*)f.Get(graphname);
#    gr->SetLineWidth(3);
#    gr->SetLineColor(2);
#    gr->SetLineStyle(2);
#    gr->Draw("same");
#    leg4->AddEntry(gr,"ES total","l");
   
#    graphname=fluxname+"_nue_es_"+exptconfig+"_smeared";
#    gr = (TGraph*)f.Get(graphname);
#    gr->SetLineWidth(3);
#    gr->SetLineColor(3);
#    gr->SetLineStyle(2);
#    gr->Draw("same");
#    leg4->AddEntry(gr,"#nu_{e}","l");

#    graphname=fluxname+"_nuebar_es_"+exptconfig+"_smeared";
#    gr = (TGraph*)f.Get(graphname);
#    gr->SetLineWidth(3);
#    gr->SetLineColor(4);
#    gr->SetLineStyle(2);
#    gr->Draw("same");
#    leg4->AddEntry(gr,"#bar{#nu}_{e}","l");

#    graphname=fluxname+"_nux_es_"+exptconfig+"_smeared";
#    gr = (TGraph*)f.Get(graphname);
#    gr->SetLineWidth(3);
#    gr->SetLineColor(7);
#    gr->SetLineStyle(2);
#    gr->Draw("same");
#    leg4->AddEntry(gr,"#nu_{x}","l");

#    leg4->Draw();

#    printfilename = "true_flavor_content_ES_"+fluxname+"_"+exptconfig+".pdf";

#    canv4->Print(printfilename);

# }
