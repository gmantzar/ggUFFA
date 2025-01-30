import ROOT

import sys
sys.path.append('/Users/georgios/UFFA')

import FileReader as FR



class UffaReader(FR.FileReader):
    def __init__(self, ifile, thedirectory = None, therebin = False, thedifferential = None):
        """
        if directory == "":
            directory = "femto-dream-pair-task-track-track_std"
        elif directory[0] == '_':
            directory =="femto-dream-pair-task-track-track" + directory
        """

        self.rebin = therebin
        self.differential = thedifferential
        FR.FileReader.__init__(self, ifile, thedirectory)
        nevents = self.getNevents()



    ### Getter Functions for Uffa ###

    def get_zvtx(self):
        return self.get_histo("zvtxhist", "Event")

    def get_multNTr(self):
        return self.get_histo("MultNTr", "Event")

    def old_get_se(self):
        if self.rebin == False:
            return self.get_histo("SE")
        elif type(self.rebin) is int:
            return self.get_histo("SE", "rebin: "+str(self.rebin))
        else:
            print("Please enter either no rebin or an integer as rebin factor")
            return None

    def old_get_se_mc(self):
        if self.rebin == False:
            return self.get_histo("SE", "mc")
        elif type(self.rebin) is int:
            return self.get_histo("SE", "mc/rebin: "+str(self.rebin))
        else:
            print("Please enter either no rebin or an integer as rebin factor")
            return None
    
    def old_get_me(self):

        if self.rebin == False:
            return self.get_histo("ME")
        elif type(self.rebin) is int:
            return self.get_histo("ME", "rebin: "+str(self.rebin))
        else:
            print("Please enter either no rebin or an integer as rebin factor")
            return None

    def old_get_me_mc(self):
        if self.rebin == False:
            return self.get_histo("ME", "mc")
        elif type(self.rebin) is int:
            return self.get_histo("ME", "mc/rebin: "+str(self.rebin))
        else:
            print("Please enter either no rebin or an integer as rebin factor")
            return None

    def old_get_cf(self):
        if self.rebin == False:
            return self.get_histo("CF")
        elif type(self.rebin) is int:
            return self.get_histo("CF", "rebin: "+str(self.rebin))
        else:
            print("Please enter either no rebin or an integer as rebin factor")
            return None

    def old_get_cf_mc(self):
        if self.rebin == False:
            return self.get_histo("CF", "mc")
        elif type(self.rebin) is int:
            return self.get_histo("CF", "mc/rebin: "+str(self.rebin))
        else:
            print("Please enter either no rebin or an integer as rebin factor")
            return None




    def get_kstar_plot_base(self, name, mc=False):

        mcDir = ""
        rebinDir = ""
        if(mc):
            mcDir = "mc"
        if type(self.rebin) is int:
            rebinDir = "/rebin: "+str(self.rebin)

        if self.differential == None:

            return self.get_histo(name, mcDir+rebinDir)
        
        elif type(self.differential) == int:
            return self.get_histo(name, mcDir + "bin: "+str(self.differential) + rebinDir)
            
        elif type(self.differential) == list:
            hlist = []
            
            #directory = directory + "rebin: "+str(self.rebin)
            
            for i in range(len(self.differential)):
                hlist.append(self.get_histo(name, mcDir + "bin: " + str(i) + rebinDir))
            
            return hlist
            


    def get_se(self):
        return self.get_kstar_plot_base("SE")

    def get_se_mc(self):
        return self.get_kstar_plot_base("SE", True)
    
    def get_me(self):
        return self.get_kstar_plot_base("ME")

    def get_me_mc(self):
        return self.get_kstar_plot_base("ME", True)
    
    def get_cf(self):
        return self.get_kstar_plot_base("CF")

    def get_cf_mc(self):
        return self.get_kstar_plot_base("CF", True)



    
    
    def get_pt(self):
        return self.get_histo("hPt", "Tracks_one")

    def get_pt_mc(self):
        return self.get_histo("hPt", "Tracks_one_MC")

    def get_eta(self):
        return self.get_histo("hEta", "Tracks_one")

    def get_eta_mc(self):
        return self.get_histo("hEta", "Tracks_one_MC")

    def get_phi(self):
        return self.get_histo("hPhi", "Tracks_one")

    def get_phi_mc(self):
        return self.get_histo("hPhi", "Tracks_one_MC")

    def getNevents(self):
        zvtx = self.get_zvtx()
        return zvtx.GetEntries()


    ### Groupfetch histogramms ###

    def fetchCorrelationHistos(self):

        histos = []
        histos.append(self.get_se())
        histos.append(self.get_me())
        histos.append(self.get_cf())

        return histos

    def fetchQAHistos(self):

        histos = []
        histos.append(self.get_pt())
        histos.append(self.get_eta())
        histos.append(self.get_phi())

        return histos
