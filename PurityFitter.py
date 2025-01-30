from enum import Enum
import ROOT

def pidBkgPrefitFunction(x, par):
    lower = par[4] 
    upper = par[5] 
    
    if (x[0] > lower) and (x[0] < upper):
        ROOT.TF1.RejectPoint()
        return 0
    return par[0] + par[1]*x[0] + par[2]*ROOT.TMath.Exp(-par[3]*x[0])

def invMassLambdaBkgPrefitFunction(x, par):
    lower = par[0] 
    upper = par[1] 
    
    if (x[0] > lower) and (x[0] < upper):
        ROOT.TF1.RejectPoint()
        return 0
    return par[2] + par[3]*x[0] + par[4]*x[0]*x[0] + par[5]*x[0]*x[0]*x[0]



class PurityFitter:

    class Particles(Enum):
        proton = 1
        Lambda = 2
        Xi = 3

    class Histograms(Enum):
        tpcpid = 1
        tpcpidtail = 2
        tofpid = 3
        combpid_simple = 4
        combpid_gausBkg = 5
        combpid = 6
        invMass = 7

    def __init__(self, part, hist, fitvarVSpt, ptConfig, fittername='Fitter'):

        self.fitvarVSpt = fitvarVSpt
        self.ptConfig = ptConfig
        self.fittername = fittername
        
        self.particleType = None
        self.histType = None
        self.signalEdges = []
        self.bkgEdges = []
        self.bkgInnerEdges = []

        self.fitSlices = []
        self.set_fitSlices()
        self.NfitSlices = len(self.fitSlices)

        self.set_particleType(part) #for example "proton" or "Lambda"
        self.set_histType(hist)     # "tpc", "tof", "invMass"

        self.fCombined = []
        self.fSignal = []
        self.fSignalPrefit = []
        self.fBkg = []
        self.fBkgPrefit = []
        self.fBkgLeft = []
        self.fBkgRight = []

        self.swapBkgOrder = False #used to swap the order of the two background functions

        self.purities = []
        self.chi2_comb = []
        self.ndf_comb = []
        self.chi2_sig = []
        self.ndf_sig = []
        
        #configurations of the fit
        self.discretizeResults = False
        self.hCombined = []
        self.hSignal = []
        
        self.useLimits = False
        self.limitValues = []

        self.fitResult = []

        self.set_fitFunctions()



        #variables to configure the fit settings, etc.
        self.prefitSeparatedBkg = False
        self.rejectSignal = False 

    def set_particleType(self, partType):
        if partType == "p":
            self.particleType =  self.Particles.proton.value
        elif partType == "Lambda":
            self.particleType = self.Particles.Lambda.value
        elif partType == "Xi":
            self.particleType = self.Particles.Xi.value
        else:
            print('Please select either "Proton" or "Lambda"')

    def set_histType(self, histType):

        if histType == "tpc":
            self.histType = self.Histograms.tpcpid.value
        
        elif histType == "tpctail":
            self.histType = self.Histograms.tpcpidtail.value

        elif histType == "tof":
            self.histType = self.Histograms.tofpid.value

        elif histType == "comb_simple":
            self.histType = self.Histograms.combpid_simple.value
        
        elif histType == "comb_gausBkg":
            self.histType = self.Histograms.combpid_gausBkg.value
        
        elif histType == "comb":
            self.histType = self.Histograms.combpid.value

        elif histType == "invMass":
            self.histType = self.Histograms.invMass.value
        else:
            print("Histogram Type is not selected!")

    def set_fitSlices(self):

        if self.ptConfig[0] == "allBins":
        # make the projections for all bins starting at bin number self.ptConfig[1] up to self.ptConfig[2]

            for i in range(self.ptConfig[1], self.ptConfig[2]+1):
                self.fitSlices.append( self.fitvarVSpt.ProjectionY( self.fittername+"_slice_"+str(i)+"_"+str(i), i, i) )

        elif self.ptConfig[0] == "binnumber":
        # make the projections for each pair of bin edges, as specified in the self.ptConfig

            for i in range(1, len(self.ptConfig)):
                self.fitSlices.append( self.fitvarVSpt.ProjectionY( self.fittername+"_slice_"+str(self.ptConfig[i][0])+"_"+str(self.ptConfig[i][1]), self.ptConfig[i][0], self.ptConfig[i][1]) )
        else:
            for i in range(len(self.ptConfig)-1):
                self.fitSlices.append( self.fitvarVSpt.ProjectionY( self.fittername+"_slice_"+str(self.ptConfig[i])+"_"+str(self.ptConfig[i+1]),\
                                                    self.fitvarVSpt.GetXaxis().FindBin( self.ptConfig[i] ),\
                                                    self.fitvarVSpt.GetXaxis().FindBin( self.ptConfig[i+1] )) )
                
            #self.fitSlices[i].SetName("slice_"+str(self.ptConfig[i])+"_"+str(self.ptConfig[i+1]))





    #TODO: write functions then enable the manuall setting of the fit function/ bkg functoin and signal function

    
    ###############################
    #Define here the functions that shall be available

    def set_tpcFitFnc(self):

        self.signalEdges = []
        self.bkgInnerEdges = []
        self.bkgEdges = []
        
        self.fCombined = []
        self.fSignal = []
        self.fBkg = []

        for i in range(len(self.fitSlices)):
            
            # set the default Ranges:
            self.bkgEdges.append([-10., 10.])
            self.bkgInnerEdges.append([-3.5, 3.5])
            self.signalEdges.append([-3., 3.])
            
            # set the fit functions and the corresponding QA functions
            self.fCombined.append( ROOT.TF1("fCombined_"+str(i),
                                            "[0] * TMath::Gaus(x, [1], [2]) + [3] + [4]*x + [5]*exp(-[6]*x)",
                                            self.bkgEdges[i][0],
                                            self.bkgEdges[i][1],
                                            7)
                                 )
                                
            self.fSignal.append( ROOT.TF1("fSignal_"+str(i),
                                          "[0] * TMath::Gaus(x, [1], [2])",
                                          self.signalEdges[i][0],
                                          self.signalEdges[i][1],
                                          3)
                               )
            self.fBkg.append( ROOT.TF1("fBkg_"+str(i),
                                        "[0] + [1]*x + [2]*exp(-[3]*x)",
                                        self.bkgEdges[i][0],
                                        self.bkgEdges[i][1],
                                        4)
                            )
 
    def set_tpcFitFncTail(self):

        self.signalEdges = []
        self.bkgInnerEdges = []
        self.bkgEdges = []
        
        self.fCombined = []
        self.fSignal = []
        self.fBkg = []
        self.fBkgPrefit = []
        self.fBkgLeft = []
        self.fBkgRight = []

        for i in range(len(self.fitSlices)):
            
            # set the default Ranges:
            self.bkgEdges.append([-10., 10.])
            self.bkgInnerEdges.append([-3.5, 3.5])
            self.signalEdges.append([-3., 3.])
            
            # set the fit functions and the corresponding QA functions
            self.fCombined.append( ROOT.TF1(self.fittername+"_fCombined_"+str(i),
                                            "(x <= ([3] + [1])) * [0] * TMath::Gaus(x, [1], [2]) + (x > ([3] + [1])) * [0] * TMath::Gaus([3] + [1], [1], [2]) * exp(-([3]) * (x - [3] - [1]) / ([2] * [2])) + [4] + [5]*x + [6]*exp(-[7]*x)",
                                            self.bkgEdges[i][0],
                                            self.bkgEdges[i][1],
                                            8)
                                 )
                                
            self.fSignal.append( ROOT.TF1(self.fittername+"_fSignal_"+str(i),
                                          "(x <= ([3] + [1])) * [0] * TMath::Gaus(x, [1], [2]) + (x > ([3] + [1])) * [0] * TMath::Gaus([3] + [1], [1], [2]) * exp(-([3]) * (x - [3] - [1]) / ([2] * [2]))",
                                          self.signalEdges[i][0],
                                          self.signalEdges[i][1],
                                          4)
                               )
            self.fBkg.append( ROOT.TF1(self.fittername+"_fBkg_"+str(i),
                                        "[0] + [1]*x + [2]*exp(-[3]*x)",
                                        self.bkgEdges[i][0],
                                        self.bkgEdges[i][1],
                                        4)
                            )

            self.fBkgPrefit.append( ROOT.TF1(self.fittername+"_fBkgPrefit_"+str(i),
                                    pidBkgPrefitFunction,
                                    self.bkgEdges[i][0],
                                    self.bkgEdges[i][1],
                                    6)
                                  )


            self.fBkgLeft.append( ROOT.TF1(self.fittername+"_fBkgLeft_"+str(i),
                                          "[0]*exp(-[1]*x)",
                                          self.bkgEdges[i][0],
                                          self.bkgInnerEdges[i][0],
                                          2)
                                )

            self.fBkgRight.append( ROOT.TF1(self.fittername+"_fBkgRight_"+str(i),
                                          "[0] + [1]*x",
                                          self.bkgInnerEdges[i][1],
                                          self.bkgEdges[i][1],
                                          2)
                                )
            
    def set_tofFitFncTail(self):

        self.bkgInnerEdges = []
        self.bkgEdges = []
        self.signalEdges = []

        self.fCombined = []
        self.fSignal = []
        self.fBkg = []
        self.fBkgPrefit = []
        self.fBkgLeft = []
        self.fBkgRight = []

        for i in range(len(self.fitSlices)):

            # set the default Ranges:
            self.bkgEdges.append([-10., 10.])
            self.bkgInnerEdges.append([-3.5, 3.5])
            self.signalEdges.append([-3., 3.])

            # set the fit functions and the corresponding QA functions
            self.fCombined.append( ROOT.TF1("fCombined_"+str(i),
                                            "(x <= ([3] + [1])) * [0] * TMath::Gaus(x, [1], [2]) + (x > ([3] + [1])) * [0] * TMath::Gaus([3] + [1], [1], [2]) * exp(-([3]) * (x - [3] - [1]) / ([2] * [2])) + [4] + [5]*x + [6]*exp(-[7]*x)",
                                            self.bkgEdges[i][0],
                                            self.bkgEdges[i][1],
                                            8)
                                 )

            self.fSignal.append( ROOT.TF1("fSignal_"+str(i),
                                          "(x <= ([3] + [1])) * [0] * TMath::Gaus(x, [1], [2]) + (x > ([3] + [1])) * [0] * TMath::Gaus([3] + [1], [1], [2]) * exp(-([3]) * (x - [3] - [1]) / ([2] * [2]))",
                                          self.signalEdges[i][0],
                                          self.signalEdges[i][1],
                                          4)
                               )
            
            self.fBkg.append( ROOT.TF1("fBkg_"+str(i),
                                        "[0] + [1]*x + [2]*exp(-[3]*x)",
                                        self.bkgEdges[i][0],
                                        self.bkgEdges[i][1],
                                        4)
                            )
            
            self.fBkgPrefit.append( ROOT.TF1("fBkgPrefit_"+str(i),
                                    pidBkgPrefitFunction,
                                    self.bkgEdges[i][0],
                                    self.bkgEdges[i][1],
                                    4)
                                  )

            self.fBkgLeft.append( ROOT.TF1("fBkgLeft_"+str(i),
                                          "[0]*exp(-[3]*x)",
                                          self.bkgInnerEdges[i][1],
                                          self.bkgEdges[i][1],
                                          2)
                                )

            self.fBkgRight.append( ROOT.TF1("fBkgRight_"+str(i),
                                          "[0] + [1]*x",
                                          self.bkgEdges[i][0],
                                          self.bkgInnerEdges[i][0],
                                          2)
                                )
    
    def set_tofFitFncTailCustom(self):
    
        self.bkgInnerEdges = []
        self.bkgEdges = []
        self.signalEdges = []

        self.fCombined = []
        self.fSignal = []
        self.fBkg = []
        self.fBkgPrefit = []
        self.fBkgLeft = []
        self.fBkgRight = []

        for i in range(len(self.fitSlices)):

            # set the default Ranges:
            self.bkgEdges.append([-10., 10.])
            self.bkgInnerEdges.append([-3.5, 3.5])
            self.signalEdges.append([-3., 3.])

            # set the fit functions and the corresponding QA functions
            # "(x <= ([3] + [1])) * [0] * TMath::Gaus(x, [1], [2]) + (x > ([3] + [1])) * [0] * TMath::Gaus([3] + [1], [1], [2]) * exp(-([3]) * (x - [3] - [1]) / ([2] * [2])) + [4] + [5]*x + [6]*exp(-[7]*x)",
            self.fCombined.append( ROOT.TF1("fCombined_"+str(i),
                                            "[0] * exp(-0.5 * ((x - [1]) / [2])**2) + [3] * exp([4] * (x - [1])) + [5] * exp(-[6] * (x - [1])) + [7]*x + [8]*exp(-[9]*x)",
                                            self.bkgEdges[i][0],
                                            self.bkgEdges[i][1],
                                            8)
                                 )

            # "(x <= ([3] + [1])) * [0] * TMath::Gaus(x, [1], [2]) + (x > ([3] + [1])) * [0] * TMath::Gaus([3] + [1], [1], [2]) * exp(-([3]) * (x - [3] - [1]) / ([2] * [2]))",
            # "[0] * exp(-0.5 * ((x - [1]) / [2])**2) + [3] * exp([4] * (x - [1])) + [5] * exp(-[6] * (x - [1]))",
            self.fSignal.append( ROOT.TF1("fSignal_"+str(i),
                                          "(x <= ([4] - [1])) * [0] * TMath::Gaus([4] - [1], [1], [2]) * exp(([4]) * (x - [4] + [1]) / ([2] * [2]))  +  (([4] - [1]) < x <= ([3] + [1])) * [0] * TMath::Gaus(x, [1], [2]) + (x > ([3] + [1])) * [0] * TMath::Gaus([3] + [1], [1], [2]) * exp(-([3]) * (x - [3] - [1]) / ([2] * [2]))",
                                          self.signalEdges[i][0],
                                          self.signalEdges[i][1],
                                          5)
                               )
            
            self.fBkg.append( ROOT.TF1("fBkg_"+str(i),
                                        "[0] + [1]*x + [2]*exp(-[3]*x)",
                                        self.bkgEdges[i][0],
                                        self.bkgEdges[i][1],
                                        4)
                            )
            
            self.fBkgPrefit.append( ROOT.TF1("fBkgPrefit_"+str(i),
                                    pidBkgPrefitFunction,
                                    self.bkgEdges[i][0],
                                    self.bkgEdges[i][1],
                                    6)
                                  )

            self.fBkgLeft.append( ROOT.TF1("fBkgLeft_"+str(i),
                                          "[0]*exp(-[3]*x)",
                                          self.bkgInnerEdges[i][1],
                                          self.bkgEdges[i][1],
                                          2)
                                )

            self.fBkgRight.append( ROOT.TF1("fBkgRight_"+str(i),
                                          "[0] + [1]*x",
                                          self.bkgEdges[i][0],
                                          self.bkgInnerEdges[i][0],
                                          2)
                                )
    
    def set_combFitFnc_simple(self):
        
        self.bkgInnerEdges = []
        self.bkgEdges = []
        self.signalEdges = []

        self.fCombined = []
        self.fSignal = []
        self.fBkg = []
        
        self.fBkgLeft = None 
        self.fBkgRight = None 
           
        for i in range(len(self.fitSlices)):

            # set the default Ranges:
            self.bkgEdges.append([0., 15.])
            self.bkgInnerEdges.append([0., 0.])
            self.signalEdges.append([0., 5.])

            self.fCombined.append( ROOT.TF1("fCombined_"+str(i),
                                            #"[0] * TMath::Gaus(x, [1], [2]) + [3] +  [4] * exp([5]*x)",
                                            #"[0] * x/([1]*[1]) * exp( -x*x/(2.*[1]*[1]) ) + [2] +  [3] * exp([4]*x)",
                                            "(x <= [2]) * [0] * x/([1]*[1]) * exp( -x*x/(2.*[1]*[1]) ) + (x>[2]) * [0] * x/([1]*[1]) * exp( -([2]*[2])/(2.*[1]*[1]) ) * exp(-[2] * (x-[2])/([1]*[1])) + [3] + [4]*exp([5]*x) + [6]*x",
                                            self.bkgEdges[i][0],
                                            self.bkgEdges[i][1],
                                            7)
                                 )
            self.fSignal.append( ROOT.TF1("fSignal_"+str(i),
                                          #"[0] * TMath::Gaus(x, [1], [2])",
                                          #"[0] * x/([1]*[1]) * exp( -x*x/(2.*[1]*[1]) )",
                                          "(x <= [2]) * [0] * x/([1]*[1]) * exp( -x*x/(2.*[1]*[1]) ) + (x>[2]) * [0] * x/([1]*[1]) * exp( -([2]*[2])/(2.*[1]*[1]) ) * exp(-[2] * (x-[2])/([1]*[1]))",
                                          self.signalEdges[i][0],
                                          self.signalEdges[i][1],
                                          3)
                               )

            self.fBkg.append( ROOT.TF1("fBkg_"+str(i),
                                        "[0] + [1] * exp([2]*x) + [3]*x",
                                        self.bkgInnerEdges[i][0],
                                        self.bkgInnerEdges[i][1],
                                        4)
                            )
    
    def set_combFitFnc_gausBkg(self):
        """
        function to be used for the low pT bins for fitting the signal distribution. This is in order to account for (unknown) structures in the PID signal
        usage for the background:
        left Bkg: additional gauss used to describe the shoulder
        right Bkg: exponential background

        usage of the fit ranges:
        total fit range:
        left background: upper signal - lower bkg
        right bkg" lower bkg - upper bkg
        """ 
        print("Doing your thins!!")
        self.bkgInnerEdges = []
        self.bkgEdges = []
        self.signalEdges = []

        self.fCombined = []
        self.fSignal = []
        self.fBkg = []
        
        self.fBkgLeft = [] 
        self.fBkgRight = [] 
           
        for i in range(len(self.fitSlices)):

            print(i)
            # set the default Ranges:
            self.bkgEdges.append([0., 15.])
            self.bkgInnerEdges.append([0., 0.])
            self.signalEdges.append([0., 5.])

            self.fCombined.append( ROOT.TF1("fCombined_"+str(i),
                                            #"[0] * TMath::Gaus(x, [1], [2]) + [3] +  [4] * exp([5]*x)",
                                            #"[0] * x/([1]*[1]) * exp( -x*x/(2.*[1]*[1]) ) + [2] +  [3] * exp([4]*x)",
                                            "(x <= [2]) * [0] * x/([1]*[1]) * exp( -x*x/(2.*[1]*[1]) ) + (x>[2]) * [0] * x/([1]*[1]) * exp( -([2]*[2])/(2.*[1]*[1]) ) * exp(-[2] * (x-[2])/([1]*[1])) + [3]*exp([4]*x) + [5] * TMath::Gaus(x, [6], [7])",
                                            self.bkgEdges[i][0],
                                            self.bkgEdges[i][1],
                                            8)
                                 )
            self.fSignal.append( ROOT.TF1("fSignal_"+str(i),
                                          #"[0] * TMath::Gaus(x, [1], [2])",
                                          #"[0] * x/([1]*[1]) * exp( -x*x/(2.*[1]*[1]) )",
                                          "(x <= [2]) * [0] * x/([1]*[1]) * exp( -x*x/(2.*[1]*[1]) ) + (x>[2]) * [0] * x/([1]*[1]) * exp( -([2]*[2])/(2.*[1]*[1]) ) * exp(-[2] * (x-[2])/([1]*[1]))",
                                          self.signalEdges[i][0],
                                          self.signalEdges[i][1],
                                          3)
                               )

            self.fBkg.append( ROOT.TF1("fBkg_"+str(i),
                                        "[0]*exp([1]*x) + [2] * TMath::Gaus(x, [3], [4])",
                                        self.bkgEdges[i][0],
                                        self.bkgEdges[i][1],
                                        5)
                            )
            lowerlimit_leftBkg = 0.
            upperlimit_leftBkg = 0.
            lowerlimit_rightBkg = 0.
            upperlimit_rightBkg = 0.

            if self.swapBkgOrder is False:
                lowerlimit_leftBkg = self.signalEdges[i][1]
                upperlimit_leftBkg = self.bkgInnerEdges[i][0]
                lowerlimit_rightBkg = self.bkgInnerEdges[i][1]
                upperlimit_rightBkg = self.bkgEdges[i][1]
            else:
                lowerlimit_rightBkg = self.signalEdges[i][1]
                upperlimit_rightBkg = self.bkgInnerEdges[i][0]
                lowerlimit_leftBkg = self.bkgInnerEdges[i][1]
                upperlimit_leftBkg = self.bkgEdges[i][1]

            self.fBkgLeft.append( ROOT.TF1("fBkg_left_"+str(i),
                                           "[0] * TMath::Gaus(x, [1], [2])",
                                           lowerlimit_leftBkg,
                                           upperlimit_leftBkg,
                                           3)
                                )
            
            self.fBkgRight.append( ROOT.TF1("fBkg_left_"+str(i),
                                           "[0]*exp([1]*x)",
                                           lowerlimit_rightBkg,
                                           upperlimit_rightBkg,
                                           2)
                                )

    def set_combFitFnc(self):
        
        self.bkgInnerEdges = []
        self.bkgEdges = []
        self.signalEdges = []

        self.fCombined = []
        self.fSignal = []
        self.fBkg = []
        
        self.fBkgLeft = None 
        self.fBkgRight = None 
           
        for i in range(len(self.fitSlices)):

            # set the default Ranges:
            self.bkgEdges.append([0., 15.])
            self.bkgInnerEdges.append([0., 0.])
            self.signalEdges.append([0., 5.])

            self.fCombined.append( ROOT.TF1("fCombined_"+str(i),
                                            "(x <= ([3] + [1])) * [0] * TMath::Gaus(x, [1], [2]) + (x > ([3] + [1])) * [0] * TMath::Gaus([3] + [1], [1], [2]) * exp(-([3]) * (x - [3] - [1]) / ([2] * [2])) + [4] + [5] * exp([6]*x)", 
                                            self.bkgEdges[i][0],
                                            self.bkgEdges[i][1],
                                            7)
                                 )
            self.fSignal.append( ROOT.TF1(self.fittername+"_fSignal_"+str(i),
                                          "(x <= ([3] + [1])) * [0] * TMath::Gaus(x, [1], [2]) + (x > ([3] + [1])) * [0] * TMath::Gaus([3] + [1], [1], [2]) * exp(-([3]) * (x - [3] - [1]) / ([2] * [2]))",
                                          self.signalEdges[i][0],
                                          self.signalEdges[i][1],
                                          4)
                               )

            self.fBkg.append( ROOT.TF1("fBkg_"+str(i),
                                        "[0] + [1] * exp([2]*x)",
                                        self.bkgEdges[i][0],
                                        self.bkgEdges[i][1],
                                        3)
                            )

    def set_DoubleGauss_invMassFnc(self):
        
        self.signalEdges = []
        self.bkgInnerEdges = []
        self.bkgEdges = []
        
        self.fCombined = []
        self.fSignal = []
        self.fSignalPrefit = []
        self.fBkg = []
        self.fBkgPrefit = []
        self.fBkgLeft = []
        self.fBkgRight = []

        for i in range(len(self.fitSlices)):
            
            # set the default Ranges:
            if self.Particles(self.particleType) is self.Particles.Lambda:
                self.bkgEdges.append([1.284, 1.36])
                self.bkgInnerEdges.append([1.310, 1.330])
                self.signalEdges.append([1.310, 1.330])
            elif self.Particles(self.particleType) is self.Particles.Xi:
                self.bkgEdges.append([1.09, 1.14])
                self.bkgInnerEdges.append([1.110, 1.121])
                self.signalEdges.append([1.1117, 1.1197])
            else:
                print("This histogram is not supported yet")

            # set the fit functions and the corresponding QA functions
            self.fCombined.append( ROOT.TF1("fCombined_"+str(i),
                                            "[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [5]) + [6] + [7]*x + [8]*x*x + [9]*x*x*x",
                                            self.bkgEdges[i][0],
                                            self.bkgEdges[i][1],
                                            10)
                                 )
            self.fCombined[i].SetParLimits(0, 0., 1e10)
            self.fCombined[i].SetParLimits(3, 0., 1e10)

            self.fSignalPrefit.append([])  
            self.fSignalPrefit[i].append( ROOT.TF1("fSignal_"+str(i)+"_narrowGauss",
                                          "[0]*TMath::Gaus(x, [1], [2])",
                                          self.signalEdges[i][0],
                                          self.signalEdges[i][1],
                                          3)
                               )
            self.fSignalPrefit[i].append( ROOT.TF1("fSignal_"+str(i)+"_wideGauss",
                                          "[0]*TMath::Gaus(x, [1], [2])",
                                          self.signalEdges[i][0],
                                          self.signalEdges[i][1],
                                          3)
                               )
            self.fSignalPrefit[i][0].SetParLimits(0, 0., 1e10)
            self.fSignalPrefit[i][0].SetParLimits(3, 0., 1e10)
            self.fSignalPrefit[i][1].SetParLimits(0, 0., 1e10)
            self.fSignalPrefit[i][1].SetParLimits(3, 0., 1e10)

            self.fSignal.append( ROOT.TF1("fSignal_"+str(i),
                                          "[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [5])",
                                          self.signalEdges[i][0],
                                          self.signalEdges[i][1],
                                          6)
                               )
            self.fSignal[i].SetParLimits(0, 0., 1e10)
            self.fSignal[i].SetParLimits(3, 0., 1e10)
            
            self.fBkg.append( ROOT.TF1("fBkg_"+str(i),
                                        "[0] + [1]*x + [2]*x*x + [3]*x*x*x",
                                        self.bkgEdges[i][0],
                                        self.bkgEdges[i][1],
                                        4)
                            )

            self.fBkgPrefit.append( ROOT.TF1("fBkgPrefit_"+str(i),
                                    invMassLambdaBkgPrefitFunction,
                                    self.bkgEdges[i][0],
                                    self.bkgEdges[i][1],
                                    6)
                                  )


            self.fBkgLeft.append( ROOT.TF1("fBkgLeft_"+str(i),
                                          "[0] + [1]*x + [2]*x*x + [3]*x*x*x",
                                          self.bkgEdges[i][0],
                                          self.bkgInnerEdges[i][0],
                                          4)
                                )

            self.fBkgRight.append( ROOT.TF1("fBkgRight_"+str(i),
                                          "[0] + [1]*x + [2]*x*x + [3]*x*x*x",
                                          self.bkgInnerEdges[i][1],
                                          self.bkgEdges[i][1],
                                          4)
                                )
    
    def set_Lambda_invMassFnc(self):
        """
        depricated and kept for compatibility
        use the general invarian mass functions instead
        TODO: try to remove this as soon as possible
        """
        
        self.signalEdges = []
        self.bkgInnerEdges = []
        self.bkgEdges = []
        
        self.fCombined = []
        self.fSignal = []
        self.fSignalPrefit = []
        self.fBkg = []
        self.fBkgPrefit = []
        self.fBkgLeft = []
        self.fBkgRight = []

        for i in range(len(self.fitSlices)):
            
            # set the default Ranges:
            self.bkgEdges.append([1.09, 1.14])
            self.bkgInnerEdges.append([1.110, 1.121])
            self.signalEdges.append([1.1117, 1.1197])

            # set the fit functions and the corresponding QA functions
            self.fCombined.append( ROOT.TF1("fCombined_"+str(i),
                                            "[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [5]) + [6] + [7]*x + [8]*x*x + [9]*x*x*x",
                                            self.bkgEdges[i][0],
                                            self.bkgEdges[i][1],
                                            10)
                                 )
            self.fCombined[i].SetParLimits(0, 0., 1e10)
            self.fCombined[i].SetParLimits(3, 0., 1e10)

            self.fSignalPrefit.append([])  
            self.fSignalPrefit[i].append( ROOT.TF1("fSignal_"+str(i)+"_narrowGauss",
                                          "[0]*TMath::Gaus(x, [1], [2])",
                                          self.signalEdges[i][0],
                                          self.signalEdges[i][1],
                                          3)
                               )
            self.fSignalPrefit[i].append( ROOT.TF1("fSignal_"+str(i)+"_wideGauss",
                                          "[0]*TMath::Gaus(x, [1], [2])",
                                          self.signalEdges[i][0],
                                          self.signalEdges[i][1],
                                          3)
                               )
            self.fSignalPrefit[i][0].SetParLimits(0, 0., 1e10)
            self.fSignalPrefit[i][0].SetParLimits(3, 0., 1e10)
            self.fSignalPrefit[i][1].SetParLimits(0, 0., 1e10)
            self.fSignalPrefit[i][1].SetParLimits(3, 0., 1e10)

            self.fSignal.append( ROOT.TF1("fSignal_"+str(i),
                                          "[0]*TMath::Gaus(x, [1], [2]) + [3]*TMath::Gaus(x, [4], [5])",
                                          self.signalEdges[i][0],
                                          self.signalEdges[i][1],
                                          6)
                               )
            self.fSignal[i].SetParLimits(0, 0., 1e10)
            self.fSignal[i].SetParLimits(3, 0., 1e10)
            
            self.fBkg.append( ROOT.TF1("fBkg_"+str(i),
                                        "[0] + [1]*x + [2]*x*x + [3]*x*x*x",
                                        self.bkgEdges[i][0],
                                        self.bkgEdges[i][1],
                                        4)
                            )

            self.fBkgPrefit.append( ROOT.TF1("fBkgPrefit_"+str(i),
                                    invMassLambdaBkgPrefitFunction,
                                    self.bkgEdges[i][0],
                                    self.bkgEdges[i][1],
                                    6)
                                  )


            self.fBkgLeft.append( ROOT.TF1("fBkgLeft_"+str(i),
                                          "[0] + [1]*x + [2]*x*x + [3]*x*x*x",
                                          self.bkgEdges[i][0],
                                          self.bkgInnerEdges[i][0],
                                          4)
                                )

            self.fBkgRight.append( ROOT.TF1("fBkgRight_"+str(i),
                                          "[0] + [1]*x + [2]*x*x + [3]*x*x*x",
                                          self.bkgInnerEdges[i][1],
                                          self.bkgEdges[i][1],
                                          4)
                                )

    # Automatically sets one of the above defined functions
    def set_fitFunctions(self):

        #TODO: implement checks so that the fitting is only executed if all parameters are available
        #TODO: print the configurations of the fit when executing

        if self.Histograms(self.histType) is self.Histograms.tpcpid:
            self.set_tpcFitFnc()
        
        elif self.Histograms(self.histType) is self.Histograms.tpcpidtail:
            self.set_tpcFitFncTail()
        
        elif self.Histograms(self.histType) is self.Histograms.tofpid:
            self.set_tofFitFncTail()

        elif self.Histograms(self.histType) is self.Histograms.combpid_simple:
            self.set_combFitFnc_simple()
        
        elif self.Histograms(self.histType) is self.Histograms.combpid_gausBkg:
            self.set_combFitFnc_gausBkg()
        
        elif self.Histograms(self.histType) is self.Histograms.combpid:
            self.set_combFitFnc()

        elif self.Histograms(self.histType) is self.Histograms.invMass:
            self.set_DoubleGauss_invMassFnc()
            #self.set_Lambda_invMassFnc()

        #TODO: implement support for TOF fit function
        #TODO: implement support for invMass fit function
        else:
            print("this histogramm is not supported yet")

        
        for i, [fcomb, fsign, fbkg] in enumerate(zip(self.fCombined, self.fSignal, self.fBkg)):
            fcomb.SetLineColor(1)

            fsign.SetLineColor(ROOT.kGreen)
            fsign.SetLineWidth(2)

            fbkg.SetLineColor(2)

            if self.fBkgLeft is not None:
                self.fBkgLeft[i].SetLineColor(2)
                self.fBkgLeft[i].SetLineStyle(2)
            if self.fBkgRight is not None:
                self.fBkgRight[i].SetLineColor(2)
                self.fBkgRight[i].SetLineStyle(2)
                


    ############################################
    ########    set parameter limits ###########
    ############################################

    def set_fitrange(self, rangebkg, rangebkgInner, rangesign, ifit):
        
        if rangebkg is not None:
            self.bkgEdges[ifit] = rangebkg
            self.fCombined[ifit].SetRange(self.bkgEdges[ifit][0], self.bkgEdges[ifit][1])
            self.fBkg[ifit].SetRange(self.bkgEdges[ifit][0], self.bkgEdges[ifit][1])
            if self.Histograms(self.histType) is not self.Histograms.combpid_simple and self.Histograms(self.histType) is not self.Histograms.combpid and self.Histograms(self.histType) is not self.Histograms.combpid_gausBkg: #TODO: maybe need ot change the or into an and
                self.fBkgPrefit[ifit].SetRange(self.bkgEdges[ifit][0], self.bkgEdges[ifit][1])
        
        if rangebkgInner is not None:
            self.bkgInnerEdges[ifit] = rangebkgInner
        
        if rangesign is not None:
            self.signalEdges[ifit] = rangesign
            self.fSignal[ifit].SetRange(self.signalEdges[ifit][0], self.signalEdges[ifit][1])

            if self.Histograms(self.histType) is self.Histograms.invMass:
                self.fSignalPrefit[ifit][0].SetRange(self.signalEdges[ifit][0], self.signalEdges[ifit][1])
                self.fSignalPrefit[ifit][1].SetRange(self.signalEdges[ifit][0], self.signalEdges[ifit][1])


    ###############################
    ###############################


    def update_tpcQAFncs(self):
        for i in range(len(self.fitSlices)):
            self.fSignal[i].SetParameter(0, self.fCombined[i].GetParameter(0))
            self.fSignal[i].SetParameter(1, self.fCombined[i].GetParameter(1))
            self.fSignal[i].SetParameter(2, self.fCombined[i].GetParameter(2))

            self.fBkg[i].SetParameter(0, self.fCombined[i].GetParameter(3))
            self.fBkg[i].SetParameter(1, self.fCombined[i].GetParameter(4))
            self.fBkg[i].SetParameter(2, self.fCombined[i].GetParameter(5))
            self.fBkg[i].SetParameter(3, self.fCombined[i].GetParameter(6))
    
    def update_tpcQAFncsTail(self):
        for i in range(len(self.fitSlices)):
            self.fSignal[i].SetParameter(0, self.fCombined[i].GetParameter(0))
            self.fSignal[i].SetParameter(1, self.fCombined[i].GetParameter(1))
            self.fSignal[i].SetParameter(2, self.fCombined[i].GetParameter(2))
            self.fSignal[i].SetParameter(3, self.fCombined[i].GetParameter(3))

            self.fBkg[i].SetParameter(0, self.fCombined[i].GetParameter(4))
            self.fBkg[i].SetParameter(1, self.fCombined[i].GetParameter(5))
            self.fBkg[i].SetParameter(2, self.fCombined[i].GetParameter(6))
            self.fBkg[i].SetParameter(3, self.fCombined[i].GetParameter(7))
    
    def update_tofQAFncsTail(self):
    
        for i in range(len(self.fitSlices)):
            self.fSignal[i].SetParameter(0, self.fCombined[i].GetParameter(0))
            self.fSignal[i].SetParameter(1, self.fCombined[i].GetParameter(1))
            self.fSignal[i].SetParameter(2, self.fCombined[i].GetParameter(2))
            self.fSignal[i].SetParameter(3, self.fCombined[i].GetParameter(3))

            self.fBkg[i].SetParameter(0, self.fCombined[i].GetParameter(4))
            self.fBkg[i].SetParameter(1, self.fCombined[i].GetParameter(5))
            self.fBkg[i].SetParameter(2, self.fCombined[i].GetParameter(6))
            self.fBkg[i].SetParameter(3, self.fCombined[i].GetParameter(7))
            
            """
            self.fBkgRight[i].SetParameter(0, self.fBkg[i].GetParameter(0))
            self.fBkgRight[i].SetParameter(1, self.fBkg[i].GetParameter(1))
            self.fBkgLeft[i].SetParameter(0, self.fBkg[i].GetParameter(2))
            self.fBkgLeft[i].SetParameter(1, self.fBkg[i].GetParameter(3))
            """

            #self.fBkg[i].SetParameter(4, self.fCombined[i].GetParameter(8))
    
    def update_combQAFncs_simple(self):

        for i in range(len(self.fitSlices)):
            self.fSignal[i].SetParameter(0, self.fCombined[i].GetParameter(0))
            self.fSignal[i].SetParameter(1, self.fCombined[i].GetParameter(1))
            self.fSignal[i].SetParameter(2, self.fCombined[i].GetParameter(2))

            self.fBkg[i].SetParameter(0, self.fCombined[i].GetParameter(3))
            self.fBkg[i].SetParameter(1, self.fCombined[i].GetParameter(4))
            self.fBkg[i].SetParameter(2, self.fCombined[i].GetParameter(5))
            self.fBkg[i].SetParameter(3, self.fCombined[i].GetParameter(6))
    
    def update_combQAFncs_gausBkg(self):

        for i in range(len(self.fitSlices)):
            self.fSignal[i].SetParameter(0, self.fCombined[i].GetParameter(0))
            self.fSignal[i].SetParameter(1, self.fCombined[i].GetParameter(1))
            self.fSignal[i].SetParameter(2, self.fCombined[i].GetParameter(2))

            self.fBkg[i].SetParameter(0, self.fCombined[i].GetParameter(3))
            self.fBkg[i].SetParameter(1, self.fCombined[i].GetParameter(4))
            self.fBkg[i].SetParameter(2, self.fCombined[i].GetParameter(5))
            self.fBkg[i].SetParameter(3, self.fCombined[i].GetParameter(6))
            self.fBkg[i].SetParameter(4, self.fCombined[i].GetParameter(7))
            self.fBkg[i].SetParameter(5, self.fCombined[i].GetParameter(8))

    def update_combQAFncs(self):

        for i in range(len(self.fitSlices)):
            self.fSignal[i].SetParameter(0, self.fCombined[i].GetParameter(0))
            self.fSignal[i].SetParameter(1, self.fCombined[i].GetParameter(1))
            self.fSignal[i].SetParameter(2, self.fCombined[i].GetParameter(2))
            self.fSignal[i].SetParameter(3, self.fCombined[i].GetParameter(3))
            self.fSignal[i].SetParameter(4, self.fCombined[i].GetParameter(4))
            self.fSignal[i].SetParameter(5, self.fCombined[i].GetParameter(5))

            self.fBkg[i].SetParameter(0, self.fCombined[i].GetParameter(6))
            self.fBkg[i].SetParameter(1, self.fCombined[i].GetParameter(7))
            self.fBkg[i].SetParameter(2, self.fCombined[i].GetParameter(8))

    def update_invMassLambdaFnc(self):

        for i in range(len(self.fitSlices)):
            self.fSignal[i].SetParameter(0, self.fCombined[i].GetParameter(0))
            self.fSignal[i].SetParameter(1, self.fCombined[i].GetParameter(1))
            self.fSignal[i].SetParameter(2, self.fCombined[i].GetParameter(2))
            self.fSignal[i].SetParameter(3, self.fCombined[i].GetParameter(3))
            self.fSignal[i].SetParameter(4, self.fCombined[i].GetParameter(4))
            self.fSignal[i].SetParameter(5, self.fCombined[i].GetParameter(5))

            self.fBkg[i].SetParameter(0, self.fCombined[i].GetParameter(6))
            self.fBkg[i].SetParameter(1, self.fCombined[i].GetParameter(7))
            self.fBkg[i].SetParameter(2, self.fCombined[i].GetParameter(8))
            self.fBkg[i].SetParameter(3, self.fCombined[i].GetParameter(9))



            self.fSignalPrefit[i][0].SetParameter(0, self.fCombined[i].GetParameter(0))
            self.fSignalPrefit[i][0].SetParameter(1, self.fCombined[i].GetParameter(1))
            self.fSignalPrefit[i][0].SetParameter(2, self.fCombined[i].GetParameter(2))
            
            self.fSignalPrefit[i][1].SetParameter(0, self.fCombined[i].GetParameter(3))
            self.fSignalPrefit[i][1].SetParameter(1, self.fCombined[i].GetParameter(4))
            self.fSignalPrefit[i][1].SetParameter(2, self.fCombined[i].GetParameter(5))

    
    
    def update_qaFncs(self):
        """
        - this is a bit unecessary. After the paper proposal it should be updates such that it automatically updates the parameters by getting the number of parameters
        - optionally include an array of parameters to be skipped
        """
        if self.Histograms(self.histType) is self.Histograms.tpcpid:
            self.update_tpcQAFncs()

        elif self.Histograms(self.histType) is self.Histograms.tpcpidtail:
            self.update_tpcQAFncsTail()

        elif self.Histograms(self.histType) is self.Histograms.tofpid:
            self.update_tofQAFncsTail()

        elif self.Histograms(self.histType) is self.Histograms.combpid_simple:
            self.update_combQAFncs_simple()
        
        elif self.Histograms(self.histType) is self.Histograms.combpid_gausBkg:
            self.update_combQAFncs_gausBkg()
        
        elif self.Histograms(self.histType) is self.Histograms.combpid:
            self.update_combQAFncs()

        elif self.Histograms(self.histType) is self.Histograms.invMass:
            self.update_invMassLambdaFnc()
    ###############################

    def fixSinglebkgpar(self, ifit, par, val):
        offset = 0
        if self.Histograms(self.histType) is self.Histograms.combpid_simple or self.Histograms(self.histType) is self.Histograms.combpid_gausBkg:
            offset = 3 #TODO: in future this value should be retrieved automatically from the function --> Save it as an class atribute
        elif self.Histograms(self.histType) is self.Histograms.invMass:
            offset = 6  
        self.fBkg[ifit].FixParameter(par, val)
        self.fCombined[ifit].FixParameter(offset+par, val)
    
    def setSinglebkgparlimits(self, ifit, par, vallow, valup):
        offset = 0
        if self.Histograms(self.histType) is self.Histograms.combpid_simple or self.Histograms(self.histType) is self.Histograms.combpid_gausBkg:
            offset = 3 #TODO: in future this value should be retrieved automatically from the function --> Save it as an class atribute
        
        self.fBkg[ifit].SetParLimits(par, vallow, valup)
        self.fCombined[ifit].SetParLimits(offset+par, vallow, valup)
        

    def fixAllbkgpar(self, par, val):
        offset = 0
        if self.Histograms(self.histType) is self.Histograms.combpid_simple or self.Histograms(self.histType) is self.Histograms.combpid_gausBkg:
            offset = 3 #TODO: in future this value should be retrieved automatically from the function --> Save it as an class atribute
        elif self.Histograms(self.histType) is self.Histograms.invMass:
            offset = 6  
        
        for fnc, bkgfnc, bkgfncPrefit in zip(self.fCombined, self.fBkg, self.fBkgPrefit):
            bkgfnc.FixParameter(par, val)
            bkgfncPrefit.FixParameter(par+2, val)
            fnc.FixParameter(offset+par, val)

    
    def setAllbkgparlimits(self, par, vallow, valup):
        offset = 0
        if self.Histograms(self.histType) is self.Histograms.combpid_simple or self.Histograms(self.histType) is self.Histograms.combpid_gausBkg:
            offset = 3 #TODO: in future this value should be retrieved automatically from the function --> Save it as an class atribute
        
        for fnc, bkgfnc in zip(self.fCombined, self.fBkg):
            bkgfnc.SetParLimits(par, vallow, valup)
            fnc.SetParLimits(offset+par, vallow, valup)

    def prefitAllsignal(self, option="QLRN", fix=False):

        for i, [ hist, fnc ] in enumerate(zip(self.fitSlices, self.fSignal)):

            if self.Histograms(self.histType) is self.Histograms.tpcpid:
                print("TPCPID not supported yet") 

            elif self.Histograms(self.histType) is self.Histograms.tpcpidtail or self.Histograms(self.histType) is self.Histograms.combpid:
                
                fnc.SetParameter(0, self.fitSlices[i].GetBinContent(self.fitSlices[i].GetXaxis().FindBin(0.)))
                fnc.SetParameter(1, 0.)
                fnc.SetParameter(2, 1.2)
                hist.Fit(fnc, option, "", self.signalEdges[i][0], self.signalEdges[i][1])

                if fix:
                    self.fCombined[i].FixParameter(0, fnc.GetParameter(0))
                    self.fCombined[i].FixParameter(1, fnc.GetParameter(1))
                    self.fCombined[i].FixParameter(2, fnc.GetParameter(2))
                    self.fCombined[i].FixParameter(3, fnc.GetParameter(3))
                else: 
                    self.fCombined[i].SetParameter(0, fnc.GetParameter(0))
                    self.fCombined[i].SetParameter(1, fnc.GetParameter(1))
                    self.fCombined[i].SetParameter(2, fnc.GetParameter(2))
                    self.fCombined[i].SetParameter(3, fnc.GetParameter(3))
            
            elif self.Histograms(self.histType) is self.Histograms.tofpid:
                
                fnc.SetParameter(0, self.fitSlices[i].GetBinContent(self.fitSlices[i].GetXaxis().FindBin(1.0)))
                fnc.SetParameter(1, 1.0)
                fnc.SetParameter(2, 2.)
                hist.Fit(fnc, option, "", self.signalEdges[i][0], self.signalEdges[i][1])

                self.fCombined[i].SetParameter(0, fnc.GetParameter(0))
                self.fCombined[i].SetParameter(1, fnc.GetParameter(1))
                self.fCombined[i].SetParameter(2, fnc.GetParameter(2))
                self.fCombined[i].SetParameter(3, fnc.GetParameter(3))
                #self.fCombined[i].SetParameter(4, fnc.GetParameter(4))
                #self.fCombined[i].SetParameter(5, fnc.GetParameter(5))

            elif self.Histograms(self.histType) is self.Histograms.combpid_simple or self.Histograms(self.histType) is self.Histograms.combpid_gausBkg:

                # double exponential gauss:
                fnc.SetParameter(0, 0.5*self.fitSlices[i].GetBinContent(self.fitSlices[i].GetXaxis().FindBin(1.0)))
                fnc.SetParameter(1, 1.2)
                fnc.SetParameter(2, 3.)
                
                hist.Fit(fnc, option, "", self.signalEdges[i][0], self.signalEdges[i][1])

                self.fCombined[i].SetParameter(0, fnc.GetParameter(0))
                self.fCombined[i].SetParameter(1, fnc.GetParameter(1))
                self.fCombined[i].SetParameter(2, fnc.GetParameter(2))

            elif self.Histograms(self.histType) is self.Histograms.invMass:
                
                if self.Particles(self.particleType) is self.Particles.Lambda:
                    mPDG = 1.1157
                    self.fSignalPrefit[i][0].SetParameter(0, hist.GetBinContent( hist.GetXaxis().FindBin(mPDG) )) 
                    self.fSignalPrefit[i][0].SetParameter(1, mPDG) 
                    self.fSignalPrefit[i][0].SetParameter(2, 0.004) 
                elif self.Particles(self.particleType) is self.Particles.Lambda:
                    mPDG = 1.3217
                    self.fSignalPrefit[i][0].SetParameter(0, hist.GetBinContent( hist.GetXaxis().FindBin(mPDG) )) 
                    self.fSignalPrefit[i][0].SetParameter(1, mPDG) 
                    self.fSignalPrefit[i][0].SetParameter(2, 0.005) 

                hist.Fit(self.fSignalPrefit[i][0], "QLRN", "", self.signalEdges[i][0], self.signalEdges[i][1])

                fnc.SetParameter(0, 0.8*self.fSignalPrefit[i][0].GetParameter(0))
                fnc.SetParameter(1, self.fSignalPrefit[i][0].GetParameter(1))
                fnc.SetParameter(2, 0.8*self.fSignalPrefit[i][0].GetParameter(2))
                
                if fnc.GetParameter(5) == 1234.567:
                    fnc.SetParameter(3, 0.)
                    fnc.SetParameter(4, 0.)
                    fnc.SetParameter(5, 1.)
                else:
                    fnc.SetParameter(3, 0.2*self.fSignalPrefit[i][0].GetParameter(0))
                    fnc.SetParameter(4, self.fSignalPrefit[i][0].GetParameter(1))
                    fnc.SetParameter(5, 2.0*self.fSignalPrefit[i][0].GetParameter(2))

                hist.Fit(fnc, option, "", self.signalEdges[i][0], self.signalEdges[i][1])

                self.fCombined[i].SetParameter(0, fnc.GetParameter(0))
                self.fCombined[i].SetParameter(1, fnc.GetParameter(1))
                self.fCombined[i].SetParameter(2, fnc.GetParameter(2))
                self.fCombined[i].SetParameter(3, fnc.GetParameter(3))
                self.fCombined[i].SetParameter(4, fnc.GetParameter(4))
                self.fCombined[i].SetParameter(5, fnc.GetParameter(5))

                self.fSignalPrefit[i][0].SetParameter(0, fnc.GetParameter(0))
                self.fSignalPrefit[i][0].SetParameter(1, fnc.GetParameter(1))
                self.fSignalPrefit[i][0].SetParameter(2, fnc.GetParameter(2))
                self.fSignalPrefit[i][1].SetParameter(3, fnc.GetParameter(3))
                self.fSignalPrefit[i][1].SetParameter(4, fnc.GetParameter(4))
                self.fSignalPrefit[i][1].SetParameter(5, fnc.GetParameter(5))

            """
            elif self.Histograms(self.histType) is self.Histograms.combpid:

                # double exponential gauss:
                fnc.SetParameter(0, 0.5*self.fitSlices[i].GetBinContent(self.fitSlices[i].GetXaxis().FindBin(1.0)))
                fnc.SetParameter(1, 1.2)
                fnc.SetParameter(2, 0.5)
                fnc.SetParameter(3, 0.5*self.fitSlices[i].GetBinContent(self.fitSlices[i].GetXaxis().FindBin(1.0)))
                fnc.SetParameter(4, 0.5)
                fnc.SetParameter(5, 0.5)
                
                hist.Fit(fnc, option, "", self.signalEdges[i][0], self.signalEdges[i][1])

                self.fCombined[i].SetParameter(0, fnc.GetParameter(0))
                self.fCombined[i].SetParameter(1, fnc.GetParameter(1))
                self.fCombined[i].SetParameter(2, fnc.GetParameter(2))
                self.fCombined[i].SetParameter(3, fnc.GetParameter(3))
                self.fCombined[i].SetParameter(4, fnc.GetParameter(4))
                self.fCombined[i].SetParameter(5, fnc.GetParameter(4))
                
            """

    def prefitAllbkg(self, option="QLRN"):
        
        for i, [ hist, fnc ] in enumerate(zip(self.fitSlices, self.fBkg)):

            if self.Histograms(self.histType) is self.Histograms.tpcpid:
                print("TPCPID not supported yet") 
            
            elif self.Histograms(self.histType) is self.Histograms.tpcpidtail:
                
                if self.prefitSeparatedBkg is True:

                    hist.Fit(self.fBkgLeft[i], option, "", self.bkgEdges[i][0], self.bkgInnerEdges[i][0])
                    hist.Fit(self.fBkgRight[i], option, "", self.bkgInnerEdges[i][1], self.bkgEdges[i][1])
                    
                    self.fCombined[i].SetParameter(4, self.fBkgRight[i].GetParameter(0))
                    self.fCombined[i].SetParameter(5, self.fBkgRight[i].GetParameter(1))
                    self.fCombined[i].SetParameter(6, self.fBkgLeft[i].GetParameter(0))
                    self.fCombined[i].SetParameter(7, self.fBkgLeft[i].GetParameter(1))
                    
                    fnc.SetParameter(0, self.fBkgRight[i].GetParameter(0))
                    fnc.SetParameter(1, self.fBkgRight[i].GetParameter(1))
                    fnc.SetParameter(2, self.fBkgLeft[i].GetParameter(0))
                    fnc.SetParameter(3, self.fBkgLeft[i].GetParameter(1))
                
                else:

                    self.fBkgPrefit[i].FixParameter(4, self.bkgInnerEdges[i][0])
                    self.fBkgPrefit[i].FixParameter(5, self.bkgInnerEdges[i][1])
                    
                    hist.Fit(self.fBkgPrefit[i], option, "", self.bkgEdges[i][0], self.bkgEdges[i][1])
                    
                    self.fCombined[i].SetParameter(4, self.fBkgPrefit[i].GetParameter(0))
                    self.fCombined[i].SetParameter(5, self.fBkgPrefit[i].GetParameter(1))
                    self.fCombined[i].SetParameter(6, self.fBkgPrefit[i].GetParameter(2))
                    self.fCombined[i].SetParameter(7, self.fBkgPrefit[i].GetParameter(3))
                    
                    fnc.SetParameter(0, self.fBkgPrefit[i].GetParameter(0))
                    fnc.SetParameter(1, self.fBkgPrefit[i].GetParameter(1))
                    fnc.SetParameter(2, self.fBkgPrefit[i].GetParameter(2))
                    fnc.SetParameter(3, self.fBkgPrefit[i].GetParameter(3))

            elif self.Histograms(self.histType) is self.Histograms.tofpid:
            
                # self.fBkgPrefit[i].FixParameter(4, self.bkgInnerEdges[i][0])
                # self.fBkgPrefit[i].FixParameter(5, self.bkgInnerEdges[i][1])
                 
                hist.Fit(fnc, option, "", self.bkgEdges[i][0], self.bkgEdges[i][1])
                
                self.fCombined[i].SetParameter(5, fnc.GetParameter(0))
                self.fCombined[i].SetParameter(6, fnc.GetParameter(1))
                self.fCombined[i].SetParameter(7, fnc.GetParameter(2))
                self.fCombined[i].SetParameter(8, fnc.GetParameter(3))
                
                
                """
                fnc.RejectPoint(self.bkgInnerEdges[i][0])
                fnc.RejectPoint(self.bkgInnerEdges[i][1])

                hist.Fit(fnc, "QLRN", "", self.bkgEdges[i][0], self.bkgEdges[i][1])

                self.fCombined[i].SetParameter(4, fnc.GetParameter(0))
                self.fCombined[i].SetParameter(5, fnc.GetParameter(1))
                self.fCombined[i].SetParameter(6, fnc.GetParameter(2))
                self.fCombined[i].SetParameter(7, fnc.GetParameter(3))
                """
            
            elif self.Histograms(self.histType) is self.Histograms.combpid_simple:

                #fnc.FixParameter(3, 0.)
                #fnc.SetParLimits(0, 0., 99999999.)
                #fnc.SetParLimits(0, -99999999., 0.)

                hist.Fit(fnc, option, "", self.bkgInnerEdges[i][0], self.bkgInnerEdges[i][1])
                self.fCombined[i].SetParameter(3, fnc.GetParameter(0))
                self.fCombined[i].SetParameter(4, fnc.GetParameter(1))
                self.fCombined[i].SetParameter(5, fnc.GetParameter(2))
                self.fCombined[i].SetParameter(6, fnc.GetParameter(3))

                #self.fCombined[i].FixParameter(6, 0.) #fix the parameter to use only a constant term + exponential
                #self.fCombined[i].SetParLimits(3, 0., 99999999.)
                #self.fCombined[i].SetParLimits(3, -99999999., 0.)
            
            elif self.Histograms(self.histType) is self.Histograms.combpid_gausBkg:
                    
                if self.prefitSeparatedBkg is True:

                    if self.swapBkgOrder is False:
                        hist.Fit(self.fBkgLeft[i], option, "", self.signalEdges[i][1], self.bkgInnerEdges[i][0])
                        hist.Fit(self.fBkgRight[i], option, "", self.bkgInnerEdges[i][0], self.bkgEdges[i][1])
                    else:
                        hist.Fit(self.fBkgLeft[i], option, "", self.bkgInnerEdges[i][0], self.bkgEdges[i][1])
                        hist.Fit(self.fBkgRight[i], option, "", self.signalEdges[i][1], self.bkgInnerEdges[i][0])

                    fnc.SetParameter(0, self.fBkgRight[i].GetParameter(0)) #in the definition of the function: fBkg = exp + gaus
                    fnc.SetParameter(1, self.fBkgRight[i].GetParameter(1))
                    fnc.SetParameter(2, self.fBkgLeft[i].GetParameter(0))
                    fnc.SetParameter(3, self.fBkgLeft[i].GetParameter(1))
                    fnc.SetParameter(4, self.fBkgLeft[i].GetParameter(2))

                    self.fCombined[i].SetParameter(3, fnc.GetParameter(0))
                    self.fCombined[i].SetParameter(4, fnc.GetParameter(1))
                    self.fCombined[i].SetParameter(5, fnc.GetParameter(2))
                    self.fCombined[i].SetParameter(6, fnc.GetParameter(3))
                    self.fCombined[i].SetParameter(6, fnc.GetParameter(4))
                else:
                    hist.Fit(fnc, option, "", self.bkgInnerEdges[i][0], self.bkgInnerEdges[i][1])
                    self.fCombined[i].SetParameter(3, fnc.GetParameter(0))
                    self.fCombined[i].SetParameter(4, fnc.GetParameter(1))
                    self.fCombined[i].SetParameter(5, fnc.GetParameter(2))
                    self.fCombined[i].SetParameter(6, fnc.GetParameter(3))
                    self.fCombined[i].SetParameter(6, fnc.GetParameter(4))
                
                #self.fCombined[i].FixParameter(6, 0.)


            elif self.Histograms(self.histType) is self.Histograms.combpid:

                hist.Fit(fnc, option, "", self.bkgEdges[i][0], self.bkgEdges[i][1])
                self.fCombined[i].SetParameter(4, fnc.GetParameter(0))
                self.fCombined[i].SetParameter(5, fnc.GetParameter(1))
                self.fCombined[i].SetParameter(6, fnc.GetParameter(2))

            elif self.Histograms(self.histType) is self.Histograms.invMass:

                self.fBkgPrefit[i].FixParameter(0, self.bkgInnerEdges[i][0])
                self.fBkgPrefit[i].FixParameter(1, self.bkgInnerEdges[i][1])
                
                hist.Fit(self.fBkgPrefit[i], option, "", self.bkgEdges[i][0], self.bkgEdges[i][1])
                
                self.fCombined[i].SetParameter(6, self.fBkgPrefit[i].GetParameter(2))
                self.fCombined[i].SetParameter(7, self.fBkgPrefit[i].GetParameter(3))
                self.fCombined[i].SetParameter(8, self.fBkgPrefit[i].GetParameter(4))
                self.fCombined[i].SetParameter(9, self.fBkgPrefit[i].GetParameter(5))

                self.fBkgPrefit[i].SetLineStyle(2)



    def fitAllSlices(self, option="LRN"):

        #self.fitResult = []
        for i, [ hist, fnc ] in enumerate(zip(self.fitSlices, self.fCombined)):
            print("Fitting slice "+str(i))
            hist.Fit(fnc, option, "", self.bkgEdges[i][0], self.bkgEdges[i][1])
        
        self.update_qaFncs() 

    def set_prefitSeparatedBkg(self, setting):
        """
        TODO: For now this has to be called AFTER setting the fit ranges
        """
        self.prefitSeparatedBkg = setting

        if(setting is True):
            for ifit in range(len(self.fitSlices)):
                self.fBkgLeft[ifit].SetRange(self.bkgEdges[ifit][0], self.bkgInnerEdges[ifit][0])
                self.fBkgRight[ifit].SetRange(self.bkgInnerEdges[ifit][1], self.bkgEdges[ifit][1])

    
    def get_NfitSlices(self):
        return self.NfitSlices

    def get_entries(self, lowerval, upperval):

        entries = []
        for i, fslice in enumerate(self.fitSlices):
            
            entries.append( fslice.Integral(fslice.GetXaxis().FindBin(lowerval), fslice.GetXaxis().FindBin(upperval)) )
        
        return entries

    def calc_purities(self, lowerval, upperval):

        self.purities = []
        for i, [fsig, fcomb] in enumerate(zip(self.fSignal, self.fCombined)):
            
            self.purities.append( fsig.Integral(lowerval, upperval)/fcomb.Integral(lowerval, upperval))
        
        return self.purities


    def calc_sigmas(self, lowerval, upperval):

        self.sigmas = []
        if self.Histograms(self.histType) is self.Histograms.invMass:
            tempIntegral1 = 0.
            tempIntegral2 = 0.
            for fsigprefit in self.fSignalPrefit:
                 #tempIntegral1 = fsigprefit[0].Integral(fsigprefit[0].FindBin(lowerval), fsigprefit[0].FindBin(upperval))
                 #tempIntegral2 = fsigprefit[1].Integral(fsigprefit[1].FindBin(lowerval), fsigprefit[1].FindBin(upperval))
                 tempIntegral1 = fsigprefit[0].Integral(lowerval, upperval)
                 tempIntegral2 = fsigprefit[1].Integral(lowerval, upperval)
                 
                 self.sigmas.append( (fsigprefit[0].GetParameter(2)*tempIntegral1 + fsigprefit[1].GetParameter(2)*tempIntegral2)/(tempIntegral1+tempIntegral2) )
        else:
            for fsig in self.fSignal:
                self.sigmas.append( fsig.GetParameter(2))
        
        return self.sigmas


    def calc_chi2s(self, lowerval_signal, upperval_signal):

        self.chi2_comb=[]
        self.ndf_comb=[]
        self.chi2_sig=[]
        self.ndf_sig=[]

        for i, [hist, fcomb, fsig] in enumerate(zip(self.fitSlices, self.fCombined, self.fSignal)):

            chi2_comb_val = 0.
            chi2_sig_val = 0.
            nbinsComb = 0
            nbinsSig = 0

            lowerbinBkg = hist.GetXaxis().FindBin(self.bkgEdges[i][0])
            upperbinBkg = hist.GetXaxis().FindBin(self.bkgEdges[i][1])
            lowerbinSig = hist.GetXaxis().FindBin(lowerval_signal)
            upperbinSig = hist.GetXaxis().FindBin(upperval_signal)

            for ibin in range(lowerbinBkg, upperbinBkg):
                chi2_comb_val += ((hist.GetBinContent(ibin)-fcomb.Eval(hist.GetXaxis().GetBinCenter(ibin)))**2)/((hist.GetBinError(ibin))**2)
                nbinsComb += 1

            for ibin in range(lowerbinSig, upperbinSig):
                chi2_sig_val += ((hist.GetBinContent(ibin)-fcomb.Eval(hist.GetXaxis().GetBinCenter(ibin)))**2)/((hist.GetBinError(ibin))**2)
                nbinsSig += 1

            self.chi2_comb.append( chi2_comb_val )
            self.ndf_comb.append(nbinsComb - fcomb.GetNpar())
            self.chi2_sig.append( chi2_sig_val ) 
            self.ndf_sig.append( (nbinsSig - fcomb.GetNpar()) )

        return [self.chi2_comb, self.ndf_comb, self.chi2_sig, self.ndf_sig]



    def printFitParameters(self, ifit):
        for ipar in range(self.fCombined[ifit].GetNpar()):
            print(f"Parameter {ipar}: {self.fCombined[ifit].GetParameter(ipar):.2e}")

    def rebin_fitSlices(self, rebinlist):
        for i, hist in enumerate(self.fitSlices):
            hist.Rebin(rebinlist[i])

    def rebin_fitSlice(self, islice, rebinfactor):
        self.fitSlices[islice].Rebin(rebinfactor)

    def set_swapBkgOrder(self, setting):
        self.swapBkgOrder = setting

   # def calc_chi2(self):





####################################
########## Draw Functions ##########
####################################


    def set_discretize(self, discretize, nbinsComb, nbinsSignal):
        self.discretizeResults = discretize

        self.hCombined = []
        self.hSignal = []

        if self.discretizeResults:
            for i, [fcomb, fsig] in enumerate(zip(self.fCombined, self.fSignal)):
                self.hCombined.append( ROOT.TH1F("hCombined_"+str(i), "hCombined_"+str(i), nbinsComb, self.bkgEdges[i][0], self.bkgEdges[i][1]) )
                self.hSignal.append( ROOT.TH1F("hSignal_"+str(i), "hSignal_"+str(i), nbinsSignal, self.signalEdges[i][0], self.signalEdges[i][1]) )

                for ibin in range(1, nbinsComb+1):
                    self.hCombined[i].SetBinContent( ibin, fcomb.Eval( self.hCombined[i].GetXaxis().GetBinCenter(ibin) ) )
                    self.hCombined[i].SetLineColor(1)
                
                for ibin in range(1, nbinsSignal):
                    self.hSignal[i].SetBinContent( ibin, fsig.Eval( self.hSignal[i].GetXaxis().GetBinCenter(ibin) ) )
                    self.hSignal[i].SetLineColor(5)
                    self.hSignal[i].SetLineWidth(2)


    def draw_fitSlices(self, canvas, option = "", howDivide=None, drawOffset=0):

        print(canvas.GetListOfPrimitives().GetSize())
        if canvas.GetListOfPrimitives().GetSize() == 0:
            print("Dividing Canvas")
            if howDivide is None:
                canvas.Divide(len(self.fitSlices))
            else:
                canvas.Divide(howDivide[0], howDivide[1])

        for i, hist in enumerate(self.fitSlices, 1):
            canvas.cd(i+drawOffset)
            hist.DrawClone(option)

    def draw_fitResults(self, canvas, option="same", skip = [ False, False, False ], drawoffset=0):

        print(canvas.GetListOfPrimitives().GetSize())
        if canvas.GetListOfPrimitives().GetSize() == 0:
            print("Dividing Canvas")
            canvas.Divide(len(self.fitSlices))

        for i, [fcomb, fsign, fbkg] in enumerate(zip(self.fCombined, self.fSignal, self.fBkg)):
            canvas.cd(i+1+drawoffset)
            if not skip[0]:
                if self.discretizeResults:
                    self.hCombined[i].DrawClone(option)
                else:
                    fcomb.DrawClone(option)

            if not skip[1]:
                if self.discretizeResults:
                    self.hSignal[i].DrawClone(option)
                else:
                    fsign.DrawClone(option)

            if not skip[2]:
                fbkg.DrawClone(option)

    def draw_bkgPrefit(self, canvas, option="same", drawOffset=0):
        for i, fbkgprefit in enumerate(self.fBkgPrefit):
            canvas.cd(i+1+drawOffset)
            fbkgprefit.DrawClone(option)

    def draw_signalPrefit(self, canvas, option="same", drawOffset=0):

        for ifit, fsigprefit in enumerate(self.fSignalPrefit):
            canvas.cd(ifit+1+drawOffset)
            for i in range(len(fsigprefit)):
                fsigprefit[i].DrawClone(option)

    def draw_SeparateBkgPrefits(self, canvas, option=""):
        for i, [fbkgLeft, fbkgRight] in enumerate(zip(self.fBkgLeft, self.fBkgRight)):
            canvas.cd(i+1) 
            fbkgLeft.DrawClone("same "+option)
            fbkgRight.DrawClone("same "+option)

    def draw_SomeNumbers(self, canvas, interval, pos=[0.13, 0.2, 0.53, 0.3], drawOffset=0):

        thepurities = self.calc_purities(interval[0], interval[1])
        thesigmas = self.calc_sigmas(interval[0], interval[1])

        for i, [ipur, isig] in enumerate(zip(thepurities, thesigmas)):
            
            canvas.cd(i+1+drawOffset)
            text = ROOT.TPaveText(pos[0], pos[1], pos[2], pos[3], "NDC")
            text.SetBorderSize(0)
            text.SetFillColor(0)
            text.SetFillStyle(0)
            text.SetTextFont(42)
            text.AddText("Purity: {:.4f}".format(ipur))
            text.AddText("Sigma: {:.2e}".format(isig))
            text.DrawClone()




    def drawSignalLines(self, canvas, lowerlimit, upperlimit, drawOffset=0):

        for i, hist in enumerate(self.fitSlices):
            canvas.cd(i+1+drawOffset)
            
            lowerline = ROOT.TLine(lowerlimit, 0, lowerlimit, hist.GetMaximum())
            lowerline.SetLineWidth(2)
            lowerline.SetLineColor(ROOT.kBlack)
            lowerline.SetLineStyle(2)
            upperline = ROOT.TLine(upperlimit, 0, upperlimit, hist.GetMaximum())
            upperline.SetLineWidth(2)
            upperline.SetLineColor(ROOT.kBlack)
            upperline.SetLineStyle(2)
            

            lowerline.DrawClone()
            upperline.DrawClone()
            


####################################
#######   Getter Functions #########
####################################
    def get_ptConfig(self):
        return self.ptConfig
    
    def get_ptEdges(self):
        ptBins = []
        for ibin in range(1, len(self.ptConfig)):
            ptBins.append(self.ptConfig[ibin][0])
        ptBins.append(self.ptConfig[-1][1])
        
        return ptBins

    def get_ptValues(self):

        ptBins = self.get_ptEdges()
        ptValues = []
        for ibin in ptBins:
            ptValues.append(self.fitvarVSpt.GetXaxis().GetBinLowEdge(ibin))

        return ptValues

    def get_fitvarVSpt(self):
        return self.fitvarVSpt

    def get_particleType(self):
        return self.Particles(self.particleType)

    def get_fitSlice(self, ifit):
        return self.fitSlices[ifit]

    def get_fitSlices(self):
        return self.fitSlices

    def get_fCombined(self, i):
        return self.fCombined[i]

    def get_fCombinedAll(self):
        return self.fCombined
    
    def get_fSignal(self, i):
        return self.fSignal[i]
    
    def get_fSignalAll(self):
        return self.fSignal
    
    def get_fBkg(self, i):
        return self.fBkg[i]

    def get_fBkgAll(self):
        return self.fBkg
    
    def get_fBkgPrefit(self, i):
        return self.fBkgPrefit[i]

    def get_fBkgPrefitAll(self):
        return self.fBkgPrefit
    
    
    def get_fBkgLeft(self, i):
        return self.fBkgLeft[i]
    
    def get_fBkgLeftAll(self):
        return self.fBkgLeft
    
    def get_fBkgRight(self, i):
        return self.fBkgRight[i]
    
    def get_fBkgRightAll(self):
        return self.fBkgRight

    def get_nSlices(self):
        return len(self.fitSlices)




























    def set_fitLimits(self, singalLimits, bkgLimits):

        self.signalEdges = []
        self.bkgEdges = []
        self.bkgInnerEdges = []

        if (len(singalLimits) != len(self.fitSlices) ) or (len(bkgLimits) != (len(self.fitSlices))):
            print("You have not the mathing number of bin edges")
        for signal, bkg in zip(singalLimits, bkgLimits):
            self.signalEdges.append(signal)
            self.bkgEdges.append(bkg)


# This is not really needed anymore because with the getter functions one can get a pointer to the object and set the limits in the analysis script/ notebook 
    ############################################
    ########    set parameter limits ###########
    ############################################
    def set_limitValue(self, limitVal, ifit, ipar):
        self.fCombined[ifit].SetParLimits(ipar, limitVal[0], limitVal[1])

    def set_limitValues(self, limitVals, ifit):
        for ipar in range(self.fCombined[ifit].GetNpar()):
            self.set_limitValue(limitVals[ipar], ifit, ipar)

    def set_limitValuesAll(self, limitValuesAll):
        for ifit in range(len(self.fitSlices)):
            self.set_limitValues(limitValuesAll[ifit], ifit)

    ############################################
    #######    initialize parameters ###########
    ############################################
    def set_initValue(self, initVal, ifit, ipar):
        self.fCombined[ifit].SetParameter(ipar, initVal)
    
    def set_initValues(self, initVals, ifit):
        #print("Number of parameters: "+str(self.fCombined[ifit].GetNpar()))
        for ipar in range(self.fCombined[ifit].GetNpar()):
            #print(ipar)
            if initVals[ipar] is not None:
                self.set_initValue(initVals[ipar], ifit, ipar)
    
    def set_initValuesAll(self, initValsAll):
        #print("Number of fitSlices: "+str(len(self.fitSlices)))
        for ifit in range(len(self.fitSlices)):
            #print(ifit)
            self.set_initValues(initValsAll[ifit], ifit)







"""
for ibin in range(1, temphist.GetNbinsX()):
    if ibin > lowerbin and ibin < upperbin:
        temphist.SetBinError(ibin, hist.GetBinContent(hist.GetXaxis().FindBin(0.))*10.)
    else:
        temphist.SetBinError(ibin, 0.)

lin_x1 = self.bkgEdges[i][0]
lin_x2 = self.bkgEdges[i][1]
lin_y1 = hist.GetBinContent(hist.GetXaxis().FindBin(lin_x1))
lin_y2 = hist.GetBinContent(hist.GetXaxis().FindBin(lin_x2))

#fnc.SetParameter(0, 0.)
#fnc.SetParameter(1, 0.)
fnc.SetParameter(2, lin_y2 - lin_y1)
fnc.SetParameter(3, ROOT.TMath.Exp((lin_y2 - lin_y1)/(lin_x2 - lin_x1)) )

#fnc.SetParameter(2, hist.GetBinContent( hist.FindBin(self.signalEdges[i][0]) ) )
#fnc.SetParameter(3, 0.)

temphist.Fit(fnc, option, "", self.bkgEdges[i][0], self.bkgEdges[i][1])
self.fCombined[i].SetParameter(4, fnc.GetParameter(0))
self.fCombined[i].SetParameter(5, fnc.GetParameter(1))
self.fCombined[i].SetParameter(6, fnc.GetParameter(2))
self.fCombined[i].SetParameter(7, fnc.GetParameter(3))
"""