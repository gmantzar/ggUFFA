import ROOT

class Comparer:
    
    def __init__(self, nMult, nMt, nParallelFits, compname=""):

        self.compname = compname 
        self.nMult = nMult   #bins in multiplicity
        self.nMt = nMt       #bins in mt
        self.nParallelFits = nParallelFits  #how many fits are executed in parallel (e.g. pp and apap) They should be in common for all lists of fit results
        
        
        #TODO: These two offset values might be usefull to become parameters
        self.nMultOffset = 0
        self.nMtOffset = 2
        
        #self.data = [[]]*self.nParallelFits
        self.data = []
        self.data_extra = []
            
        self.counter = 0
        
        self.namelist = []
        self.descriptionlist = []
        
        #fit results
        self.fitlist = []  #fit results to compare
        self.Bllist = []
        self.radii = []
        self.radiiErr = []
        self.mTscalingList = []
        
        self.legend = []
        
        
        #some configurations that might be changable in the future (hardcoded for now)
        self.lineWidth_fitter = 3
        self.lineWidth_data = 2
        self.color_data = 4
        self.myColorPalette = [ ROOT.kBlue, ROOT.kRed, ROOT.kBlack, ROOT.kBlue+7, ROOT.kOrange, ROOT.kCyan+1, ROOT.kGray+2, ROOT.kYellow-7  ]
        #self.myColorPalette = [ ROOT.kBlue, ROOT.kRed, ROOT.kBlack, ROOT.kOrange, ROOT.kBlue+7, ROOT.kCyan+1, ROOT.kGray+2, ROOT.kYellow-7  ]
    
        self.mTNames = ["[0., 1.02)", "[1.02, 1.14)", "[1.14, 1.20)", "[1.20, 1.26)", "[1.26, 1.38)", "[1.38, 1.56)", "[1.56, 1.86)", "[1.86, 2.21)"]
        self.multNames = ["0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"]

    ##############################################
    ################ functions ###################
    ##############################################
        
    ### cosmetics ###
    def SetDataLineWidth(self, linewidth):
        self.lineWidth_data = linewidth

    ### usefull stuff ###
    def ReadData(self, path, dataname="CF", name1="CATS_input_22all_apap_mTBin_", name2="_Cent_", name3=""):
        """
        name1, name2, name3: strings in order to match the name of the input file
        """
        self.data = []
        for imult in range(self.nMultOffset, self.nMult+self.nMultOffset):
            self.data.append( [] )
            #print(f"mult {imult}") 
            for imt in range(self.nMtOffset, self.nMt+self.nMtOffset):
                
                #print(f"mt {imt}") 
                
                file = ROOT.TFile.Open(path+name1+str(imt)+name2+str(imult)+name3+".root")
                hist = file.Get(dataname)
                self.data[imult-self.nMultOffset].append( hist.Clone() )
                self.data[imult-self.nMultOffset][imt-self.nMtOffset].SetDirectory(0)
                self.data[imult-self.nMultOffset][imt-self.nMtOffset].SetLineWidth(self.lineWidth_data)
                self.data[imult-self.nMultOffset][imt-self.nMtOffset].SetLineColor(self.color_data)
                self.data[imult-self.nMultOffset][imt-self.nMtOffset].SetTitle(f"Centrality: {self.multNames[imult-self.nMultOffset]}, mT in {self.mTNames[imt-self.nMtOffset]}")
                
    def ReadMultipleData(self, path, dataname, name1, name2, name3):
        """
        name1, name2, name3: strings in order to match the name of the input file
        """
        self.data_extra.append( [] )

        for imult in range(self.nMultOffset, self.nMult+self.nMultOffset):
            self.data_extra[len(self.data_extra)-1].append( [] )
            #print(f"mult {imult}") 
            for imt in range(self.nMtOffset, self.nMt+self.nMtOffset):

                file = ROOT.TFile.Open(path+name1+str(imt)+name2+str(imult)+name3+".root")
                hist = file.Get(dataname)
                self.data_extra[len(self.data_extra)-1][imult-self.nMultOffset].append( hist.Clone() )
                self.data_extra[len(self.data_extra)-1][imult-self.nMultOffset][imt-self.nMtOffset].SetDirectory(0)
                self.data_extra[len(self.data_extra)-1][imult-self.nMultOffset][imt-self.nMtOffset].SetLineWidth(self.lineWidth_data)
                self.data_extra[len(self.data_extra)-1][imult-self.nMultOffset][imt-self.nMtOffset].SetLineColor(self.myColorPalette[self.counter-1])
                self.data_extra[len(self.data_extra)-1][imult-self.nMultOffset][imt-self.nMtOffset].SetTitle(f"Centrality: {self.multNames[imult-self.nMultOffset]}, mT in {self.mTNames[imt-self.nMtOffset]}")

    def AddFitResults(self, path, name, description, UseIndividualData):
 
        self.namelist.append( name )
        self.descriptionlist.append( description )
           
        self.fitlist.append( [] )
        self.Bllist.append( [] )
        self.data_extra.append( [] )

        self.radii.append( [] )
        self.radiiErr.append( [] )
        self.mTscalingList.append( [] )

        self.counter = len(self.namelist)
        for imult in range(self.nMultOffset, self.nMult+self.nMultOffset):
            self.fitlist[len(self.namelist)-1].append( [] )
            self.Bllist[len(self.namelist)-1].append( [] )
            self.data_extra[len(self.namelist)-1].append( [] )
            
            self.radii[self.counter-1].append( [] )
            self.radiiErr[self.counter-1].append( [] )
            self.mTscalingList[self.counter-1].append( ROOT.TH1F("mT_scaling_"+self.namelist[self.counter-1]+f"_mult_{imult}", "mT_scaling_"+self.namelist[self.counter-1]+f"_mult_{imult}", 8, 0, 8) )
            self.mTscalingList[self.counter-1][imult].SetDirectory(0)

            for imt in range(self.nMtOffset, self.nMt+self.nMtOffset):
                
                file = file = ROOT.TFile.Open(path+f"FitResults_mTBin_{imt}_Cent_{imult}_smearing_0_femtorange_0_fitrange_0_lambda_0_bl_0.root")
                #Fit result
                fit = file.Get("fFitResult")
                self.fitlist[self.counter-1][imult-self.nMultOffset].append( fit.Clone() )
                self.fitlist[self.counter-1][imult-self.nMultOffset][imt-self.nMtOffset].SetName("fFitResult_"+name[self.counter-1])
                self.fitlist[self.counter-1][imult-self.nMultOffset][imt-self.nMtOffset].SetLineColor(self.myColorPalette[self.counter-1])
                self.fitlist[self.counter-1][imult-self.nMultOffset][imt-self.nMtOffset].SetLineWidth(self.lineWidth_fitter)

                self.radii[self.counter-1][imult-self.nMultOffset].append( fit.GetParameter(6) )
                self.radiiErr[self.counter-1][imult-self.nMultOffset].append( fit.GetParError(6) )
                self.mTscalingList[self.counter-1][imult-self.nMultOffset].SetBinContent(imt-self.nMtOffset, fit.GetParameter(6))
                self.mTscalingList[self.counter-1][imult-self.nMultOffset].SetBinError(imt-self.nMtOffset, fit.GetParError(6))

                #Baseline
                bl = file.Get("fBl")
                self.Bllist[self.counter-1][imult-self.nMultOffset].append( bl.Clone() )
                self.Bllist[self.counter-1][imult-self.nMultOffset][imt-self.nMtOffset].SetName("fBl_"+name[self.counter-1])
                self.Bllist[self.counter-1][imult-self.nMultOffset][imt-self.nMtOffset].SetLineColor(self.myColorPalette[self.counter-1])
                self.Bllist[self.counter-1][imult-self.nMultOffset][imt-self.nMtOffset].SetLineWidth(self.lineWidth_fitter)

                # data
                if UseIndividualData:
                    data = file.Get("hData")
                    self.data_extra[self.counter-1][imult-self.nMultOffset].append( data.Clone() )
                    self.data_extra[self.counter-1][imult-self.nMultOffset][imt-self.nMtOffset].SetName("hData_"+name[self.counter-1])
                    self.data_extra[self.counter-1][imult-self.nMultOffset][imt-self.nMtOffset].SetLineColor(self.myColorPalette[self.counter-1])
                    self.data_extra[self.counter-1][imult-self.nMultOffset][imt-self.nMtOffset].SetLineWidth(self.lineWidth_data)

            
            self.mTscalingList[self.counter-1][imult-self.nMultOffset].SetLineColor(self.myColorPalette[self.counter-1])
            self.mTscalingList[self.counter-1][imult-self.nMultOffset].SetLineWidth(self.lineWidth_fitter)
    
    def MakeLegend(self, lx, ly, ux, uy):
        
        self.legend = ROOT.TLegend(lx, ly, ux, uy)
        #self.legend.AddEntry(self.data[0], "CF")
        for i in range(self.counter):
            self.legend.AddEntry(self.fitlist[i][0][0], self.namelist[i])
    
    def GetAllData(self):
        return self.data
    
    def GetAllFitResults(self):
        return self.fitlist
    
    def GetAllBl(self):
        return self.Bllist
    
    def DrawAll(self, whichmult, whichmt, xrange=None, yrange=None, option="", UseIndividualData=False):
        
        if UseIndividualData:
            if xrange is not None:
                self.data_extra[0][whichmult][whichmt].GetXaxis().SetRangeUser(xrange[0], xrange[1])
            if yrange is not None:
                self.data_extra[0][whichmult][whichmt].GetYaxis().SetRangeUser(yrange[0], yrange[1])
            
            self.data_extra[0][whichmult][whichmt].SetStats(0)
            self.data_extra[0][whichmult][whichmt].Draw(option)

        else:
            if xrange is not None:
                self.data[whichmult][whichmt].GetXaxis().SetRangeUser(xrange[0], xrange[1])
            if yrange is not None:
                self.data[whichmult][whichmt].GetYaxis().SetRangeUser(yrange[0], yrange[1])
            
            self.data[whichmult][whichmt].SetStats(0)
            self.data[whichmult][whichmt].Draw(option)
        
        for i in range(self.counter):

            if UseIndividualData and i>0:
                self.data_extra[i][whichmult][whichmt].Draw(option+" same") 
            self.fitlist[i][whichmult][whichmt].Draw(option+" same")
            self.Bllist[i][whichmult][whichmt].Draw(option+" same")
        
        self.legend.Draw("same")

    def DrawMtScaling(self, whichmult, title="dummy", xrange=None, yrange=None):

        hdummy = ROOT.TH1F(title, title, 100, 0, 100)
        if xrange is not None:
            hdummy.GetXaxis().SetRangeUser(xrange[0], xrange[1])
        if yrange is not None:
            hdummy.GetYaxis().SetRangeUser(yrange[0], yrange[1])

        hdummy.SetStats(0)
        hdummy.DrawClone()
        
        for i in range(self.counter):
            self.mTscalingList[i][whichmult].Draw("same")
        
        #for i in range(self.counter):
        #    if i==0: 
        #        self.mTscalingList[i][whichmult].SetStats(0)
        #        self.mTscalingList[i][whichmult].Draw()
        #    else:
        #        self.mTscalingList[i][whichmult].Draw("same")





    ### compare data only ###

    def AddData(self, data_array, name, internal_name=""):

        self.namelist.append( name )
        self.counter = len(self.namelist)
        
        self.data.append([])
        for i in range(self.nMult):

            self.data[self.counter-1].append( [] )
            for j in range(self.nMt):
                self.data[self.counter-1][i].append( data_array[i][j].Clone() )
                self.data[self.counter-1][i][j].SetDirectory(0)
                self.data[self.counter-1][i][j].SetName(self.compname+"_"+internal_name+"_"+name)
                self.data[self.counter-1][i][j].SetLineColor(self.myColorPalette[self.counter-1])
                self.data[self.counter-1][i][j].SetLineWidth(self.lineWidth_data)

    def MakeLegendAllData(self, lx, ly, ux, uy):
        
        self.legend = ROOT.TLegend(lx, ly, ux, uy)
        #self.legend.AddEntry(self.data[0], "CF")
        for i in range(self.counter):
            self.legend.AddEntry(self.data[i][0][0], self.namelist[i])

    def DrawAllData(self, whichmult, whichpt, xrange=None, yrange=None):
        
         
        full_integral = 0.
        integrals = []
        for i in range(0,self.counter):
            self.data[i][whichmult][whichpt].GetXaxis().SetRangeUser(-0.1, 0.1)
            
            integrals.append( self.data[i][whichmult][whichpt].Integral() ) 
            full_integral += self.data[i][whichmult][whichpt].Integral()
            self.data[i][whichmult][whichpt].GetXaxis().SetRange(0,0) # reset range
        
        
        if xrange is not None:
            self.data[0][whichmult][whichpt].GetXaxis().SetRangeUser(xrange[0], xrange[1])
        if yrange is not None:
            self.data[0][whichmult][whichpt].GetYaxis().SetRangeUser(yrange[0], yrange[1])
            
        self.data[0][whichmult][whichpt].SetStats(0)
        self.data[0][whichmult][whichpt].Draw()

        print(" ")
        print(f"### pT range {whichpt} -> full integral: {full_integral} ###")
        for i in range(1,self.counter):
            
            self.data[i][whichmult][whichpt].Draw("same")
            print(f"Integral {self.namelist[i]} --> {integrals[i]}")
            print(f"  --> ratio = {integrals[i]/full_integral}")

         
        self.legend.Draw("same")