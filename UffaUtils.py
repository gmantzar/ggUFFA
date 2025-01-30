import ROOT
import UffaReader as UR
#import numpy as np

import sys
sys.path.append('/Users/georgios/UFFA')
import FileReader as FR


def normaliseToValue(histo, factor=1., lowerlimit=1, upperlimit=-1):
    if(upperlimit==-1):
        upperlimit = histo.GetNbinsX()

    if type(lowerlimit) is float:
        lowerlimit = histo.FindBin(lowerlimit)
    elif lowerlimit == None:
        lowerlimit=1
    if type(upperlimit) is float:
        upperlimit = histo.FindBin(upperlimit)
    elif upperlimit == None:
        upperlimit=-1

    #print("Lowerlimit "+str(lowerlimit)+ " and Upperlimit "+str(upperlimit))
    integral = histo.Integral(lowerlimit, upperlimit)
    histo.Scale(factor/integral)
    #return histo

def setuphistos(hist, color, width, rebin, lower_x, upper_x, lower_y, upper_y, normarray=None):
    """
    Function to return the histogram "hist" with a specified color, width, rebin and axis ranges

    Parameters:

    hist (TH1F) = input histogram that shall be modified

    color (int) = linecolor for the histogram

    width (int) = linewidth of the histogram

    rebin (int) = rebin factor. If given 0 the rebin is not done

    lower_x, upper_x (float) = minimum and maximum value for the displaying range in the x axis
    lower_y, upper_y (float) = minimum and maximum value for the displaying range in the y axis

    Output:
    histogram with the modifications applied according to the parameters

    """

    hist.SetStats(0)
    hist.SetLineColor(color)
    hist.SetLineWidth(width)

    if rebin != 0:
        hist.Rebin(rebin)

    if normarray !=None:
        hist = normaliseToValue(hist, normarray[0], normarray[1], normarray[2])

    if lower_x != None:
         hist.GetXaxis().SetRangeUser(lower_x,upper_x)

    if lower_y != None:
        hist.GetYaxis().SetRangeUser(lower_y,upper_y)

def fetchCorrelationHistos(file, directory, rebin=False, integrated= [False, None]):
    
    histos = []
    isintegrated = integrated[0]
    if(isintegrated):
        ur = UR.UffaReader(file, directory, rebin)

        histos.append(ur.get_se())
        histos.append(ur.get_me())
        histos.append(ur.get_cf())
    else:
        fr = FR.FileReader(file)
        histos.append(fr.get_histo("SE", directory))
        histos.append(fr.get_histo("ME", directory))
        histos.append(fr.get_histo("CF", directory))

    return histos

def setupHistosCosmetics(histo, config, drawrange):
    color = config[0]
    width = config[1]
    markercolor = config[2]
    markerstyle = config[3]
    markersize = config[4]

    histo.SetLineColor(color)
    histo.SetLineWidth(width)
    histo.SetMarkerColor(markercolor)
    histo.SetMarkerStyle(markerstyle)
    histo.SetMarkerSize(markersize)
    xl = drawrange[0]
    xh = drawrange[1]
    yl = drawrange[2]
    yh = drawrange[3]

    if(xl != None and xh != None):
        histo.GetXaxis().SetRangeUser(xl, xh)
    if(yl != None and yh != None):
        histo.GetYaxis().SetRangeUser(yl, yh)

    return histo

def pimpMyHistos(histos, conf):

    """
    function to set cosmetics of histograms
    for an explanation of the list elemets in conf, see below
    """
    color = int(conf[0])
    width = conf[1]
    markerstyle = conf[2]
    markersize = conf[3]
    markercolor = conf[4]

    for hist in histos:
        hist.SetLineColor(color)
        hist.SetLineWidth(width)
        hist.SetMarkerStyle(markerstyle)
        hist.SetMarkerSize(markersize)
        hist.SetMarkerColor(markercolor)
        hist.SetStats(0)
        hist.SetDirectory(0)

def compareUFFA_results(uffa_list, dir_list, rebin_list, name_list):
    """
    uffa_list: list of FileReaders
    dir_list:
    rebin_list:
    """
    pt = []
    eta = []
    phi = []
    dca = []
    dcapt = []


    se = []
    me = []
    me_unw = []
    cf = []
    cf_unw = []

    for i in range(len(uffa_list)):
        #fetch the histogramms
        pt[i] = uffa_list.get_histo("hPt", dir_list[i]+"/Tracks_one")
        eta[i] = uffa_list.get_histo("hEta", dir_list[i]+"/Tracks_one")
        phi[i] = uffa_list.get_histo("hPhi", dir_list[i]+"/Tracks_one")
        dcapt[i] = uffa_list.get_histo("hDCAxy", dir_list[i]+"/Tracks_one")
        dca[i] = dcapt[i].ProjectionY("hdca_"+i, 1, -1)

        se[i] = uffa_list.get_histo("SE", dir_list[i]+rebin_list[i])
        me[i] = uffa_list.get_histo("ME", dir_list[i]+rebin_list[i])
        me_unw[i] = uffa_list.get_histo("ME unw", dir_list[i]+rebin_list[i])
        cf[i] = uffa_list.get_histo("CF", dir_list[i]+rebin_list[i])
        cf_unw[i] = uffa_list.get_histo("CF_unw", dir_list[i]+rebin_list[i])

        #rename 
        pt[i].SetName(name_list[i]+"hPt")
        eta[i].SetName(name_list[i]+"hEta")
        phi[i].SetName(name_list[i]+"hPhi")
        dca[i].SetName(name_list[i]+"hdca")
        
        se[i].SetName(name_list[i]+"hSE")
        me[i].SetName(name_list[i]+"hME")
        me_unw[i].SetName(name_list[i]+"hME_unw")
        cf[i].SetName(name_list[i]+"hCF")
        cf_unw[i].SetName(name_list[i]+"hCF_unw")

        histlist = [ pt[i], eta[i], phi[i], dca[i], se[i], me[i], me_unw[i], cf[i], cf_unw[i]]
        pimpMyHistos(histlist, config) 

def setupHistosForComparison(infile, directory, config, histconfig, cfrebin, differential=None):

    """
    differential: UFFA output has differential SE, ME and CF
                  None: no differential plots
                  int:  take the differential bin int
                  list: CF will be a list with the CFs in the bins of the list
    """
    color = config[0]
    width = config[1]
    markercolor = config[2]
    markerstyle = config[3]
    markersize = config[4]
    #normalizeByEvents = config[5]

    uffareader = UR.UffaReader(infile, directory, cfrebin, differential)

    hist_corr = []
    hist_corr.extend(uffareader.fetchCorrelationHistos())

    hist_qa = []
    hist_qa.extend(uffareader.fetchQAHistos())


    nevents = uffareader.getNevents()

    histlist = hist_corr
    histlist.extend(hist_qa)

    for i in range(len(histlist)):
       
        
        
        #normalize by the number of events
        if(histconfig[i][2] == 1):
            histlist[i].Scale(1./nevents)
        elif(histconfig[i][2] == 2):
            histlist[i].Scale
            normaliseToValue(histlist[i], 1., histconfig[i][0], histconfig[i][1])
        elif(histconfig[i][2] ==0):
            print("No normalization applied for "+histconfig[i][3])
        else:
            print("Value unknown: Please choose between 1, 2 or 0")
        
        #x-axis range
        if(histconfig[i][0]!=None and histconfig[i][1]!=None):
            histlist[i].GetXaxis().SetRangeUser(histconfig[i][0],histconfig[i][1])


        #SetTitle
        # histlist[i].SetTitle(histconfig[i][3])
        # cosmetics
        histlist[i].SetLineColor(color)
        histlist[i].SetLineWidth(width)
        histlist[i].SetMarkerStyle(markerstyle)
        histlist[i].SetMarkerSize(markersize)
        histlist[i].SetMarkerColor(markercolor)
        histlist[i].SetStats(0)

        histlist[i].SetDirectory(0)

    return histlist

def getYdrawLimits(histolist, initval=0., sharpcut=[None, None], useDifference=[True, True], xlimits=None, excludeThreshold=[None, None], debug=False):

    """
    initval:            set the initval, which will be used to start the comparison. In case of C(k*) is is advisable to use initval = 1
                        for distributions better use 0

    sharpcut:           if None, then use the smallest/ largest value without any rescaling for the drawing
                        if it is a number, then lowerlimit = lowerlimit*sharpcut[0] and upperlimit = upperlimit*sharpcut[1]

    useDifference:      instead of up/downscaling the upper/ lower limit with the value in sharpcut, use the difference max-min, scale this with sharpcut and add/ substract it from the limits
    
    xlimits:            None or list with 2 entries
                        - if "None", the maximum across the full x-range will be searched, otherwise only within the specified range
                        - if list, the maximum and minimum will be searched between xlimits[0] < x < xlimits[1]
    
    excludeThresholds:  if a histogram has entries between excludeThreshold[0] and excludeThreshold[1], then it will not be used to compute the limits
    """
    lowerlimit = initval
    upperlimit = initval

    """
    for i in range(len(histolist)):
        if(histolist[i].GetBinContent(histolist[i].GetMaximumBin()) > upperlimit):
           upperlimit = histolist[i].GetBinContent(histolist[i].GetMaximumBin())

        if(histolist[i].GetBinContent(histolist[i].GetMinimumBin()) < lowerlimit):
           lowerlimit = histolist[i].GetBinContent(histolist[i].GetMinimumBin())
    """
    
    for hist in histolist:
        #print("New hist")
        #print("lowerlimit: "+str(lowerlimit))
        #print("upperlimit: "+str(upperlimit))
        if xlimits is None:
            # skip the comparison for outliers
            # this might not be working properly. Test this in a clean scenario
            if excludeThreshold[0] is not None and excludeThreshold[0]>hist.GetBinContent(hist.GetMaximumBin()):
                continue
            if excludeThreshold[1] is not None and excludeThreshold[1]<hist.GetBinContent(hist.GetMaximumBin()):
                continue
            
            if(hist.GetBinContent(hist.GetMaximumBin()) > upperlimit):
                upperlimit = hist.GetBinContent(hist.GetMaximumBin())
            if(hist.GetBinContent(hist.GetMinimumBin()) < lowerlimit):
                lowerlimit = hist.GetBinContent(hist.GetMinimumBin())
        else:
            for xbin in range(hist.FindBin(xlimits[0]), hist.FindBin(xlimits[1])):
                val = hist.GetBinContent(xbin)
                if excludeThreshold[0] is not None and excludeThreshold[0]>val:
                    #print("Enter exlcude for val: "+str(val))
                    continue
                if excludeThreshold[1] is not None and excludeThreshold[1]<val:
                    continue

                #print(val)
                if debug:
                    print(val)
                if val > upperlimit:
                    upperlimit = hist.GetBinContent(xbin)
                if val < lowerlimit:
                    lowerlimit = hist.GetBinContent(xbin)
    # up to here the maximum and minimum value of the histogram have been computed
    # now it is time to set the limits for drawing

    difference = upperlimit-lowerlimit
    if (sharpcut[0] is not None):
        if useDifference[0]: 
            lowerlimit = lowerlimit-(difference*sharpcut[0])
        else:
            lowerlimit = lowerlimit*sharpcut[0]
    if (sharpcut[1] is not None):
        if useDifference[1]:
            upperlimit = upperlimit+(difference*sharpcut[1])
        else:
            upperlimit = upperlimit*sharpcut[1]
    return [lowerlimit, upperlimit]

def setLegendPosition(legpos, size=1):
    if(legpos == 0):
        if(size==1):
            lx = 0.5
            ly = 0.65
            ux = 0.9
            uy = 0.9
        else:
            lx = 0.6
            ly = 0.75
            ux = 0.9
            uy = 0.9
    elif(legpos == 1):
        lx = 0.1
        ly = 0.75
        ux = 0.4
        uy = 0.9
    elif(legpos == 2):
        lx = 0.1
        ly = 0.1
        ux = 0.4
        uy = 0.25
    elif(legpos == 3):
        lx = 0.6
        ly = 0.1
        ux = 0.9
        uy = 0.25
    elif(legpos == 4):
        lx = 0.4
        ly = 0.25
        ux = 0.75
        uy = 0.9
    else:
        print("legend position not supported; doing the default: 0")
        lx = 0.6
        ly = 0.75
        ux = 0.9
        uy = 0.9

    return [lx, ly, ux, uy]

def rebinToCustomBinning(histo, bin_edges):
    
    rebHist = ROOT.TH1F(histo.GetName()+"_rebinned", histo.GetName()+"_rebinned", len(bin_edges)-1, bin_edges)
    #todo: inlcude a warning, in case the upper and lower bin edges do not coincide!
    for i in range(len(bin_edges)-1):
        startbin = histo.FindBin(bin_edges[i])
        endbin = histo.FindBin(bin_edges[i+1])
        buffer = 0.
        for jbin in range(startbin, endbin):
            buffer = buffer + histo.GetBinContent(jbin)
        
        rebHist.SetBinContent(i+1, buffer)
    
    integralOriginal = histo.Integral()
    integralRebinned = rebHist.Integral()
    frac = float(integralOriginal)/float(integralRebinned)
    
    if(frac != 1.0):
        print("WARNING: YOU ARE POTENTIALLY LOOSING ENTRIES WHILE REBINNING!!!")
        print("Integral original: "+str(integralOriginal)+"integral rebinned: "+str(integralRebinned))
        print("fraction: "+str(frac))
        return None
    else:
        return rebHist

def drawInSameCanvas(drawlist, drawoption = "", xDrawLimits=None, setautorange=True, initval=0., sharpcut=[None, None], useDifference=[True, True], xlimits=None, excludeThreshold=[None, None], forceYlimits=[None, None]):
    """
    legpos:         position of the legend (0: upper right, 1: upper left, 2: lower left, 3: lower right)

    initval, sharpcut, useDifference, xlimits, excludeThreshold: see documentation of function "getYdrawLimits" (above)

    forceYlimits    Use this to set the limits to a value after calculating the drawing limits automatically. Usecase: When only one of the two
                    upper limits should be calculated and the other should be fixed
    """

    # position of the legend
    #lx, ly, ux, uy = setLegendPosition(legpos)

    #find the largest bins and set the y range appropreately
    if setautorange:
        limits = getYdrawLimits(drawlist, initval=initval, sharpcut=sharpcut, useDifference=useDifference, xlimits=xlimits, excludeThreshold=excludeThreshold)
        #print("Limits:")
        #print(limits)

    if forceYlimits[0] is not None:
        limits[0] = forceYlimits[0]
    if forceYlimits[1] is not None:
        limits[1] = forceYlimits[1]

    for i in range(len(drawlist)):
        if xDrawLimits is not None:
            drawlist[i].GetXaxis().SetRangeUser(xDrawLimits[0], xDrawLimits[1])
        if(i==0):
            if setautorange:
                drawlist[i].GetYaxis().SetRangeUser(limits[0], limits[1])
            drawlist[i].Draw(drawoption)
        else:
            drawlist[i].Draw("same "+drawoption)

        #leg.AddEntry(drawlist[i], [i])

    #leg.Draw("same")
    #if(savecanvas[0]):
    #    canvas.SaveAs(savecanvas[1])

def drawSameCanvas(drawlist, titlelist, savecanvas, legpos = 0, drawoption = ""):
    """
    legpos:     position of the legend (0: upper right, 1: upper left, 2: lower left, 3: lower right)
    """

    # position of the legend
    lx, ly, ux, uy = setLegendPosition(legpos)

    #find the largest bins and set the y range appropreately
    limits = getYdrawLimits(drawlist)

    canvas = ROOT.TCanvas("canvas_"+savecanvas[2], "canvas_"+savecanvas[2])

    leg = ROOT.TLegend(lx, ly, ux, uy)

    for i in range(len(drawlist)):
        if(i==0):
            drawlist[i].GetYaxis().SetRangeUser(limits[0], limits[1])
            drawlist[i].Draw(drawoption)
        else:
            drawlist[i].Draw("same "+drawoption)

        leg.AddEntry(drawlist[i], titlelist[i])

    leg.Draw("same")
    if(savecanvas[0]):
        canvas.SaveAs(savecanvas[1])
    #canvas.SetDirectory(0)
    return canvas

def drawWithSubCanvas(canvasname, dimensions, histlist, histoconfig):

    canvas = ROOT.TCanvas(canvasname, canvasname)
    canvas.Divide(dimensions[0], dimensions[1])

    for j in range(len(histlist[0])):
        
        drawoption = histoconfig[j][5]
        canvas.cd(j+1)
        for i in range(len(histlist)):

            if(i==0):
                histlist[i][j].Draw(drawoption)
            else:
                histlist[i][j].Draw("same "+drawoption)

    return canvas

def setPositionObject(obj, x1=None, x2=None, y1=None, y2=None):

    if x1 is not None:
        obj.SetX1NDC(x1)
    if x2 is not None:
        obj.SetX2NDC(x2)
    if y1 is not None:
        obj.SetY1NDC(y1)
    if y2 is not None:
        obj.SetY2NDC(y2)

def setPositionPalette(hist, x1=None, x2=None, y1=None, y2=None):

    if x1 is not None:
        hist.GetListOfFunctions().FindObject("palette").SetX1NDC(x1)
    if x2 is not None:
        hist.GetListOfFunctions().FindObject("palette").SetX2NDC(x2)
    if y1 is not None:
        hist.GetListOfFunctions().FindObject("palette").SetY1NDC(y1)
    if y2 is not None:
        hist.GetListOfFunctions().FindObject("palette").SetY2NDC(y2)
    

################################################
#                                              #
#              UFFA SETUPPER                   #
#                                              #
################################################
def moveToGeV(histo):

    nbins = histo.GetNbinsX()
    low_edge = histo.GetBinLowEdge(1)
    up_edge = histo.GetBinLowEdge(nbins+1)
    name =  histo.GetName()

    histo.SetName(name+"_old")
    histo_new = ROOT.TH1F(f'{name}', r';k* (MeV/#it{c}); C(k*)', nbins, low_edge*0.001, up_edge*0.001)
    for i in range(0, nbins+2):
        histo_new.SetBinContent(i, histo.GetBinContent(i))
        histo_new.SetBinError(i, histo.GetBinError(i))

    histo_new.SetDirectory(0) 
    return histo_new

def moveToMeV(histo):

    nbins = histo.GetNbinsX()
    low_edge = histo.GetBinLowEdge(1)
    up_edge = histo.GetBinLowEdge(nbins+1)
    name =  histo.GetName()

    histo.SetName(name+"_old")
    histo_new = ROOT.TH1F(f'{name}', r';k* (MeV/#it{c}); C(k*)', nbins, low_edge*1000, up_edge*1000)
    for i in range(0, nbins+2):
        if histo.GetBinContent(i) == 0:
            continue
        histo_new.SetBinContent(i, histo.GetBinContent(i))
        histo_new.SetBinError(i, histo.GetBinError(i))

    histo_new.SetDirectory(0) 
    return histo_new

def moveToMeV_2D(histo):

    nbinsX = histo.GetNbinsX()
    low_edgeX = histo.GetXaxis().GetBinLowEdge(1)
    up_edgeX = histo.GetXaxis().GetBinLowEdge(nbinsX+1)
    
    nbinsY = histo.GetNbinsY()
    low_edgeY = histo.GetYaxis().GetBinLowEdge(1)
    up_edgeY = histo.GetYaxis().GetBinLowEdge(nbinsX+1)

    name = histo.GetName()
    histo.SetName(name+"_old")
    
    histo_new = ROOT.TH2F(f'{name}', r';k* (MeV/#it{c}); k* (MeV/#it{c}); Entries', nbinsX, low_edgeX*1000, up_edgeX*1000, nbinsY, low_edgeY*1000, up_edgeY*1000)
    for ix in range(0, nbinsX):
        for iy in range(0, nbinsY):
            histo_new.SetBinContent(ix, iy, histo.GetBinContent(ix, iy))
            histo_new.SetBinError(ix, iy, histo.GetBinError(ix, iy))
    
    histo_new.SetDirectory(0)

    return histo_new

def getCatsCF(filename, directory):

    theFile = FR.FileReader(filename, directory)

    uffa_cf = theFile.get_histo("CF", directory)
    uffa_cf_mev = moveToMeV(uffa_cf)

    return uffa_cf_mev

def getCatsME(filename, directory):

    theFile = FR.FileReader(filename, directory)

    uffa_me = theFile.get_histo("ME", directory)
    uffa_me_mev = moveToMeV(uffa_me)

    return uffa_me_mev

def getCatsResMatrix(filename, directory):

    theFile = FR.FileReader(filename, directory)

    resmat = theFile.get_histo("kstar_resolution", directory)
    resmat_mev = moveToMeV_2D(resmat)

    return resmat_mev 

def getCatsResMatrixCustom(filename, directory, name, GeV2MeV):

    theFile = FR.FileReader(filename, directory)

    resmat = theFile.GetHisto(name, directory)
    if(GeV2MeV):
        resmat_mev = moveToMeV_2D(resmat)
        return resmat_mev 
    else: 
        resmat.SetDirectory(0)
        return resmat

def saveForCats(histos, outputfile):

    saveoption = ""
    so = outputfile[1]
    if so == 1:
        saveoption = "new"
    elif so == 2:
        saveoption = "recreate"
    elif so == 3:
        saveoption = "update"
    else:
        print("Save option not supported! Please select 1 (new), 2 (recreate) or 3 (update)")
        return


    outfile = ROOT.TFile(outputfile[0],saveoption)
    for hist in histos:
        hist.Write()
    
    outfile.Close()

def catsSetUpper(file_data, file_mc, outputfile, customResMat=False):

    """
    file_data:      array with 3 entries
                    [0]: path+filename of the input data file
                    [1]: directory inside the input file (e.g.: .../ap-base)

    file_mc:        array with 2 entries
                    [0]: path+filename of the input monte carlo file
                    [1]: directory inside the input file (e.g.: .../ap-base)

    outputfile:     array with 2 entries
                    [0]: path+name of the output file 
                    [1]: Save Option: 1 (new), 2 (recreate) or 3 (update)
    """

    cf_mev = getCatsCF(file_data[0], file_data[1])
    me_mev = getCatsME(file_data[0], file_data[1])

    if(customResMat == False):
        pp_smearing = getCatsResMatrix(file_mc[0], file_mc[1])
    else:
        pp_smearing = getCatsResMatrixCustom(file_mc[0], file_mc[1], file_mc[2], False)

    histos = []
    histos.append(cf_mev)
    histos.append(me_mev)
    histos.append(pp_smearing)

    saveForCats(histos, outputfile)



################################################
#                                              #
#             draw CATS output                 #
#                                              #
################################################
def fetchCatsOutput(catsoutput):

    cfile = FR.FileReader(catsoutput)

    histos = []
    histos.append(cfile.get_histos())
    return histos

def drawCatsOutput(catsoutput, dataset, savecanvas):

    cfile = FR.FileReader(catsoutput)

    #histos = cfile.get_histos()
    histos = []
    histos.append(cfile.get_histo("FitResultLinearBl"))
    histos.append(cfile.get_histo("FitResultCubicBl"))
    histos.append(cfile.get_histo("data"))


    histos[0].SetLineColor(4)
    histos[0].SetLineWidth(3)
    histos[1].SetLineColor(2)
    histos[1].SetLineWidth(2)
    histos[2].SetLineColor(1)
    histos[2].SetLineWidth(1)

    titlelist = []
    titlelist.append("Linear Baseline")
    titlelist.append("Cubic Baseline")
    titlelist.append(dataset)
    drawSameCanvas(histos, titlelist, savecanvas)



################################################
#                                              #
#  functions used in jupiter for comparisons   #
#                                              #
################################################



def drawList(hlist, logaxis=[False, False, False], option=""):
    first = True
    for h in hlist:
        if(first):
            h.Draw(option+"")
            first = False
        else:
            h.Draw(option+"same")
        if(logaxis[0]):
            ROOT.gPad.SetLogx()
        if(logaxis[1]):
            ROOT.gPad.SetLogy()
        if(logaxis[2]):
            ROOT.gPad.SetLogz()
            
def drawSplitCanvas(canvas, hlist, rangex, rangey, logaxis=[False, False, False], option=""):

    """
    Draws the histograms in hlist in a canvas, where each entry of hlist is drawn in a the next subcanvas of canvas
    --> the canvas needs to be defined manually before calling this function
    """
    for i in range(len(hlist)):

        if rangex is not None:
            hlist[i].GetXaxis().SetRangeUser(rangex[0], rangex[1])
        if rangey is not None:
            hlist[i].GetYaxis().SetRangeUser(rangey[0], rangey[1])
        canvas.cd(i+1)

        hlist[i].SetTitle(hlist[i].GetName())
        hlist[i].Draw(option)

        if(logaxis[0]):
            ROOT.gPad.SetLogx()
        if(logaxis[1]):
            ROOT.gPad.SetLogy()
        if(logaxis[2]):
            ROOT.gPad.SetLogz()
            
def DrawAllQA_and_CreateCanvas(hist_list, name_list, divide, logaxis, ranges_x, ranges_y, option):
    canvas_list = []
    
    for htype in range(len(hist_list)):
        
        print(htype)
        canvas_list.append(ROOT.TCanvas("c_"+name_list[htype], "c_"+name_list[htype]))
        canvas_list[htype].Divide(divide[0],divide[1])
        
        drawSplitCanvas(canvas_list[htype], hist_list[htype], ranges_x[htype], ranges_y[htype], logaxis, option)
        canvas_list[htype].Draw()
    return canvas_list

def CreateCanvas(nCanvas, divide, name_list, name_suffix, width=200, height=200):
    canvas_list = []

    if type(nCanvas) is not int:
        raise Exception("Histogram list not supported anymore as input!!! please give an integer corresponding to the number of canvases you want to create")
    for htype in range(nCanvas):
        canvas_list.append(ROOT.TCanvas("c_"+name_list[htype]+name_suffix, "c_"+name_list[htype]+name_suffix, width, height))
        canvas_list[htype].Divide(divide[0], divide[1])
    return canvas_list

def DrawAllQA(canvas_list, hist_list, ranges_x, ranges_y, logaxis, option):
    
    for htype in range(len(canvas_list)):
        drawSplitCanvas(canvas_list[htype], hist_list[htype], ranges_x[htype], ranges_y[htype], logaxis, option)   

def fetchHistos_event(uffa_file, name, uffa_dir, conf, subdir = "/Event"):
    
    hzvtx = uffa_file.get_histo("zvtxhist", uffa_dir+subdir)
    if hzvtx == None:
        hzvtx = uffa_file.get_histo("hZvtx", uffa_dir+subdir)

    hV0Mult = uffa_file.get_histo("MultV0M", uffa_dir+subdir)
    if hV0Mult == None:
        hV0Mult = uffa_file.get_histo("hMultV0M", uffa_dir+subdir)
    
    hmultNTr = uffa_file.get_histo("MultNTr", uffa_dir+subdir)
    if hmultNTr == None:
        hmultNTr = uffa_file.get_histo("hMultNTr", uffa_dir+subdir)

    hzvtx.SetName(name+"hzvtx") 
    hV0Mult.SetName(name+"hV0Mult") 
    hmultNTr.SetName(name+"hmultNTr")
    
    nevents = hzvtx.GetEntries()
    eventhistos = [ hzvtx, hV0Mult, hmultNTr ]
    pimpMyHistos(eventhistos, conf)
    eventhistos.append(nevents)

    return eventhistos 




def fetchHistos_tracks(uffa_file, name, uffa_dir, conf, subdir = "/Tracks_one"):

    hpt = uffa_file.get_histo("hPt", uffa_dir+subdir)
    heta = uffa_file.get_histo("hEta", uffa_dir+subdir)
    hphi = uffa_file.get_histo("hPhi", uffa_dir+subdir)
    hDCAxypt = uffa_file.get_histo("hDCAxy", uffa_dir+subdir)
    hDCAxy = hDCAxypt.ProjectionY("hDCA", 1, -1)
    
    hpt.SetName(name+"hPt")
    heta.SetName(name+"hEta")
    hphi.SetName(name+"hPhi")
    hDCAxypt.SetName(name+"hDCAxypt")
    hDCAxy.SetName(name+"hDCAxy")

    trackhistos =  [ hpt, heta, hphi, hDCAxypt, hDCAxy ] 
    pimpMyHistos(trackhistos, conf)
    
    return trackhistos

def fetchHistos_v0(uffa_file, name, uffa_dir, conf, subdir = "/V0_two"):

    hpt = uffa_file.get_histo("hPt", uffa_dir+subdir)
    heta = uffa_file.get_histo("hEta", uffa_dir+subdir)
    hphi = uffa_file.get_histo("hPhi", uffa_dir+subdir)
    hCPApt = uffa_file.get_histo("hCPA", uffa_dir+subdir)
    hCPA = hCPApt.ProjectionY("hCPA", 1, -1)
    
    hpt.SetName(name+"hPt")
    heta.SetName(name+"hEta")
    hphi.SetName(name+"hPhi")
    hCPApt.SetName(name+"hCPApt")
    hCPA.SetName(name+"hCPA")

    v0histos = [ hpt, heta, hphi, hCPApt, hCPA ]

    #positive tracks
    # childPos = fetchHistos_tracks(uffa_file, name+"_ChildPos", uffa_dir, conf, "/V0Child_pos")
    #negative tracks
    # childNeg = fetchHistos_tracks(uffa_file, name+"_ChildNeg", uffa_dir, conf, "/V0Child_neg")

    pimpMyHistos(v0histos, conf)
    # pimpMyHistos(childPos, conf)
    # pimpMyHistos(childNeg, conf)
    
    # return [ v0histos, childPos, childNeg ] 
    return v0histos

def fetchHistos_cf(uffa_file, name, uffa_dir, conf, rebin_dir=""):
    
    se = uffa_file.get_histo("SE", uffa_dir)
    me = uffa_file.get_histo("ME", uffa_dir)
    me_unw = uffa_file.get_histo("ME unw", uffa_dir)
    cf = uffa_file.get_histo("CF", uffa_dir+rebin_dir)
    cf_unw = uffa_file.get_histo("CF unw", uffa_dir+rebin_dir)

    se.SetName(name+"hSE")
    me.SetName(name+"hME")
    me_unw.SetName(name+"hME_unw")
    cf.SetName(name+"hCF")
    cf_unw.SetName(name+"hCF_unw")
    
    cfhistos = [ se, me, me_unw, cf, cf_unw ]
    pimpMyHistos(cfhistos, conf)

    return cfhistos

def fetchHistos_cf_3d(uffa_file, name, uffa_dir, theconf, nmt, nmult, rebin_dir="", reweightQA=True, theconf_unw=None, normalizeME=False):

    se = []
    me = []
    cf = []

    me_unw = []
    cf_unw = []

    se_mult = []
    me_mult = []
    me_mult_unw = []

    for imt in range(nmt):

        se.append( [] )
        me.append( [] )
        cf.append( [] )

        me_unw.append( [] )
        cf_unw.append( [] )

        se_mult.append( uffa_file.GetHisto("SE kmult", uffa_dir+"bin mt: "+str(imt+1)))
        me_mult.append( uffa_file.GetHisto("ME kmult", uffa_dir+"bin mt: "+str(imt+1)))
        me_mult_unw.append( uffa_file.GetHisto("ME kmult unw", uffa_dir+"bin mt: "+str(imt+1)))

        for imult in range(nmult):
            se[imt].append( uffa_file.GetHisto("SE", uffa_dir+"bin mt: "+str(imt+1)+"/bin: "+str(imult+1)+rebin_dir) )
            me[imt].append( uffa_file.GetHisto("ME", uffa_dir+"bin mt: "+str(imt+1)+"/bin: "+str(imult+1)+rebin_dir) )
            cf[imt].append( uffa_file.GetHisto("CF", uffa_dir+"bin mt: "+str(imt+1)+"/bin: "+str(imult+1)+rebin_dir) )
            
            if reweightQA: 
                me_unw[imt].append( uffa_file.GetHisto("ME_unw", uffa_dir+"bin mt: "+str(imt+1)+"/bin: "+str(imult+1)+rebin_dir) )
                cf_unw[imt].append( uffa_file.GetHisto("CF_unw", uffa_dir+"bin mt: "+str(imt+1)+"/bin: "+str(imult+1)+rebin_dir) )

            se[imt][imult].SetName( name+"hSE_mt_"+str(imt)+"_imult_"+str(imult) )
            me[imt][imult].SetName( name+"hME_mt_"+str(imt)+"_imult_"+str(imult) )
            cf[imt][imult].SetName( name+"hCF_mt_"+str(imt)+"_imult_"+str(imult) )
            if reweightQA: 
                me_unw[imt][imult].SetName( name+"hME_unw_mt_"+str(imt)+"_imult_"+str(imult) )
                cf_unw[imt][imult].SetName( name+"hCF_unw_mt_"+str(imt)+"_imult_"+str(imult) )

                if normalizeME:
                    normaliseToValue(me_unw[imt][imult])
                    normaliseToValue(me[imt][imult])

        pimpMyHistos(se[imt], theconf)
        pimpMyHistos(me[imt], theconf)
        pimpMyHistos(cf[imt], theconf)

        if reweightQA:
            if theconf_unw is None:
                theconfunw = theconf
            else:
                theconfunw = theconf_unw
            pimpMyHistos(me_unw[imt], theconfunw)
            pimpMyHistos(cf_unw[imt], theconfunw)

    if reweightQA:
        return [ [se, me, cf, me_unw, cf_unw ], [se_mult, me_mult, me_mult_unw] ]
    else:
        return [ [se, me, cf], [se_mult, me_mult] ]

def projectMult(kmult, multBins, name, normalize=False):
    
    mult = []
    for i in range(len(multBins)-1):
        mult.append( kmult.ProjectionY(name+"_mult_"+str(i+1)) ) 
        mult[i].Reset()
        
        lowerbin = mult[i].FindBin(multBins[i])
        upperbin = mult[i].FindBin(multBins[i+1])
        
        #print("Lower bin: "+str(lowerbin))
        #print("Upper bin: "+str(upperbin))
        for ibin in range(lowerbin, upperbin):
            imult = kmult.ProjectionX("imult",ibin, ibin)
            imultInt = imult.Integral()
            #print(imultInt)
            mult[i].SetBinContent(ibin, imultInt )
            
        if normalize:
            UU.normaliseToValue(mult[i])
        #mult[i].GetXaxis().SetRangeUser(multBins[i], multBins[i+1])
        if lowerbin ==1:
            mult[i].GetXaxis().SetRange(lowerbin, upperbin)
        elif upperbin == mult[i].GetXaxis().GetNbins():
            mult[i].GetXaxis().SetRange(lowerbin-1, upperbin-1)
        else:
            mult[i].GetXaxis().SetRange(lowerbin-1, upperbin)
            
    return mult

def fetchUffaHistos_pp(uffa_list, name_list, dir_list, rebin_list, conf):

    eventhistos = []
    trackhistos = []
    cfhistos = []     

    for i in range(len(uffa_list)):

        print(i)
        eventhistos.append( fetchHistos_event(uffa_list[i], name_list[i], dir_list[i], conf[i]) )
        trackhistos.append( fetchHistos_tracks(uffa_list[i], name_list[i], dir_list[i], conf[i]) )
        cfhistos.append( fetchHistos_cf(uffa_list[i], name_list[i], dir_list[i], conf[i], rebin_list[i]) )

    return [ eventhistos, cfhistos, trackhistos ]

def fetchUffaHistos_pL(uffa_list, name_list, dir_list, rebin_list, conf):

    eventhistos = []
    v0histos = []
    childPos = []
    childNeg = []
    trackhistos = []
    cfhistos = []

    for i in range(len(uffa_list)):

        eventhistos.append( fetchHistos_event(uffa_list[i], name_list[i], dir_list[i], conf[i]) )
        cfhistos.append( fetchHistos_cf(uffa_list[i], name_list[i], dir_list[i], conf[i], rebin_list[i]) )
        trackhistos.append( fetchHistos_tracks(uffa_list[i], name_list[i], dir_list[i], conf[i]) )
        childPos.append( fetchHistos_tracks(uffa_list[i], name_list[i]+"_ChildPos", dir_list[i], conf[i], "/V0Child_pos") )
        childNeg.append( fetchHistos_tracks(uffa_list[i], name_list[i]+"_ChildNeg", dir_list[i], conf[i], "/V0Child_neg") )
        v0histos.append( fetchHistos_v0(uffa_list[i], name_list[i], dir_list[i], conf[i]) )

    return [ eventhistos, cfhistos, trackhistos, v0histos, childPos, childNeg ]


def fetchDebugHistos_tracks(uffa_file, name, uffa_dir, conf, subdir="/Tracks"):
    
    pid_tpc = uffa_file.get_histo("nSigmaTPC_p", uffa_dir+subdir)
    pid_tof = uffa_file.get_histo("nSigmaTOF_p", uffa_dir+subdir)
    pid_comb = uffa_file.get_histo("nSigmaComb_p", uffa_dir+subdir)
    
    pid_tpc.SetName(name+"_nSigmaTPC_p")
    pid_tof.SetName(name+"_nSigmaTOF_p")
    pid_comb.SetName(name+"_nSigmaComb_p")
    
    hpt = uffa_file.get_histo("hPt", uffa_dir+subdir)
    heta = uffa_file.get_histo("hEta", uffa_dir+subdir)
    hphi = uffa_file.get_histo("hPhi", uffa_dir+subdir)
    hpt.SetName(name+"hPt_debug")
    heta.SetName(name+"hEta_debug")
    hphi.SetName(name+"hPhi_debug")
    
    hDCApt = uffa_file.get_histo("hDCA", uffa_dir+subdir)
    hDCAxypt = uffa_file.get_histo("hDCAxy", uffa_dir+subdir)
    hDCAzpt = uffa_file.get_histo("hDCAz", uffa_dir+subdir)
    hDCApt.SetName(name+"hDCApt_debug")
    hDCAxypt.SetName(name+"hDCAxypt_debug")
    hDCAzpt.SetName(name+"hDCAzpt_debug")
    
    hDCA = hDCApt.ProjectionY("hDCA", 1, -1)
    hDCAxy = hDCAxypt.ProjectionY("hDCAxy", 1, -1)
    hDCAz = hDCAzpt.ProjectionY("hDCAz", 1, -1)
    hDCA.SetName(name+"hDCA_debug")
    hDCAxy.SetName(name+"hDCAxy_debug")
    hDCAz.SetName(name+"hDCAz_debug")
    
    hITSclusters = uffa_file.get_histo("hITSclusters", uffa_dir+subdir)
    hITSclustersIB = uffa_file.get_histo("hITSclustersIB", uffa_dir+subdir)
    hTPCcrossedOverFindable = uffa_file.get_histo("hTPCcrossedOverFindable", uffa_dir+subdir)
    hTPCcrossedRows = uffa_file.get_histo("hTPCcrossedRows", uffa_dir+subdir)
    hTPCfindable = uffa_file.get_histo("hTPCfindable", uffa_dir+subdir)
    hTPCfound = uffa_file.get_histo("hTPCfound", uffa_dir+subdir)
    hTPCshared = uffa_file.get_histo("hTPCshared", uffa_dir+subdir)
    
    hITSclusters.SetName(name+"hITSclusters")
    hITSclustersIB.SetName(name+"hITSclustersIB")
    hTPCcrossedOverFindable.SetName(name+"hTPCcrossedOverFindable")
    hTPCcrossedRows.SetName(name+"hTPCcrossedRows")
    hTPCfindable.SetName(name+"hTPCfindable")
    hTPCfound.SetName(name+"hTPCfound")
    hTPCshared.SetName(name+"hTPCshared")

    debughistos = [ hpt, heta, hphi, hDCA, hDCAxy, hDCAz, hITSclusters, hITSclustersIB, hTPCcrossedOverFindable, hTPCcrossedRows, hTPCfindable, hTPCfound, hTPCshared ]
    pimpMyHistos(debughistos, conf)

    return [ [pid_tpc, pid_tof, pid_comb], [hDCApt, hDCAxypt, hDCAzpt], debughistos ]

def fetchDebugHistos_v0s(uffa_file, name, uffa_dir, conf, subdir="/V0"):

    hDCAdaug = []

    hpt = uffa_file.get_histo("hPt", uffa_dir+subdir)
    heta = uffa_file.get_histo("hEta", uffa_dir+subdir)
    hphi = uffa_file.get_histo("hPhi", uffa_dir+subdir)
    hpt.SetName(name+"hPt_debug")
    heta.SetName(name+"hEta_debug")
    hphi.SetName(name+"hPhi_debug")
    
    hCPApt = uffa_file.get_histo("hCPA", uffa_dir+subdir)
    hCPApt.SetName(name+"hCPApt_debug")
    hCPA = hCPApt.ProjectionY("hCPA", 1, -1)
    hCPA.SetName(name+"hCPA_debug")

    hDCAdaug = uffa_file.get_histo("hDaughDCA_debug", uffa_dir+subdir)
    hinvMassLambda = uffa_file.get_histo("hInvMassLambda_debug", uffa_dir+subdir)
    hinvMassAntiLambda = uffa_file.get_histo("hInvMassAntiLambda_debug", uffa_dir+subdir)
    hTransRadius = uffa_file.get_histo("hTransRadius_debug", uffa_dir+subdir)
    hDCAdaug.SetName(name+"hDaughDCA_debug")
    hinvMassLambda.SetName(name+"hinvMassLambda_debug")
    hinvMassAntiLambda.SetName(name+"hinvMassAntiLambda_debug")
    hTransRadius.SetName(name+"hTransRadius_debug")

    debughistosv0 = [ hpt, heta, hphi, hDCAdaug, hinvMassLambda, hinvMassAntiLambda, hTransRadius, hDCAdaug, hCPA ]
    pimpMyHistos(debughistosv0, conf)
    debughistosv0.append(hCPApt)

    return debughistosv0


def fetchDebugHistos_pp(uffa_list, name_list, dir_list, conf):
    
    debughistos = []
    debughistosPID = []
    debughistos2D = []
    debughistos1D = []

    eventhistos = []
    
    for i in range(len(uffa_list)):
        debughistos.append( fetchDebugHistos_tracks( uffa_list[i], name_list[i], dir_list[i], conf[i]))
        debughistosPID.append( debughistos[i][0] )
        debughistos2D.append( debughistos[i][1] )
        debughistos1D.append( debughistos[i][2] )
        
        eventhistos.append( fetchHistos_event(uffa_list[i], name_list[i], dir_list[i], conf[i]) )
    return [ debughistosPID, debughistos2D, debughistos1D, eventhistos ]

def fetchDebugHistos_pL(uffa_list, name_list, dir_list, conf):
    
    eventhistos = []

    v0debughistos = []

    debughistos_childPos = []
    debughistosPID_childPos = []
    debughistos2D_childPos = []
    debughistos1D_childPos = []
    
    debughistos_childNeg = []
    debughistosPID_childNeg = []
    debughistos2D_childNeg = []
    debughistos1D_childNeg = []

    
    for i in range(len(uffa_list)):
        eventhistos.append( fetchHistos_event(uffa_list[i], name_list[i], dir_list[i], conf[i]) )

        v0debughistos.apppend( fetchDebugHistos_v0s( uffa_list[i], name_list[i], dir_list[i], conf[i]) )

        debughistos_childPos.append( fetchDebugHistos_tracks( uffa_list[i], name_list[i], dir_list[i], conf[i], "/V0Child_pos"))

        debughistosPID_childPos.append( debughistos_childPos[i][0] )
        debughistos2D_childPos.append( debughistos_childPos[i][1] )
        debughistos1D_childPos.append( debughistos_childPos[i][2] )
        
        debughistos_childNeg.append( fetchDebugHistos_tracks( uffa_list[i], name_list[i], dir_list[i], conf[i], "/V0Child_neg"))

        debughistosPID_childNeg.append( debughistos_childNeg[i][0] )
        debughistos2D_childNeg.append( debughistos_childNeg[i][1] )
        debughistos1D_childNeg.append( debughistos_childNeg[i][2] )
        
    return [ eventhistos, v0debughistos,debughistosPID_childPos, debughistos2D_childPos, debughistos1D_childPos, debughistosPID_childNeg, debughistos2D_childNeg, debughistos1D_childNeg  ]




"""
def fetchUffaHistos_pL(uffa_list, name_list, dir_list, rebin_list, conf):
"""

"""
def fetchUffaHistos(uffa_list, name_list, dir_list, rebin_list, conf):


    pt = []
    eta = []
    phi = []
    dca = []
    dcapt = []

    se = []
    me = []
    me_unw = []
    cf = []
    cf_unw = []

    hmult = []
    mult = []




    for i in range(len(uffa_list)):
        #print(i)
        #fetch the histogramms
        pt.append(uffa_list[i].get_histo("hPt", dir_list[i]+"/Tracks_one"))
        eta.append(uffa_list[i].get_histo("hEta", dir_list[i]+"/Tracks_one"))
        phi.append(uffa_list[i].get_histo("hPhi", dir_list[i]+"/Tracks_one"))
        #dcapt.append(uffa_list[i].get_histo("hDCAxy", dir_list[i]+"/Tracks_one"))
        #dca.append(dcapt[i].ProjectionY("hdca_"+str(i), 1, -1))

        se.append(uffa_list[i].get_histo("SE", dir_list[i]))
        me.append(uffa_list[i].get_histo("ME", dir_list[i]))
        me_unw.append(uffa_list[i].get_histo("ME unw", dir_list[i]))
        # cf.append(uffa_list[i].get_histo("CF", dir_list[i]))
        # cf_unw.append(uffa_list[i].get_histo("CF unw", dir_list[i]))    
    
        #rebinned histogramms
        # se.append(uffa_list[i].get_histo("SE", dir_list[i]+rebin_list[i]))
        # me.append(uffa_list[i].get_histo("ME", dir_list[i]+rebin_list[i]))
        # me_unw.append(uffa_list[i].get_histo("ME unw", dir_list[i]+rebin_list[i]))
        cf.append(uffa_list[i].get_histo("CF", dir_list[i]+rebin_list[i]))
        cf_unw.append(uffa_list[i].get_histo("CF unw", dir_list[i]+rebin_list[i]))
    

    
        #rename 
        pt[i].SetName(name_list[i]+"hPt")
        eta[i].SetName(name_list[i]+"hEta")
        phi[i].SetName(name_list[i]+"hPhi")
        #dca[i].SetName(name_list[i]+"hdca")
    
        se[i].SetName(name_list[i]+"hSE")
        me[i].SetName(name_list[i]+"hME")
        me_unw[i].SetName(name_list[i]+"hME_unw")
        cf[i].SetName(name_list[i]+"hCF")
        cf_unw[i].SetName(name_list[i]+"hCF_unw")

        #histlist = [ pt[i], eta[i], phi[i], dca[i], se[i], me[i], me_unw[i], cf[i], cf_unw[i]]
        histlist = [ pt[i], eta[i], phi[i], se[i], me[i], me_unw[i], cf[i], cf_unw[i]]
    
        pimpMyHistos(histlist, conf[i])
    
        #get the vtx plots to get the number of events
        hmult.append(uffa_list[i].get_histo("MultNTr", dir_list[i]+"/Event"))
        #mult.append(hmult[i].Integral())
        mult.append(hmult[i].GetEntries())
        
    return [ pt, eta, phi, dca, dcapt, se, me, me_unw, cf,cf_unw, hmult, mult]
"""


#insert here the fetchPIDPlots function
def fetchPIDPlots(pid_list, name_list, dir_list):
    pid_tpc = []
    pid_tof = []
    pid_comb = []
    
    for i in range(len(pid_list)):
        pid_tpc.append(pid_list[i].get_histo("nSigmaTPC_p", dir_list[i]+"/Tracks"))
        pid_tof.append(pid_list[i].get_histo("nSigmaTOF_p", dir_list[i]+"/Tracks"))
        pid_comb.append(pid_list[i].get_histo("nSigmaComb_p", dir_list[i]+"/Tracks"))
        
        pid_tpc[i].SetName(name_list[i]+"_nSigmaTPC_p")
        pid_tof[i].SetName(name_list[i]+"_nSigmaTOF_p")
        pid_comb[i].SetName(name_list[i]+"_nSigmaComb_p")
    
    return [pid_tpc, pid_tof, pid_comb]

def fetchDebugQAPlots(debug_list, name_list, dir_list, conf):
    hpt = []
    heta = []
    hphi = []
    hDCApt = []
    hDCA = []
    hDCAxypt = []
    hDCAxy = []
    hDCAzpt = []
    hDCAz = []
    
    hpidTPCpt_p = []
    hpidTPCpt_K = []
    hpidTOFpt_p = []
    hpidTOFpt_K = []
    hpidCombpt_p = []
    hpidCombpt_K = []
    
    hpidTPC_p = []
    hpidTPC_K = []
    hpidTOF_p = []
    hpidTOF_K = []
    hpidComb_p = []
    hpidComb_K = []


    hmult = []
    
    mult = []
    
    subdir = "/Tracks"

    for i in range(len(debug_list)):
        hpt.append(debug_list[i].get_histo("hPt", dir_list[i]+subdir))
        heta.append(debug_list[i].get_histo("hEta", dir_list[i]+subdir))
        hphi.append(debug_list[i].get_histo("hPhi", dir_list[i]+subdir))
        
        hDCApt.append(debug_list[i].get_histo("hDCA", dir_list[i]+subdir))
        hDCA.append(hDCApt[i].ProjectionY("hDCA", 1, -1))
        
        hDCAxypt.append(debug_list[i].get_histo("hDCAxy", dir_list[i]+subdir))
        hDCAxy.append(hDCAxypt[i].ProjectionY("hDCAxy", 1, -1))
        
        hDCAzpt.append(debug_list[i].get_histo("hDCAz", dir_list[i]+subdir))
        hDCAz.append(hDCApt[i].ProjectionY("hDCAz", 1, -1))
        
        hpidTPCpt_p.append(debug_list[i].get_histo("nSigmaTPC_p", dir_list[i]+subdir))
        hpidTPC_p.append(hpidTPCpt_p[i].ProjectionY("hpidTPC_p", 1, -1))
        hpidTPCpt_K.append(debug_list[i].get_histo("nSigmaTPC_K", dir_list[i]+subdir))
        hpidTPC_K.append(hpidTPCpt_p[i].ProjectionY("hpidTPC_K", 1, -1))
        
        hpidTOFpt_p.append(debug_list[i].get_histo("nSigmaTOF_p", dir_list[i]+subdir))
        hpidTOF_p.append(hpidTOFpt_p[i].ProjectionY("hpidTOF_p", 1, -1))
        hpidTOFpt_K.append(debug_list[i].get_histo("nSigmaTOF_K", dir_list[i]+subdir))
        hpidTOF_K.append(hpidTOFpt_p[i].ProjectionY("hpidTOF_K", 1, -1))

        hpidCombpt_p.append(debug_list[i].get_histo("nSigmaComb_p", dir_list[i]+subdir))
        hpidComb_p.append(hpidTOFpt_p[i].ProjectionY("hpidComb_p", 1, -1))
        hpidCombpt_K.append(debug_list[i].get_histo("nSigmaComb_K", dir_list[i]+subdir))
        hpidComb_K.append(hpidTOFpt_p[i].ProjectionY("hpidComb_K", 1, -1))
        
        
        hmult.append(debug_list[i].get_histo("MultNTr", dir_list[i]+"/Event"))
        
        #rename the histos
        hpt[i].SetName(name_list[i]+"_hPt")
        heta[i].SetName(name_list[i]+"_hEta")
        hphi[i].SetName(name_list[i]+"_hPhi")
        
        hDCApt[i].SetName(name_list[i]+"_hDCApt")
        hDCA[i].SetName(name_list[i]+"_hDCA") 
        
        hDCAxypt[i].SetName(name_list[i]+"_hDCAxypt")
        hDCAxy[i].SetName(name_list[i]+"_hDCAxy")
        
        hDCAzpt[i].SetName(name_list[i]+"_hDCAzpt")
        hDCAz[i].SetName(name_list[i]+"_hDCAz")
        
        hpidTPCpt_p[i].SetName(name_list[i]+"_nSigmaTPCpt_p")
        hpidTPC_p[i].SetName(name_list[i]+"_hpidTPC_p")
        hpidTPCpt_K[i].SetName(name_list[i]+"_nSigmaTPCpt_K")
        hpidTPC_K[i].SetName(name_list[i]+"_hpidTPC_K")
        
        hpidTOFpt_p[i].SetName(name_list[i]+"_nSigmaTOFpt_p")
        hpidTOF_p[i].SetName(name_list[i]+"_hpidTOF_p")
        hpidTOFpt_K[i].SetName(name_list[i]+"_nSigmaTOFpt_K")
        hpidTOF_K[i].SetName(name_list[i]+"_hpidTOF_K")

        hpidCombpt_p[i].SetName(name_list[i]+"_nSigmaCombpt_p")
        hpidComb_p[i].SetName(name_list[i]+"_hpidComb_p")
        hpidCombpt_K[i].SetName(name_list[i]+"_nSigmaTCombpt_K")
        hpidComb_K[i].SetName(name_list[i]+"_hpidComb_K")
        
        hmult[i].SetName(name_list[i]+"_hmult")
        
        mult.append(hmult[i].Integral())
        
        pimpMyHistos([hpt[i], heta[i], hphi[i], hDCA[i], hDCAxy[i], hDCAz[i], hmult[i]], conf[i])
        pimpMyHistos([hpidTPC_p[i], hpidTOF_p[i], hpidComb_p[i], hpidTPC_K[i], hpidTOF_K[i], hpidComb_K[i]], conf[i])

    return [hpt, heta, hphi, [hDCApt, hDCA, hDCAxypt, hDCAxy, hDCAzpt, hDCAz], [[hpidTPCpt_p, hpidTPC_p, hpidTOFpt_p, hpidTOF_p, hpidCombpt_p, hpidComb_p], [hpidTPCpt_K, hpidTPC_K, hpidTOFpt_K, hpidTOF_K, hpidCombpt_K, hpidComb_K]], hmult, mult]

def FillLegend(legend, hist_list, name_list):

    for i in range(len(hist_list)):
        legend.AddEntry(hist_list[i], name_list[i])
        
def PrintStatisticsAR(se, pt, name_list, nevent_list, compared=[False, 0]):
    """
    compared: 1st argument: Enable/ disable a comparison
              2nd argument: to whihch should the statistics be compared
    """
    
    if(compared[0]):
    #print the relative numbers w.r.t. a variation
        ncandidates_ref = pt[compared[1]].GetEntries()
        ratio_event = 0.
        ratio_eff = 0.
        if(se is not None):
            npairs_ref = se[compared[1]].Integral(1, se[compared[1]].GetXaxis().FindBin(0.2)) 
        for i in range(len(pt)):
            ratio_event = nevent_list[i]/nevent_list[compared[1]]

            print(f"Number of events {name_list[i]}: {ratio_event:.2e}")
            ncandidates = pt[i].GetEntries()
            ratio_ncandidates = pt[i].GetEntries()/ncandidates_ref
            print(f"--> Number of candidates: {ratio_ncandidates:.5e}")
            print(f"----> per event: {(ncandidates/nevent_list[i])/(ncandidates_ref/nevent_list[compared[1]]):.5e}")
            if(se is not None):
                npairs = se[i].Integral(1, se[i].GetXaxis().FindBin(0.2))
                ratio_npairs = se[i].Integral(1, se[i].GetXaxis().FindBin(0.2))/npairs_ref
                print(f"Pairs in low k*: {ratio_npairs:.2e}")
                print(f"--> per event: {(npairs/nevent_list[i])/(npairs_ref/nevent_list[compared[1]]):.5e}")
            print(" ")

    else:
    #print the raw numbers of the list of iput files
        for i in range(len(pt)):

            print(f"Number of events {name_list[i]}: {nevent_list[i]:.2e}")
            ncandidates = pt[i].GetEntries()
            print(f"--> Number of candidates: {ncandidates:.2e}")
            print(f"----> per event: {ncandidates/nevent_list[i]:.5e}")
            if(se is not None):
                npairs = se[i].Integral(1, se[i].GetXaxis().FindBin(0.2))
                print(f"Pairs in low k*: {npairs:.2e}")
                print(f"--> per event: {npairs/nevent_list[i]:.5e}")
            print(" ")




def CalcFractions(input_hist_list, name_list, histname_list, refRun = 0):

    """
    input_hist_list: of format: [type of histogram][list of variations to compute the fraction]
    """
    fraction_list = []
    referenceRun = refRun
        
    for htype in range(len(input_hist_list)):

        fraction_list.append([])
        
        for i in range(len(input_hist_list[htype])):
            fraction_list[htype].append(input_hist_list[htype][i].Clone(name_list[i]+"_"+histname_list[htype]))
          
        #calculate the fraction w.r.t. one selected histogram
        for i in range(len(input_hist_list[htype])):
            if (i==referenceRun):
                continue
            else:
                fraction_list[htype][i].Divide(fraction_list[htype][referenceRun])
        
    return fraction_list
"""

[[
    
    
  [<cppyy.gbl.TH2F object at 0x7f85b8ed8200>, <cppyy.gbl.TH2F object at 0x7f85b8ed7a00>,<cppyy.gbl.TH2F object at 0x7f85b8ed8800>],
  [<cppyy.gbl.TH2F object at 0x7f85ba56a200>,<cppyy.gbl.TH2F object at 0x7f85ba566200>,<cppyy.gbl.TH2F object at 0x7f85ba56b600>],
  [<cppyy.gbl.TH2F object at 0x7f85ba567800>, <cppyy.gbl.TH2F object at 0x7f85ba597200>,<cppyy.gbl.TH2F object at 0x7f85ba568a00>]
  ],
 
 [[<cppyy.gbl.TH2F object at 0x7f85bb63be00>,<cppyy.gbl.TH2F object at 0x7f85bb63e200>,<cppyy.gbl.TH2F object at 0x7f85bb63e800>],
  [<cppyy.gbl.TH2F object at 0x7f85bb5dc400>,<cppyy.gbl.TH2F object at 0x7f85bb5dbc00>,<cppyy.gbl.TH2F object at 0x7f85bb5d7c00>],
  [<cppyy.gbl.TH2F object at 0x7f85a8d76a00>,<cppyy.gbl.TH2F object at 0x7f85a8d76200>,<cppyy.gbl.TH2F object at 0x7f85a8d77c00>]]
  
  
  ]

 """