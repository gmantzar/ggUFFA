import ROOT
import numpy as np

#the order of the names has to agree with the the order, in which the histogramms are saved
#this is only a temporary solution and has to be done better!
# histnames = ["hSE_rebin", "hME_rebin", "hCF_rebin", "hPt", "hEta", "hPhi", "hDCAxy"]
histnames = ["hSE", "hME", "hCF", "hPt", "hEta", "hPhi", "hDCAxy"]


def getFileContent(file, directory):
    """
    Function to read out the file created by the function FemtoAnalysis.SaveHistogramms(...)
    TODO: the list of histogram names needs to be known as a global variable for now. This is only temporary

    Parameters:

    file (string) = path+filename of the root file from where to read out the histogramms

    directory (string) = directory inside the root file, from which to get the histogramms. This corresponds to
                         the configuration of the substraint, for example "_std/" for the standart train configuration.
                         IMPORTANT: Include the dash in the directory name

    Output:

    A list of all the histogramms inside "directory" of the root file at "file"
    """

    infile = ROOT.TFile(file,"read")

    histos = []
    for i in range(len(histnames)):
        histos.append(infile.Get(directory+histnames[i]))
        histos[i].SetDirectory(0)
        histos[i].SetName(histnames[i])

    return histos

def getSingleHisto(file, directory, histname):
    infile = ROOT.TFile(file,"read")

    histo = infile.Get(directory+histname)
    histo.SetDirectory(0)
    infile.Close()
    return histo

def getNumberOfEvents(file):

    vtxhist = getSingleHisto(file, "_std/", "zvtxhist")
    entries = vtxhist.GetEntries()
    return entries

def normaliseToValue(histo, factor=1., lowerlimit=1, upperlimit=-1):
    if(upperlimit==-1):
        upperlimit = histo.GetNbinsX()

    if type(lowerlimit) is float:
        lowerlimit = histo.FindBin(lowerlimit)
    if type(upperlimit) is float:
        upperlimit = histo.FindBin(upperlimit)

    histo.Scale(factor/histo.Integral(lowerlimit,upperlimit))
    return histo


def normaliseToEvents(histo, factor=1., lowerlimit=1, upperlimit=-1):
    if(upperlimit==-1):
        upperlimit = histo.GetNbinsX()

    if type(lowerlimit) is float:
        lowerlimit = histo.FindBin(lowerlimit)
    if type(upperlimit) is float:
        upperlimit = histo.FindBin(upperlimit)

    histo.Scale(factor/histo.Integral(lowerlimit,upperlimit))
    return histo

def setuphistos(hist, color, width, rebin, lower_x, upper_x, lower_y, upper_y, normarray=None):

    """
    Function to return the histogram "hist" with a specified color, width, rebin and axis ranges.

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

    # print( str(hist.GetBinContent(hist.GetMaximumBin())) )
    # return hist

def setuphistosEventNormalized(hist, color, width, rebin, lower_x, upper_x, lower_y, upper_y, file):
    """
    Function to return the histogram "hist" with a specified color, width, rebin and axis ranges.

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

    nevents = getNumberOfEvents(file)
    # print(str(nevents))
    hist.Scale(1./nevents)
    # print( str(hist.GetBinContent(hist.GetMaximumBin())) )

    if lower_x != None:
         hist.GetXaxis().SetRangeUser(lower_x,upper_x)

    if lower_y != None:
        hist.GetYaxis().SetRangeUser(lower_y,upper_y)

    # return hist

def drawPAPInCanvas(histo_p, histo_ap, title, savename, lx, ly, ux, uy):
    """
    Function to draw the particle and antiparticle versions of a histogram in the same plot and saving the canvas
    (This function is a specialisation of drawInCanvas which is suited only for using on particle and antiparticle version of the same histogram)

    Parameters:

    histo_p (TH1F) = particle histogram

    histo_p (TH1F) = antiparticle histogram

    title (string) = name of the histogram (e.g. ME, SE, etc)

    savename (string) = path+name under which to save the canvas (it contains also the dataformat)

    lx, ly (float) = lower x and y positions of the legend in the canvas

    ux, uy (float) = upper x and y positions of the legend in the canvas

    Output:
    no output. Saves the canvas as explained above
    """

    canvas = ROOT.TCanvas("canvas", "canvas")
    histo_p.Draw()
    histo_ap.Draw("same")
    legSE = ROOT.TLegend(lx, ly, ux, uy)
    legSE.AddEntry(histo_p,title+" particle")
    legSE.AddEntry(histo_ap,title+" antiparticle")
    legSE.Draw()
    canvas.SaveAs(savename)


def CompareParticleAntiparticle(file, path, suffix):
    """
    Function to create an image file for each the SE, ME, CF, pT and eta histograms, where the particle and antiparticle versions
    are compared to each other. It Calls the drawPAPInCanvas function five times (one for each histogram)

    Parameters:

    file (string) = path+name of the file, where the histogramms are located

    path (string) = path, where to save the images

    suffix (string) = name suffix for the image files including the datatype ending (e.g. ".svg")

    Output:
    no output
    """
    histos_particle = getFileContent(file, "_std/")
    histos_antiparticle = getFileContent(file, "_ap-base/")

    setuphistos(histos_particle[0], 4, 2, 2, None, None, None, None)
    setuphistos(histos_antiparticle[0], 2, 2, 2, None, None, None, None)
    drawPAPInCanvas(histos_particle[0],histos_antiparticle[0], "SE", path+"SE"+suffix, 0.6, 0.75, 0.9, 0.9)

    setuphistos(histos_particle[1], 4, 2, 2, None, None, None, None)
    setuphistos(histos_antiparticle[1], 2, 2, 2, None, None, None, None)
    drawPAPInCanvas(histos_particle[1],histos_antiparticle[1], "ME", path+"ME"+suffix, 0.6, 0.75, 0.9, 0.9)

    setuphistos(histos_particle[2], 4, 2, 0, 0., 0.5, 0.7, 4.2)
    setuphistos(histos_antiparticle[2], 2, 2, 0, 0., 0.5, 0.7, 4.2)
    drawPAPInCanvas(histos_particle[2],histos_antiparticle[2], "CF", path+"CF"+suffix, 0.6, 0.75, 0.9, 0.9)

    setuphistos(histos_particle[3], 4, 2, 0, 0., 5., None, None)
    setuphistos(histos_antiparticle[3], 2, 2, 0, 0., 5., None, None)
    drawPAPInCanvas(histos_particle[3],histos_antiparticle[3], "pT", path+"pT"+suffix, 0.6, 0.75, 0.9, 0.9)

    setuphistos(histos_particle[4], 4, 2, 0, None, None, None, None)
    setuphistos(histos_antiparticle[4], 2, 2, 0, None, None, None, None)
    drawPAPInCanvas(histos_particle[4],histos_antiparticle[4], "eta", path+"eta"+suffix, 0.73, 0.75, 0.9, 0.9)


def CompareConfigs(path, filelist, directorylist, confignamelist, colorlist, widthlist, suffix, normarray=[None, None, None, None, None]):
    """
    Function to create an image file for each the SE, ME, CF, pT and eta histograms, where the multiple versions across multiple files
    are compared to each other. It calls the drawInCanvas function five times (one for each histogram)

    Parameters:

    path (string) = path, where to save the images

    filelist (list of string) = list of path+name of the files, where the histogramms are located
                                Note, that for each version of the histogram the filename is needed, even if the same file has been
                                already used for another variation

    directorylist (list of string) = list of directories inside the root file, from which to get the histogramms (one corresponding to each intry in filelist).
                                     This corresponds to the configuration of the substraint, for example "_std/" for the standart train configuration.
                                     IMPORTANT: Include the dash in the directory name

    confignamelist (list of string) = list of names that will be appended to the histogram name in the legend to label the version of the histogramm
                                      (one corresponding to each entry in filelist)

    colorlist (list of int) = list of integers to assign a color to the histogram (one corresponding to each entry in filelist)

    widthlist (list of int) = list of integers to assign a linewidth to the histogram (one corresponding to each entry in filelist)

    normalization (mixed list)

    suffix (string) = name suffix for the image files including the datatype ending (e.g. ".svg")

    Output:
    no output
    """

    # preloop to get the histogramms and save the maximal value
    # for better plotting
    histolist = []
    lowerlimits = np.zeros(len(histnames))
    upperlimits = np.zeros(len(histnames))

    """
    for j in range(len(filelist)):

    print("Lower Limits")
    for i in lowerlimits:
        print(str(i))
    print("Upper Limits")
    for i in upperlimits:
        print(str(i))
    """

    for j in range(len(filelist)):
        histolist.append(getFileContent(filelist[j], directorylist[j]))

        # setuphistos(histolist[j][0],colorlist[j], widthlist[0], 2, 0., 0.2, None, None, normarray[0])   #SE
        setuphistosEventNormalized(histolist[j][0],colorlist[j], widthlist[0], 2, 0., 0.2, None, None, filelist[j])   #SE

        setuphistos(histolist[j][1],colorlist[j], widthlist[1], 2, 0., 0.2, None, None, normarray[1])   #ME

        setuphistos(histolist[j][2],colorlist[j], widthlist[2], 0, 0., 0.2, 0.4, 7., normarray[2])        #CF
        # setuphistos(histolist[j][2],colorlist[j], widthlist[2], 0, 0., 0.2, 0.7, 4.2, normarray[2])        #CF

        # setuphistos(histolist[j][3],colorlist[j], widthlist[3], 0, 0., 5., None, None, normarray[3])       #pT
        setuphistosEventNormalized(histolist[j][3], colorlist[j], widthlist[3], 0., 0., 5., None, None, filelist[j])

        # setuphistos(histolist[j][4],colorlist[j], widthlist[4], 0, None, None, None, None, normarray[4])   #eta
        setuphistosEventNormalized(histolist[j][4],colorlist[j], widthlist[4], 0, None, None, None, None, filelist[j])   #eta


        iterator = 0
        for i in range(len(histolist[j])):

            if histolist[j][i].GetBinContent(histolist[j][i].GetMaximumBin()) > upperlimits[iterator]:
                upperlimits[iterator] =histolist[j][i].GetBinContent(histolist[j][i].GetMaximumBin())

            if histolist[j][i].GetBinContent(histolist[j][i].GetMinimumBin()) < lowerlimits[iterator]:
                lowerlimits[iterator] = histolist[j][i].GetBinContent(histolist[j][i].GetMinimumBin())

            iterator += 1


    drawlist = [0,1,2,3,4]

    #drawInCanvas(histolist, confignamelist, drawlist, path, suffix)

    for i in drawlist:
        canvas = ROOT.TCanvas("canvas", "canvas")
        lx = 0.6
        ly = 0.75
        ux = 0.9
        uy = 0.9

        if i==0 or i==1:
            lx = 0.1
            ux = 0.4

        if i == 4:
            lx = 0.73

        leg = ROOT.TLegend(lx,ly,ux,uy)

        for j in range(len(histolist)):
            histolist[j][i].GetYaxis().SetRangeUser(lowerlimits[i],upperlimits[i]*1.1)
            if j == 0:
                histolist[j][i].Draw("hist")
            else:
                histolist[j][i].Draw("hist same")
            leg.AddEntry(histolist[j][i],histnames[i]+confignamelist[j])
        leg.Draw()
        canvas.SaveAs(path+histnames[i]+suffix)


