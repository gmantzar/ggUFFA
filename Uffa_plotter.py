def plotMtCorrelations(datapath, firpath, mtbins):
    
    
    file_fit = uproot.open("/Users/georgios/femto/pprun3/analysis/post_processing/pass4/approvals/merged_22moprt/pTCut_2GeV_finalmTBinning/results/Mt_diff_total_fit.root")

    NumMtBins_pp = len(mtbins):
    for iMt in range(NumMtBins_pp):
        WhichMt = 0
        LineColor = ""
        WhichSyst = ""
        SystLatex = ""
        
        LineColor = "mediumblue"
        SystLatex = 'pp'

        WhichMt = iMt
        WhichSyst = "pp"
        
        #extract the data
        file_data = uproot.open(filepath_data+"CATS_input_22moprt_merged_mTBin_"+str(iMt+1)+"_multInt.root")
        hist = file_data["CF"]

        bin_edges = hist.axis().edges()
        bin_contents = hist.values()
        bin_errors = hist.errors()
        
        bin_centers = np.zeros(len(bin_contents))
        for i in range(len(bin_edges)-1):
            bin_centers[i] = bin_edges[i]+0.5*(bin_edges[i+1]-bin_edges[i])
            
        #fit results
        graph_fit = file_fit["TotalFit_mtbin_"+str(iMt+1)]
        x_values_fit = graph_fit.all_members["fX"]
        y_values_fit = graph_fit.all_members["fY"]
        y_errors_fit = graph_fit.all_members["fEY"]

        
        #baseline
        graph_BL = file_fit["Baseline_mtbin_"+str(iMt+1)]
        x_values_BL = graph_BL.all_members["fX"]
        y_values_BL = graph_BL.all_members["fY"]
        y_errors_BL = graph_BL.all_members["fEY"]
        
        
        
        
        
        
        
        
        

        # Define figure and subplots
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True, gridspec_kw={'height_ratios': [2, 1], 'hspace': 0})

        fig.subplots_adjust(top=0.95, right=0.95, bottom=0.09, left=0.11)
        # Upper panel: plot hExp_pp_0, gFit_pp_0 and gBL_pp_0
        if WhichSyst=='pp':
            label_text = 'ALICE data ('+str(SystLatex)+')\n'+r'$m_\mathrm{T}\in['+f'{mt_range_pp[WhichMt]:.2f}'+','+f'{mt_range_pp[WhichMt+1]:.2f}'+')~\mathrm{GeV}$'


        ax1.errorbar(bin_centers, bin_contents, yerr=bin_errors, fmt='o', markersize=6, capsize=3, color='black', markerfacecolor='none', label=label_text)
        #ax1.label_size(4)

        #x = np.linspace(0, 100, 100)
        #y1 = np.sin(x*0.1)+1
        #y2 = np.cos(x*0.1)+1
        #ax1.fill_between(x, y1, y2, color='gray', alpha=0.5, label='band')
        ax1.fill_between(x_values_BL, y_values_BL - y_errors_BL, y_values_BL + y_errors_BL, color='gray', alpha=0.9, label='Baseline')
        ax1.fill_between(x_values_fit, y_values_fit - y_errors_fit, y_values_fit + y_errors_fit, color=LineColor, alpha=0.9, label='Fit')


        # set the thickness of the axis lines
        for spine in ax1.spines.values():
            spine.set_linewidth(1.2)
        for spine in ax2.spines.values():
            spine.set_linewidth(1.2)

        # set the thickness of the tick lines
        ax1.tick_params(axis='both', which='major', width=1.2)
        ax2.tick_params(axis='both', which='major', width=1.2)





        #for i in range(len(x_values_fit)):
        #    ax1.add_patch(Rectangle((x_values_fit[i]-y_errors_fit[i], y_values_fit[i]-0.02), y_errors_fit[i]*2, 0.04, color='red', alpha=0.5))

        #for i in range(len(x_values_BL)):
        #    if i==0:
        #        ax1.add_patch(Rectangle((x_values_BL[i]-y_errors_BL[i], y_values_BL[i]+0.2), y_errors_BL[i]*2, 0.04, color='gray', alpha=0.5))
        #        #ax1.add_patch(fill_between(x_values_BL[i], y_values_BL[i] - y_errors_BL[i], y_values_BL[i] + y_errors_BL[i], color='gray', alpha=0.9))
        #    else:
        #        ax1.add_patch(Rectangle((x_values_BL[i]-y_errors_BL[i], y_values_BL[i]+0.2), y_errors_BL[i]*2, 0.04, color='gray', alpha=0.5))
        #        #ax1.add_patch(fill_betweenx(y_values_BL[i],x_values_BL[i]-4,x_values_BL[i]+4, color='gray', alpha=0.9))

        #plt.fill_between(x_values_fit, y_values_fit - y_errors_fit, y_values_fit + y_errors_fit, color='blue', alpha=0.9)

        if WhichSyst=='pp':
            if WhichMt==NumMtBins_pp-1:
                ax1.set_ylim(0.8, 4.85)
            else:
                ax1.set_ylim(0.8, 4.0)
        else:
            ax1.set_ylim(0.9, 2.2)
        ax1.set_xlim(0.0, 180.0)

        if WhichSyst=="pp":
            ax1_inset = ax1.inset_axes([0.45, 0.3, 0.5, 0.35])#left top width height
            ax1_inset.set_ylim(0.94, 1.13)
            ax1_inset.set_xlim(0.0, 180.0)
            ax1_inset.errorbar(bin_centers, bin_contents, yerr=bin_errors, fmt='o', markersize=6, capsize=3, color='black', markerfacecolor='none')
            ax1_inset.fill_between(x_values_BL, y_values_BL - y_errors_BL, y_values_BL + y_errors_BL, color='gray', alpha=0.9)
            ax1_inset.fill_between(x_values_fit, y_values_fit - y_errors_fit, y_values_fit + y_errors_fit, color=LineColor, alpha=0.9)
            ax1_inset.set_xlabel(r'$k^*$ (MeV/$c$)', fontsize=13)
            ax1_inset.set_ylabel(r'$C(k^*)$', fontsize=13)
            ax1_inset.tick_params(labelsize=12)
            ax1_inset.grid(True)
            ax1_inset.grid(linestyle='--', linewidth='0.5')

        #ax1_inset.fill_betweenx([0, 1], 2, 8, alpha=0.2)
        #for i in range(len(x_values_fit)):
        #    ax1_inset.add_patch(Rectangle((x_values_fit[i]-y_errors_fit[i], y_values_fit[i]-0.02), y_errors_fit[i]*2, 0.04, color='red', alpha=0.5))
        #for i in range(len(x_values_BL)):
        #    ax1_inset.add_patch(Rectangle((x_values_BL[i]-y_errors_BL[i], y_values_BL[i]-0.02), y_errors_BL[i]*2, 0.04, color='gray', alpha=0.5))

        """
        # Lower panel: plot gRat_pp_0
        ax2.set_ylim(-4.5, 4.5)
        #ax2.errorbar(x_values_rat, y_values_rat, yerr=y_errors_rat, fmt='o', color='blue', markersize=4, label='Fit / Data')
        ax2.fill_between(x_values_rat, y_values_rat - y_errors_rat, y_values_rat + y_errors_rat, color=LineColor, alpha=0.9, label="(Fit-Data)/Error")
        ax2.axhline(y=0.0, color='gray', linestyle='--')
        """
        
        # Set axis labels and legends
        ax1.set_ylabel(r'$C(k^*)$', fontsize=16)
        ax1.tick_params(labelsize=14)
        ###ax2.set_xlabel(r'$k^*$ (MeV/$c$)', fontsize=16)
        ###ax2.set_ylabel(r'$n\sigma$', fontsize=16)
        ###ax2.tick_params(labelsize=14)



        ax1.grid(True)
        ax1.grid(linestyle='--', linewidth='0.5')
        ###ax2.grid(True)
        ###ax2.grid(linestyle='--', linewidth='0.5')

        #legend = ax1.legend()
        ## set the inverted handles and labels in the legend
        #legend.set_handles(handles)
        #legend.set_labels(labels)

        #ax1.legend(loc='upper right', fontsize=16)
        ## invert the order of the handles and labels
        handles, labels = ax1.get_legend_handles_labels()
        handles = handles[::-1]
        labels = labels[::-1]
        ax1.legend(handles, labels, fontsize=16)
        ###ax2.legend().remove()

        # Save the figure
        plt.savefig('AI_Plot_'+str(WhichSyst)+'_'+str(WhichMt)+'.pdf')
        plt.savefig('AI_Plot_'+str(WhichSyst)+'_'+str(WhichMt)+'.png')