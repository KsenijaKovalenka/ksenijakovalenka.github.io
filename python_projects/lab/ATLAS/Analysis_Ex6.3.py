'''
# Analysis.py
# Skeleton code in python provided for you
# Ksenija Kovalenka, Matthew Foale 22/03/22
# Make sure to make backups and comment as you go along :-)
####################################################################################################
Draw Histograms using format:
     -    histogram(t, weighting, variable, hist_id, n_bins, xmin, xmax, cuts="1", title=None, x_label = None, y_label = "Number of entries/bin",color = r.kBlack)
        
Define Selection cuts as Boolean String ["&&" = AND, "||" = OR] e.g.:
     -   selection_cut = "("+"sel_cut1"+"&&"+"("+"sel_cut2"+"||"+"sel_cut3"+")"+")"
    is (sel_cut1 AND (sel_cut2 OR sel_cut3))
    
Define new variable to plot using:
     -   SetAlias("name","formula")
     
####################################################################################################
THIS SCRIPT IS RUN BY RunAnalysis.py SO CODE ADJUSTMENTS SHOULD BE MADE HERE.
FOR MORE PLOTTING EXAMPLES SEE ExampleAnalysis.py 
'''
#Allow python to use ROOT
import ROOT as r 
from ShowHistogram import histogram

def Analyse(t,weighting):

#     #EXAMPLE 1: histogram of number of leptons per event
#     histogram(t=t, weighting = weighting,
#               variable = "lep_n", hist_id = "h_lep_n",
#               n_bins=10, xmin=-0.5, xmax=9.5,
#               cuts="1",
#               title = "lep_n:1", x_label = "lep_n", y_label = "lep_n",
#               color = r.kBlack)
        # Plots histogram of number of leptons between range -0.5 to 9.5 with 10 bins
        # No selection cuts are used ("1")
        # Arguments 'cuts','title','x_label','y_label' and 'color' are OPTIONAL arguments
    
    #-------------------------------------------------------------------------------
        
    #EXAMPLE 2: Using Selection Cuts
        # Selection cuts can be used to plot only event that pass certain conditions.
        # Use 'cuts' argument in histogram(...)

        #e.g. Transverse momentum of ALL leptons for events with 3 leptons:
    histogram(t=t, weighting = weighting,
              variable = "lep_pt",  hist_id = "h_lep_pt1",
              n_bins=200, xmin=0, xmax=140e3,
              cuts= "("+"(lep_n == 2)"+"&&"+"(lep_type == 11)"+"&&"+"(lep_charge == -1)"+")")
    
    histogram(t=t, weighting = weighting,
              variable = "lep_pt",  hist_id = "h_lep_pt2",
              n_bins=200, xmin=0, xmax=140e3,
              cuts= "("+"(lep_n == 2)"+"&&"+"(lep_type == 11)"+"&&"+"(lep_charge == 1)"+")")
    
    histogram(t=t, weighting = weighting,
              variable = "lep_pt",  hist_id = "h_lep_pt3",
              n_bins=200, xmin=0, xmax=140e3,
              cuts= "("+"(lep_n == 2)"+"&&"+"(lep_type == 13)"+"&&"+"(lep_charge == -1)"+")")
    
    histogram(t=t, weighting = weighting,
              variable = "lep_pt",  hist_id = "h_lep_pt4",
              n_bins=200, xmin=0, xmax=140e3,
              cuts= "("+"(lep_n == 2)"+"&&"+"(lep_type == 13)"+"&&"+"(lep_charge == 1)"+")")
    
    
#         #e.g. Transverse momentum of a SINGLE lepton (lepton [0]) for events with 3 leptons:
#     histogram(t=t, weighting = weighting,
#               variable = "lep_pt[0]",  hist_id = "h_lep_pt_0",
#               n_bins=200, xmin=0, xmax=140e3,
#               cuts= "(lep_n ==3)")
#         #Note: lep_n ==3 so lep_pt[3],lep_pt[4]...etc will create empty histograms
        
#     histogram(t=t, weighting = weighting,
#               variable = "lep_pt[0]",  hist_id = "h_lep_pt_0",
#               n_bins=200, xmin=0, xmax=140e3,
#               cuts= "(lep_n ==3)")
    
# For a list of predefined variable names see: http://opendata.atlas.cern/release/2020/documentation/datasets/dataset13.html
