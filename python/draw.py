from operator import mod
from array import array
import numpy as np
import math
import ROOT
from ROOT import TCanvas, gStyle, TTree, THStack, TH1D
from ROOT import TFile, TH1F, TLegend, TArrow, TGraph
import sys,os

def GetH1D(t):
    nb = t.GetSelectedRows()
    y = t.GetV1()
    x = t.GetV2()
    h = ROOT.TH1D('h', 'h', nb, x[0], x[nb-1])
    ymax = 0
    for i in range(nb):
        h.SetBinContent(i + 1, y[i])
        if y[i] > ymax:
            ymax = y[i]
    return ymax, h

def set_hist_style(h, xtitle, ytitle, color):
    h.GetXaxis().SetNdivisions(509)
    h.GetYaxis().SetNdivisions(504)
    h.SetLineWidth(2)
    h.SetLineWidth(2)
    h.SetStats(0)
    h.SetStats(0)
    h.GetXaxis().SetTitleSize(0.06)
    h.GetXaxis().SetTitleOffset(1.)
    h.GetXaxis().SetLabelOffset(0.01)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetYaxis().SetTitleOffset(1.)
    h.GetYaxis().SetLabelOffset(0.01)
    h.GetXaxis().SetTitle(xtitle)
    h.GetXaxis().CenterTitle()
    h.GetYaxis().SetTitle(ytitle)
    h.GetYaxis().CenterTitle()
    h.SetLineColor(color)

def temp_plt_Vmax_pin():
    simfile = ROOT.TFile("./dataout120V.root")
    datafile = ROOT.TFile("./data/edge_voltage_2019_10_24_15_12_57_HPK-EPI-W2-200-DS-SE5PINNM-01.root")
    c = ROOT.TCanvas("c", "", 100, 100, 800, 600)
    ytitle = 'Vmax/mv'
    xtitle = 'z/um'

    t_data1 = simfile.Get("tree")
    entries_data = t_data1.GetEntries()
    t_data1.Draw('(-v*1.3e2):(z)')
    tmp,h_sim = GetH1D(t_data1)
    set_hist_style(h_sim,xtitle,ytitle,1)

    t_data2 = datafile.Get("edge")
    t_data2.Draw('(-Vmin*1e3):((z*1e3)-11985)',"Vbias==-120")
    tmp,h_data = GetH1D(t_data2)
    set_hist_style(h_data,xtitle,ytitle,2)

    h_sim.Draw()
    h_data.Draw('same')
    legend = ROOT.TLegend(0.14,0.75,0.4,0.85)
    legend.AddEntry(h_sim,"RASER","l")
    legend.AddEntry(h_data,"HPK","l")
    legend.Draw("same")

    c.Update()
    c.SaveAs("pin_vmax_compre.pdf")

def draw_vmax_difV():
    sim1 = ROOT.TFile("./dataout20V.root")
    sim2 = ROOT.TFile("./dataout40V.root")
    sim3 = ROOT.TFile("./dataout60V.root")
    sim4 = ROOT.TFile("./dataout80V.root")
    sim5 = ROOT.TFile("./dataout100V.root")
    sim6 = ROOT.TFile("./dataout120V.root")
    sim7 = ROOT.TFile("./dataout140V.root")
    sim8 = ROOT.TFile("./dataout160V.root")
    sim9 = ROOT.TFile("./dataout180V.root")
    sim10 = ROOT.TFile("./dataout200V.root")

    c = ROOT.TCanvas("c", "", 100, 100, 800, 600)
    ytitle = 'Vmax/mv'
    xtitle = 'z/um'

    t_data1 = sim1.Get("tree")
    t_data1.Draw('(-v*1.3e2):(z)')
    tmp,h_sim1 = GetH1D(t_data1)
    set_hist_style(h_sim1,xtitle,ytitle,1)

    t_data2 = sim2.Get("tree")
    t_data2.Draw('(-v*1.3e2):(z)')
    tmp,h_sim2 = GetH1D(t_data2)
    set_hist_style(h_sim2,xtitle,ytitle,2)

    t_data3 = sim3.Get("tree")
    t_data3.Draw('(-v*1.3e2):(z)')
    tmp,h_sim3 = GetH1D(t_data3)
    set_hist_style(h_sim3,xtitle,ytitle,3)

    t_data4 = sim4.Get("tree")
    t_data4.Draw('(-v*1.3e2):(z)')
    tmp,h_sim4 = GetH1D(t_data4)
    set_hist_style(h_sim4,xtitle,ytitle,4)

    t_data5 = sim5.Get("tree")
    t_data5.Draw('(-v*1.3e2):(z)')
    tmp,h_sim5 = GetH1D(t_data5)
    set_hist_style(h_sim5,xtitle,ytitle,5)

    t_data6 = sim6.Get("tree")
    t_data6.Draw('(-v*1.3e2):(z)')
    tmp,h_sim6 = GetH1D(t_data6)
    set_hist_style(h_sim6,xtitle,ytitle,6)

    t_data7 = sim7.Get("tree")
    t_data7.Draw('(-v*1.3e2):(z)')
    tmp,h_sim7 = GetH1D(t_data7)
    set_hist_style(h_sim7,xtitle,ytitle,7)
    
    t_data8 = sim8.Get("tree")
    t_data8.Draw('(-v*1.3e2):(z)')
    tmp,h_sim8 = GetH1D(t_data8)
    set_hist_style(h_sim8,xtitle,ytitle,8)
    
    t_data9 = sim9.Get("tree")
    t_data9.Draw('(-v*1.3e2):(z)')
    tmp,h_sim9 = GetH1D(t_data9)
    set_hist_style(h_sim9,xtitle,ytitle,9)
    
    t_data10 = sim10.Get("tree")
    t_data10.Draw('(-v*1.3e2):(z)')
    tmp,h_sim10 = GetH1D(t_data10)
    set_hist_style(h_sim10,xtitle,ytitle,11)

    h_sim10.Draw()
    h_sim2.Draw('same')
    h_sim3.Draw('same')
    h_sim4.Draw('same')
    h_sim5.Draw('same')
    h_sim6.Draw('same')
    h_sim7.Draw('same')
    h_sim8.Draw('same')
    h_sim9.Draw('same')
    h_sim1.Draw('same')
    legend = ROOT.TLegend(0.14,0.75,0.4,0.85)
    legend.AddEntry(h_sim1,"RASER_20V","l")
    legend.AddEntry(h_sim2,"RASER_40V","l")
    legend.AddEntry(h_sim3,"RASER_60V","l")
    legend.AddEntry(h_sim4,"RASER_80V","l")
    legend.AddEntry(h_sim5,"RASER_100V","l")
    legend.AddEntry(h_sim6,"RASER_120V","l")
    legend.AddEntry(h_sim7,"RASER_140V","l")
    legend.AddEntry(h_sim8,"RASER_160V","l")
    legend.AddEntry(h_sim9,"RASER_180V","l")
    legend.AddEntry(h_sim10,"RASER_200V","l")
    legend.Draw("same")

    c.Update()
    c.SaveAs("draw_vmax_difV.pdf")

def main():
    args = sys.argv[1:]
    model = args[0]

    if model in ["pin_compre"]:
        temp_plt_Vmax_pin()

    elif model in ["pin_vbias"]:
        draw_vmax_difV()

    else:
        raise NameError(model)

if __name__ == '__main__':
    main() 