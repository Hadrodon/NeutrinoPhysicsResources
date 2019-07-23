/* 
Magdelena Allen | Summer 2018 SULI Program
Macro to read in ROOT output files of simHists_MSA6.C, pick histograms by variable, and plot comparisons in FNAL colors. Adapted from Kevin Ewart (2017). 

Notes:
	- if FSI or nuclear model unspecified in file names, assume defaults (hA and FGMBodekRitchie, respectively). 
	- pion momentum plot deactivated, since excluding pion production in simHists_MSA6.C

Instructions for use: 
	- Set the paths/files to read in from, save path/name, and legend labels for each file on lines 37-53 (2 to 4 files). 
	- Pick to plot a comparison of 2, 3, or 4 files in the main on line 372-274. 
*/

#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include "TROOT.h"
#include "TGaxis.h"

#include <iostream>

using namespace std;


// FNAL color palette
Int_t fnalBlue_i = 999997 //TColor::GetFreeColorIndex();
TColor *fnalBlue = new TColor(fnalBlue_i, 0./255., 76./255., 151./255.); //To change a color: fnalBlue->SetRGB(0., 76./255., 151./255.);
Int_t fnalLightBlue_i = 99998 //TColor::GetFreeColorIndex();
TColor *fnalLightBlue = new TColor(fnalLightBlue_i, 65./255., 182./255., 230./255.);
Int_t fLightOrange_i = 99999 //TColor::GetFreeColorIndex();
TColor *fLightOrange = new TColor(fLightOrange_i, 246./255., 141./255., 46./255.);


// --------------- CONFIGURE INPUT ---------------

// Compare FSI models (hA vs hN vs hA2015 vs hN2015) for MC in nu mode on 12C with 1e6 events: 
//string loadpath1 = "../genie_mc/nu_hA_MC/";
//string file1 = "Histograms_MSA6_nu_12C_hA_MC_n1mil_gst_weighted0.root";
//string label1 = "hA";
//string loadpath2 = "../genie_mc/nu_hA_MC/";
//string file2 = "Histograms_MSA6_nu_12C_hA2015_MC_n1mil_gst_weighted0.root";
//string label2 = "hA2015";
//string loadpath3 = "../genie_mc/nu_hN_MC/";
//string file3 = "Histograms_MSA6_nu_12C_hN_MC_n1mil_gst_weighted0.root";
//string label3 = "hN";
//string loadpath4 = "../genie_mc/nu_hN_MC/";
//string file4 = "Histograms_MSA6_nu_12C_hN2015_MC_n1mil_gst_weighted0.root";
//string label4 = "hN2015";
//string savename = "MSA6_nu_12C_MC_n1mil_hAvshNvshA2015vshN2015_"; //add var_name on later
//string savepath = "../analysis/compare_plots/nu_12C_MC_n1mil_hAvshNvshA2015vshN2015/";


// --------------- CONFIGURE PLOTS ---------------

string vars[] = {
"hvE", "hEQE", "hEQEres", "hcal", "hcalRes", "hq0", "hq3",
"hlE", "hlP", "hlpT", "hlpz", "hlTheta", "hlCosTh", "hlPhi",
"hpP", "hpT", "hppT", "hpTheta", "hpCosTh", "hpPhi", "hpNum",
//"hpiP", 
"hnP", "hknP",
"hdpt", "hdphit", "hdalphat",
"hQ2", "hW", "hbjx", "hnu"
}; //.root histogram vars to take ratio as a function of (automatically ;1)

string var_names[] = {
"Initial lepton energy (GeV)", "Quasi-elastic reconstruction (GeV)", "Quasi-elastic reconstruction resolution", "Calorimetric reconstruction (GeV)", "Calorimetric reconstruction resolution", "Energy transfer q0 (GeV)", "Momentum transfer q3 (GeV)",
"Final lepton energy (GeV)", "Final lepton momentum (GeV)", "Final lepton transverse momentum (Gev)", "Final lepton longitudinal momentum (GeV)", "Lepton scattering angle (degrees)", "Lepton cos(theta)", "Final lepton phi (degrees)",
"Proton momentum (GeV)", "Proton kinetic energy (GeV)", "Proton transverse momentum (GeV)", "Proton angle (degrees)", "Proton cos(theta) (degrees)", "Proton phi (degrees)", "Protons produced",
//"Pion Momentum (GeV) - GENIE", 
"GENIE neutron momentum (GeV)", "Kinematic neutron momentum (GeV)",
"Imbalance momentum (GeV)", "Imbalance phi (degrees)", "Imbalance alpha (degrees)",
"Momentum transfer Q2 (GeV)", "Invariant mass W (GeV)", "Bjorken x", "Kinematic nu (GeV)"
};

string legend_positions[] = {
"low_r", "low_r", "up_r", "low_r", "low_l", "low_r", "low_r",
"low_r", "low_r", "low_r", "low_r", "up_r", "up_l", "low_r",
"up_r", "low_r", "up_r", "up_r", "up_l", "low_r", "up_r",
//"up_l", 
"up_r", "up_r",
"up_r", "up_r", "low_l",
"up_r", "up_r", "low_r", "up_r",
};

int rebins[] = {
2, 2, 4, 2, 2, 4, 4,
4, 1, 1, 1, 1, 4, 1,
2, 4, 2, 4, 4, 1, 1,
//2, 
2, 2,
1, 1, 1,
1, 1, 1, 1,
};

bool norm_bool = true; //turns normalization on or off


// --------------- PLOT COMPARISON ---------------

int plot_compare_2(string var, string var_name, string file_1, string file_2, string label_1, string label_2, string load_path_1, string load_path_2, string save_name, string save_path, bool normed, string leg_pos, int rebin)
/* Overlays a comparison of two distributions on one canvas pad for the given variable. */
{
	//Configure canvas: 1 pad
	TCanvas *c1 = new TCanvas("c1", "Comparison", 1000, 750);
    float YTITLEOFFSET = 0.95;//1.07;
    float XTITLEOFFSET = 0.70;//0.55;ccc
    gStyle->SetStripDecimals(false);
    TGaxis::SetMaxDigits(3);
    //gStyle->SetOptStat(0);
    gStyle->SetPadLeftMargin(0.15);	
	gStyle->SetPadBottomMargin(0.15);
	//Load histograms from file
	TFile *f1 = new TFile((load_path_1+file_1).c_str());
	TFile *f2 = new TFile((load_path_2+file_2).c_str());
	TH1D * h1 = (TH1D*)f1->Get(var.c_str());
	TH1D * h2 = (TH1D*)f2->Get(var.c_str());
	h1->Rebin(rebin);
	h2->Rebin(rebin);
	//Normalization
	string yaxtitle;
	if(normed==true) { //norm by total entries, Poisson weighting
		yaxtitle = "Normalized counts"; 
		h1->Sumw2();
		h2->Sumw2();
		Double_t norm1 = h1->Integral();  //equivalent to h1->GetEntries(); if unweighted
		h1->Scale(1./norm1);
		Double_t norm2 = h2->Integral();
		h2->Scale(1./norm2);
	}else {
		yaxtitle = "Counts";
	}
	//Setup plot
	h1->SetStats(false);
	h1->SetTitle("");
	TLegend *leg;
	if(leg_pos == "low_r") {
		leg = new TLegend(0.58,0.19,0.885,0.33,"");
	}else if(leg_pos == "low_l") {
		leg = new TLegend(0.18,0.19,0.425,0.33,"");
	}else if(leg_pos == "up_r") {
		leg = new TLegend(0.58,0.74,0.885,0.88,"");
	}else if(leg_pos == "up_l") {
		leg = new TLegend(0.18,0.74,0.425,0.88,"");
	}else { 
		cout<<"Legend failure: default to lower right";
		TLegend *leg = new TLegend(0.58,0.19,0.885,0.33,"");
	}
	leg->AddEntry(h1, label_1.c_str());
	leg->AddEntry(h2, label_2.c_str());
	h1->GetXaxis()->SetTitle(var_name.c_str());
	h1->GetXaxis()->SetTitleSize(0.06);
	h1->GetXaxis()->CenterTitle();
	h1->GetXaxis()->SetLabelSize(0.04);
	h1->GetYaxis()->SetTitle(yaxtitle.c_str());
	h1->GetYaxis()->SetTitleSize(0.06);
	h1->GetYaxis()->CenterTitle();
	h1->GetYaxis()->SetLabelSize(0.04);
	float hymax = 1.1 * TMath::Max(h1->GetMaximum(), h2->GetMaximum());
	h1->GetYaxis()->SetRangeUser(0, hymax);
	//Plot histograms with errors		
	h1->SetLineColor(fnalBlue_i); //kBlack
	h1->SetLineWidth(2);
	h1->Draw("hist");
	h1->Draw("e1 same");
	h2->SetLineColor(fLightOrange_i); //kRed
	h2->SetLineWidth(2);
	h2->Draw("hist same");
	h2->Draw("e1 same");
	leg->Draw();
	//Save
    c1->Print((save_path+save_name+var+".png").c_str());
	c1->Close();
	return 0;
}


int plot_compare_3(string var, string var_name, string file_1, string file_2, string file_3, string label_1, string label_2, string label_3, string load_path_1, string load_path_2, string load_path_3, string save_name, string save_path, bool normed, string leg_pos, int rebin)
/* Overlays a comparison of three distributions on one canvas pad for the given variable. */
{
	//Configure canvas: 1 pad
	TCanvas *c1 = new TCanvas("c1", "Comparison", 1000, 750);
    float YTITLEOFFSET = 0.95;//1.07;
    float XTITLEOFFSET = 0.70;//0.55;
    gStyle->SetStripDecimals(false);
    TGaxis::SetMaxDigits(3);
    //gStyle->SetOptStat(0);
    gStyle->SetPadLeftMargin(0.15);	
	gStyle->SetPadBottomMargin(0.15);
	//Load histograms from file
	TFile *f1 = new TFile((load_path_1+file_1).c_str());
	TFile *f2 = new TFile((load_path_2+file_2).c_str());
	TFile *f3 = new TFile((load_path_3+file_3).c_str());
	TH1D * h1 = (TH1D*)f1->Get(var.c_str());
	TH1D * h2 = (TH1D*)f2->Get(var.c_str());
	TH1D * h3 = (TH1D*)f3->Get(var.c_str());
	h1->Rebin(rebin);
	h2->Rebin(rebin);
	h3->Rebin(rebin);
	//Normalization
	string yaxtitle;
	if(normed==true) { //norm by total entries, Poisson weighting
		yaxtitle = "Normalized counts"; 
		h1->Sumw2();
		h2->Sumw2();
		h3->Sumw2();
		Double_t norm1 = h1->Integral();  //equivalent to h1->GetEntries(); if unweighted
		h1->Scale(1./norm1);
		Double_t norm2 = h2->Integral();
		h2->Scale(1./norm2);
		Double_t norm3 = h3->Integral();
		h3->Scale(1./norm3);
	}else {
		yaxtitle = "Counts";
	}
	//Setup plot
	h1->SetStats(false);
	h1->SetTitle("");
	TLegend *leg;
	if(leg_pos == "low_r") {
		leg = new TLegend(0.58,0.19,0.885,0.33,"");
	}else if(leg_pos == "low_l") {
		leg = new TLegend(0.18,0.19,0.425,0.33,"");
	}else if(leg_pos == "up_r") {
		leg = new TLegend(0.58,0.74,0.885,0.88,"");
	}else if(leg_pos == "up_l") {
		leg = new TLegend(0.18,0.74,0.425,0.88,"");
	}else { 
		cout<<"Legend failure: default to lower right";
		TLegend *leg = new TLegend(0.58,0.19,0.885,0.33,"");
	}
	leg->AddEntry(h1, label_1.c_str());
	leg->AddEntry(h2, label_2.c_str());
	leg->AddEntry(h3, label_3.c_str());
	h1->GetXaxis()->SetTitle(var_name.c_str());
	h1->GetXaxis()->SetTitleSize(0.06);
	h1->GetXaxis()->CenterTitle();
	h1->GetXaxis()->SetLabelSize(0.04);
	h1->GetYaxis()->SetTitle(yaxtitle.c_str());
	h1->GetYaxis()->SetTitleSize(0.06);
	h1->GetYaxis()->CenterTitle();
	h1->GetYaxis()->SetLabelSize(0.04);
	float hymax = 1.1 * TMath::Max(h1->GetMaximum(), TMath::Max(h2->GetMaximum(), h3->GetMaximum()));
	h1->GetYaxis()->SetRangeUser(0, hymax);
	//Plot histograms with errors	
	h1->SetLineColor(kBlack);
	h1->SetLineWidth(2);
	h1->Draw("hist");
	h1->Draw("e1 same");
	h2->SetLineColor(fnalBlue_i); //kRed
	h2->SetLineWidth(2);
	h2->Draw("hist same");
	h2->Draw("e1 same");
	h3->SetLineColor(fLightOrange_i); //kBlue
	h3->SetLineWidth(2);
	h3->Draw("hist same");
	h3->Draw("e1 same");
	leg->Draw();
	//Save
    c1->Print((save_path+save_name+var+".png").c_str());
	c1->Close();
	return 0;
}


int plot_compare_4(string var, string var_name, string file_1, string file_2, string file_3, string file_4, string label_1, string label_2, string label_3, string label_4, string load_path_1, string load_path_2, string load_path_3, string load_path_4, string save_name, string save_path, bool normed, string leg_pos, int rebin)
/* Overlays a comparison of four distributions on one canvas pad for the given variable. */
{
	//Configure canvas: 1 pad
	TCanvas *c1 = new TCanvas("c1", "Comparison", 1000, 750);
    float YTITLEOFFSET = 0.95;//1.07;
    float XTITLEOFFSET = 0.70;//0.55;
    gStyle->SetStripDecimals(false);
    TGaxis::SetMaxDigits(3);
    //gStyle->SetOptStat(0);
    gStyle->SetPadLeftMargin(0.15);	
	gStyle->SetPadBottomMargin(0.15);
	//Load histograms from file
	TFile *f1 = new TFile((load_path_1+file_1).c_str());
	TFile *f2 = new TFile((load_path_2+file_2).c_str());
	TFile *f3 = new TFile((load_path_3+file_3).c_str());
	TFile *f4 = new TFile((load_path_4+file_4).c_str());
	TH1D * h1 = (TH1D*)f1->Get(var.c_str());
	TH1D * h2 = (TH1D*)f2->Get(var.c_str());
	TH1D * h3 = (TH1D*)f3->Get(var.c_str());
	TH1D * h4 = (TH1D*)f4->Get(var.c_str());
	h1->Rebin(rebin);
	h2->Rebin(rebin);
	h3->Rebin(rebin);
	h4->Rebin(rebin);
	//Normalization
	string yaxtitle;
	if(normed==true) { //norm by total entries, Poisson weighting
		yaxtitle = "Normalized counts"; 
		h1->Sumw2();
		h2->Sumw2();
		h3->Sumw2();
		h4->Sumw2();
		Double_t norm1 = h1->Integral();  //equivalent to h1->GetEntries(); if unweighted
		h1->Scale(1./norm1);
		Double_t norm2 = h2->Integral();
		h2->Scale(1./norm2);
		Double_t norm3 = h3->Integral();
		h3->Scale(1./norm3);
		Double_t norm4 = h4->Integral();
		h4->Scale(1./norm4);
	}else {
		yaxtitle = "Counts";
	}
	//Setup plot
	h1->SetStats(false);
	h1->SetTitle("");
	TLegend *leg;
	if(leg_pos == "low_r") {
		leg = new TLegend(0.58,0.19,0.885,0.33,"");
	}else if(leg_pos == "low_l") {
		leg = new TLegend(0.18,0.19,0.425,0.33,"");
	}else if(leg_pos == "up_r") {
		leg = new TLegend(0.58,0.74,0.885,0.88,"");
	}else if(leg_pos == "up_l") {
		leg = new TLegend(0.18,0.74,0.425,0.88,"");
	}else { 
		cout<<"Legend failure: default to lower right";
		TLegend *leg = new TLegend(0.58,0.19,0.885,0.33,"");
	}
	leg->AddEntry(h1, label_1.c_str());
	leg->AddEntry(h2, label_2.c_str());
	leg->AddEntry(h3, label_3.c_str());
	leg->AddEntry(h4, label_4.c_str());
	h1->GetXaxis()->SetTitle(var_name.c_str());
	h1->GetXaxis()->SetTitleSize(0.06);
	h1->GetXaxis()->CenterTitle();
	h1->GetXaxis()->SetLabelSize(0.04);
	h1->GetYaxis()->SetTitle(yaxtitle.c_str());
	h1->GetYaxis()->SetTitleSize(0.06);
	h1->GetYaxis()->CenterTitle();
	h1->GetYaxis()->SetLabelSize(0.04);
	float hymax = 1.1 * TMath::Max(h1->GetMaximum(), TMath::Max(h2->GetMaximum(), TMath::Max(h3->GetMaximum(),h4->GetMaximum())));
	h1->GetYaxis()->SetRangeUser(0, hymax);
	//Plot histograms with errors
	h1->SetLineColor(kBlack);
	h1->SetLineWidth(2);
	h1->Draw("hist");
	h1->Draw("e1 same");
	h2->SetLineColor(fnalBlue_i); //kRed
	h2->SetLineWidth(2);
	h2->Draw("hist same");
	h2->Draw("e1 same");
	h3->SetLineColor(fnalLightBlue_i); //kBlue
	h3->SetLineWidth(2);
	h3->Draw("hist same");
	h3->Draw("e1 same");
	h4->SetLineColor(fLightOrange_i); //kGreen
	h4->SetLineWidth(2);
	h4->Draw("hist same");
	h4->Draw("e1 same");
	leg->Draw();
	//Save
    c1->Print((save_path+save_name+var+".png").c_str());
	c1->Close();
	return 0;
}


int ComparePlotter()
{
	gROOT->ForceStyle(); //first canvas still not with forced style (padding)? 
	for(int i=0; i<(sizeof(vars)/sizeof(vars[0])); i++) { //sizeof() gives number of bytes allocated, so ratio of sizeof array to an element gives length
		//plot_compare_2(vars[i], var_names[i], file1, file2, label1, label2, loadpath1, loadpath2, savename, savepath, norm_bool, legend_positions[i], rebins[i]);
		//plot_compare_3(vars[i], var_names[i], file1, file2, file3, label1, label2, label3, loadpath1, loadpath2, loadpath3, savename, savepath, norm_bool, legend_positions[i], rebins[i]);
		plot_compare_4(vars[i], var_names[i], file1, file2, file3, file4, label1, label2, label3, label4, loadpath1, loadpath2, loadpath3, loadpath4, savename, savepath, norm_bool, legend_positions[i], rebins[i]);
	}
	return 0;
}



