#ifndef INC_H
#define INC_H

#include<iostream>
#include<iomanip>
#include<math.h>
#include<cstdlib>
#include<utility>
#include<algorithm>
#include<ctime>
#include<string>
#include<cstring>
#include<vector>

#include"TError.h"
#include"TROOT.h"
#include"TFile.h"
#include"TH1.h"
#include"TH2.h"
#include"TGraph.h"
#include"THStack.h"
#include"TF1.h"
#include"TCanvas.h"
#include"TLegend.h"
#include"TTree.h"
#include"TMath.h"
#include"TRandom.h"
#include"TColor.h"
#include"TObjArray.h"
#include"TLorentzVector.h"
#include"TStyle.h"
#include"THStack.h"
#include"TEfficiency.h"
#include"TPaveLabel.h"
#include"TFitResultPtr.h"
#include"TPaveText.h"
#include"TFitResult.h"
#include"TLatex.h"

using namespace std;

//Below are good colors to use for presentations
int nOrange = TColor::GetColor(230,159,0);
int nSkyBlue = TColor::GetColor(86,180,233);
int nBlueGreen = TColor::GetColor(0,158,115);
int nYellow = TColor::GetColor(240,228,66);
int nBlue = TColor::GetColor(0,114,178);
int nVermillion = TColor::GetColor(213,94,0);
int nRedPurple = TColor::GetColor(204,121,167);

//Used for the get_graph method
//int sel = 0; int shrthnd = 1; int bounds = 2; int path = 3;
//int gen = 0; int inter = 1; int tgt = 2; int Enu = 3;

/////////////////////////////////////////////////////////////////////////////////////////////////


void plot_title_style(TH1 *hist, const char *title="", const char *X_title="", const char *Y_title="");
void set_axis_range(TH1 *hist, float low, float high);
void standardize_gStyle();
void RWB_colors();
void hush_root();
TTree * get_tree(string InFile, string tname);
TH1 * get_graph(TCanvas *c=0x0, TTree *t=0x0, char* graphs[]={}, char* attr[]={}, string cut="", string temp="", bool save_opt=true);

/////////////////////////////////////////////////////////////////////////////////////////////////

void plot_title_style(TH1 * hist, const char * title, const char * X_title, const char * Y_title){
  /* Given a histogram and other pre-required conditions, the macro will change the titles of the
   * histogram and the axes
   */

  hist->SetTitle(title);
  hist->SetXTitle(X_title);
  hist->SetYTitle(Y_title);
  hist->SetTitleSize(0.2);
  hist->GetXaxis()->SetTitleSize(.04);
  hist->GetYaxis()->SetTitleSize(.04);
  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->SetTitleOffset(1.3);

}

/////////////////////////////////////////////////////////////////////////////////////////////////

void set_axis_range(TH1 * hist, float low, float high){
  /* When given a histogram and a high and low point, the macro will set the histogram into a
   * square histogram with bounds being the low and high point
   */

  hist->GetXaxis()->SetRangeUser(low,high);
  hist->GetYaxis()->SetRangeUser(low,high);

}

/////////////////////////////////////////////////////////////////////////////////////////////////

void standardize_gStyle(){
  /* Sets some gStyle standards so that the graph info isnt written for every graph and that
   * the default color of the background is kGray (Gray)
   */

  gStyle->SetOptStat(kFALSE);
  //gStyle->SetCanvasColor(kGray);
  RWB_colors();

}

/////////////////////////////////////////////////////////////////////////////////////////////////

void RWB_colors(){
  /* Changes the color gradient when coloring in object to a different set of colors
   */

  const Int_t NRGBs = 5;//This gives the size of the arrays below.
  const Int_t NCont = 200;//It's the number of colors between which the gradient shifts.
/* Original, but blends too well in with the white background of root
  Double_t stops[NRGBs] = { 0.00, 0.30, 0.65, 0.85, 1.00 };//Locations of colors in % of gradient
  Double_t red[NRGBs]   = { 0.99, 0.99, 0.89, 0.30, 0.00 };//Treat this like a matrix with each
  Double_t green[NRGBs] = { 1.00, 0.92, 0.32, 0.20, 0.00 };//column as the RGB of the color
  Double_t blue[NRGBs]  = { 1.00, 0.39, 0.11, 0.18, 0.00 };//at each stop.
*/

  Double_t stops[NRGBs] = { 0.00, 0.30, 0.65, 0.85, 1.00 };//Locations of colors in % of gradient
  Double_t red[NRGBs]   = { 0.99, 0.99, 0.89, 0.30, 0.00 };//Treat this like a matrix with each
  Double_t green[NRGBs] = { 0.98, 0.92, 0.32, 0.20, 0.00 };//column as the RGB of the color
  Double_t blue[NRGBs]  = { 0.76, 0.39, 0.11, 0.18, 0.00 };//at each stop.
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

}

/////////////////////////////////////////////////////////////////////////////////////////////////

void hush_root(){
  /* Makes it so that ROOT only prints out certain things that depending on what is input to
   * gErrorIgnoreLevel
   */

  gErrorIgnoreLevel = kSysError;
  /* OPTIONS:
   * kPrint - Prints all warnings
   * kInfo - Ignore printing all errors below kInfo
   * kWarning - Ignore printing all errors below kWarning
   * kError - Ignore printing all errors below kError
   * kBreak - Ignore printing all errors below kBreak
   * kSysError - Ignore pritning all errors below kSysError
   * kFatal - Ignore printing all errors below kFatal
   */

}

/////////////////////////////////////////////////////////////////////////////////////////////////

TTree * get_tree(string InFile, string tname){
  /* Gets a TTree from a file given and returns a tree
   */

  TFile * fin = new TFile(InFile.c_str(), "READ");
  TTree * tree;
  fin->GetObject(tname.c_str(), tree);
  if(tree == 0x0){cout << "Could not get tree at: " << InFile << endl; return 0x0;}
  return tree;

}

/////////////////////////////////////////////////////////////////////////////////////////////////

TH1 * get_graph(TCanvas * c, TTree * t, char* graphs[], char* attr[], string cut, string cut_title, string temp, bool save_opt){
  /* Being input all of the prerequisites, the function returns a graph based on what is input
   * and if the input save_opt set true, then it saves it in the inputted graph
   */

  char* sel = graphs[0];
  char* shrthnd = graphs[1];
  char* bounds = graphs[2];
  char* path = graphs[3];
  char* gen = attr[0];
  char* inter = attr[1];
  char* tgt = attr[2];
  char* Enu = attr[3];
  t->Draw(Form("%s >> %s%s",sel,temp.c_str(),bounds),cut.c_str(),"colz");
  TH1 * g = (TH1*)gDirectory->Get(Form("%s",temp.c_str()));
  g->SetLineWidth(2);
  if(save_opt){g->Draw(); c->SaveAs(Form("%s/%s_%s_%s_%s_%s_%s.png",path,gen,inter,Enu,tgt,shrthnd,cut_title.c_str()));}
  return g;

}

/////////////////////////////////////////////////////////////////////////////////////////////////

#endif
