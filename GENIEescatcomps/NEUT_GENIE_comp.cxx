#include<../inclusions.h>
#include<../GENIEescatcomps/inclusions/e_func.h>
#include<../GENIEescatcomps/inclusions/other_func.h>
#include<vector>
#include<iostream>
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////

int NEUT_GENIE_comp();
TTree * get_tree(string InFile);
TH1 * get_graph(TTree * t, char* sel, string cut, string temp, char* par);

/////////////////////////////////////////////////////////////////////////////////////////////////

int NEUT_GENIE_comp(){
  /*
  Compares NEUT and GENIE files made with similar conditions
  */

  cout << "-- Begining Comparison --\n";
  hush_root();
  gStyle->SetOptStat(kFALSE);
 
  vector<string> sel_ttl;
  sel_ttl.push_back("P_{perp} > 0.2");
  sel_ttl.push_back("total");

  vector<string> Enu; //Contains energy of the input files in order
  Enu.push_back("2.261"); Enu.push_back("4.4"); Enu.push_back("2.261"); Enu.push_back("4.4");

  vector<string> sel; //Selections for the graphs
  sel.push_back("(NEUTMode == 1 || NEUTMode == 2 || NEUTMode == 11 || NEUTMode == 12 || NEUTMode == 13) && Q2_true > 5E5 && (fabs(bjorken_x_rec - 1) < 0.2) && Pperp > 0.2 && NTrueProtons > 0");
  sel.push_back("(NEUTMode == 1 || NEUTMode == 2 || NEUTMode == 11 || NEUTMode == 12 || NEUTMode == 13) && Q2_true > 5E5 && (fabs(bjorken_x_rec - 1) < 0.2) && NTrueProtons > 0");

  cout << "-- Claiming Trees --\n";
  vector<TTree *> t;
  t.push_back(get_tree("../CLAS6gens/2.2_C12_GENIE_CLAS6Acc.root"));
  t.push_back(get_tree("../CLAS6gens/2.2_C12_NEUT_CLAS6Acc.root"));
  t.push_back(get_tree("../CLAS6gens/4.4_C12_GENIE_CLAS6Acc.root"));
  t.push_back(get_tree("../CLAS6gens/4.4_C12_NEUT_CLAS6Acc.root"));

  TFile * CLASdata = new TFile("../GENIEescatcomps/e2a_ep_C12_2261_neutrino6_united2.root","READ");

  TCanvas * c = new TCanvas();
  TLegend * l = new TLegend(0.1,0.7,0.5,0.9);

  for(int i = 0; i < 1; ++i){
    if(t[i]==NULL)continue;
    vector<string> attr;
    attr.push_back("dummy");
    attr.push_back("dummy");
    attr.push_back(Enu[i]);
    attr.push_back("C12");
    attr.push_back("??");
    attr.push_back("nu");
    cout << "-- Starting Comparison with events with  " << Enu[i] << " --\n";
    view_ekin_pperp(t[i],t[i+1],c,l,attr,sel,sel_ttl,CLASdata);
    
  }

  c->Close();
}

/////////////////////////////////////////////////////////////////////////////////////////////////

TTree * get_tree(string InFile){
  /*
  Gets tree from a file ran with the CLAS6AcceptSummary acceptance
  */
  TFile * fin = new TFile(InFile.c_str(), "READ");
  TTree * tree;
  fin->GetObject("CLAS6AcceptSummary", tree);
  if(tree == NULL){cout << "Returning NULL\n"; return NULL;}
  return tree;
}

/////////////////////////////////////////////////////////////////////////////////////////////////

