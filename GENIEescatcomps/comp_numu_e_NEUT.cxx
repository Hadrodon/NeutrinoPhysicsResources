#include<../inclusions.h>
#include<../GENIEescatcomps/inclusions/e_func.h>
#include<../GENIEescatcomps/inclusions/other_func.h>
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////

int comp_numu_e_NEUT(){
  /*
  Compares the energy reconstruction of muon neutrino interactions at multiple energies and selections
  Currently holds all the comparisons done with GENIE files, to which a similar file holds the same comparisons but for NEUT files
  */

  hush_root();
  standardize_gStyle();

  vector<string> mode_sel; //Holds the selections for each event. Change this to change the kind of selection you want in the graphs made
  mode_sel.push_back("((NEUTMode == 1) || (NEUTMode == 2) || (NEUTMode == 11) || (NEUTMode == 12) || (NEUTMode == 13)) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96) && (W < 2000)");
  mode_sel.push_back("((NEUTMode == 1) || (NEUTMode == 2) || (NEUTMode == 11) || (NEUTMode == 12) || (NEUTMode == 13)) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96) && (FSProton_mom > 300) && (W < 2000)");
  mode_sel.push_back("(NEUTMode == 1) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96) && (FSProton_mom > 300) && (NTrueNeutrons == 0) && (W < 2000)");
  mode_sel.push_back("(NEUTMode == 2) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96) && (FSProton_mom > 300) && (W < 2000)");
  mode_sel.push_back("((NEUTMode == 11) || (NEUTMode == 12) || (NEUTMode == 13)) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96) && (FSProton_mom > 300) && (W < 2000)");
  mode_sel.push_back("((NEUTMode == 1) || (NEUTMode == 2) || (NEUTMode == 11) || (NEUTMode == 12) || (NEUTMode == 13)) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96) && (W < 2000)");
  mode_sel.push_back("(NEUTMode == 2) && (NTrueProtons == 2) && (NTrueNeutrons == 0) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96) && (FSProton_mom > 300) && (W < 2000");
  mode_sel.push_back("(NEUTMode == 2) && (NTrueProtons == 1) && (NTrueNeutrons == 1) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96) && (FSProton_mom > 300)  && (W < 2000)");
  mode_sel.push_back("((NEUTMode == 1) || (NEUTMode == 2) || (NEUTMode == 11) || (NEUTMode == 12) || (NEUTMode == 13)) && (FSLepton_ct < 0.96) && (FSProton_mom > 300) && (fabs(bjorken_x_rec - 1) < 0.2) && (W < 2000)");
  mode_sel.push_back("((NEUTMode == 1) || (NEUTMode == 2) || (NEUTMode == 11) || (NEUTMode == 12) || (NEUTMode == 13)) && (FSLepton_ct < 0.96) && (FSProton_mom > 300) && (Q2_true > 5E5) && (W < 2000)");

  vector<string> mode_title;
  mode_title.push_back("NEUT_CC0Pi_ee_cut");
  mode_title.push_back("NEUT_CC0Pi_eep_cut");
  mode_title.push_back("NEUT_CCQE_cut");
  mode_title.push_back("NEUT_2p2h_cut");
  mode_title.push_back("NEUT_RES_cut");
  mode_title.push_back("NEUT_CC0Pi_Prot_mom_cut");
  mode_title.push_back("NEUT_2p2h_pp_cut");
  mode_title.push_back("NEUT_2p2h_np_cut");
  mode_title.push_back("NEUT_bjorken_cut");
  mode_title.push_back("NEUT_q2_cut");

  if(mode_sel.size() != mode_title.size()){
    cout << "Check cut selection and title to make sure they match" << endl;
    return 0;
  }

  vector<string> tgt; //Defines the targets in the same order that they are inserted into the input files
  tgt.push_back("Ar");
  tgt.push_back("Ar");
  tgt.push_back("Ar");
  tgt.push_back("Ar");
  tgt.push_back("Fe");
  tgt.push_back("Fe");
  tgt.push_back("Fe");
  tgt.push_back("Fe");
  tgt.push_back("C");
  tgt.push_back("C");

  vector<string> Enu; //Defines the energy of the interactions in the same order that the input files are
  Enu.push_back("1.1");
  Enu.push_back("2.261");
  Enu.push_back("3.3");
  Enu.push_back("4.4");
  Enu.push_back("1.1");
  Enu.push_back("2.261");
  Enu.push_back("3.3");
  Enu.push_back("4.4");
  Enu.push_back("2.261");
  Enu.push_back("4.4");

  vector<string> inter; //Defines the interaction type in the same order that the input files are
  inter.push_back("nu");
  inter.push_back("nu");
  inter.push_back("nu");
  inter.push_back("nu");
  inter.push_back("nu");
  inter.push_back("nu");
  inter.push_back("nu");
  inter.push_back("nu");
  inter.push_back("nu");
  inter.push_back("nu");

  cout << "--- Claiming Trees ---\n";

  vector<TTree *> t; //Holds all the ROOT files we will work with, which files 1-8 are paired with 9-16 for comparisons (this can be changed below in a for loop)
  //1-4 Are for NEUT Ar Enu from 1.1GeV to 4.4GeV
  //5-8 Are for NEUT Ar Escat from 1.1GeV to 4.4GeV
  //9-10 NEUT C12 events for 2.2GeV and 4.4GeV
  t.push_back(get_tree("../CLAS6gens/1.1_Ar_NEUT_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/2.2_Ar_NEUT_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/3.3_Ar_NEUT_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/4.4_Ar_NEUT_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/1.1_Fe_NEUT_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/2.2_Fe_NEUT_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/3.3_Fe_NEUT_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/4.4_Fe_NEUT_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/2.2_C12_NEUT_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/4.4_C12_NEUT_CLAS6Acc.root","CLAS6AcceptSummary"));

  cout << "--- Trees Claimed ---\n";

  TCanvas * c = new TCanvas();
  TLegend * l1 = new TLegend(0.1,0.6,0.6,0.9);
  TLegend * l2 = new TLegend(0.5,0.7,0.9,0.9);

  for(int i = 0; i < 4; ++i){
//  for(int i = 1; i < 2; ++i){ //use for testing purposes
    if(t[i]==0x0){continue;} //skips files that havent been made yet or cant be found
    if(t[i+4]==0x0){continue;}
    cout << "--- Looping through trees, cycle " << i+1 << "/4 ---\n";
    for(int j = 0; j < mode_sel.size(); ++j){
      vector<string> attr1;
      attr1.push_back(mode_sel[j]);
      attr1.push_back(mode_title[j]);
      attr1.push_back(Enu[i]);
      attr1.push_back(tgt[i]);
      attr1.push_back(Form("%s%s%s",tgt[i].c_str(),mode_title[j].c_str(),inter[i].c_str()));
      attr1.push_back(inter[i]);
      vector<string> attr2;
      attr2.push_back(mode_sel[j]);
      attr2.push_back(mode_title[j]);
      attr2.push_back(Enu[i+4]);
      attr2.push_back(tgt[i+4]);
      attr2.push_back(Form("%s%s%s",tgt[i+4].c_str(),mode_title[j].c_str(),inter[i].c_str()));
      attr2.push_back(inter[i+4]);
      e_comp(t[i],t[i+4],attr1,attr2,c,l1,l2);
      //other_comp(t[i],t[i+4],attr1,attr2,c);
      if(i==1){
        cout << "--- Comparing Energy Ranges ---\n";
        comp_energy_ranges(t[4],t[5],t[6],t[7],t[0],t[1],t[2],t[3],c,l1,attr2,attr1);
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////

void make_validification_plots(){
  /*
  Makes plots for validification purposes
  */

  standardize_gStyle();
  TCanvas * c = new TCanvas();

  TTree * t1 = get_tree("../CLAS6gens/0.6_Ar_NEUT_CLAS6Acc.root","CLAS6AcceptSummary");

  cout << t1 << endl;

  char * mode1 = "(NEUTMode == 1) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96)";
  char * mode_t1 = "NEUT_CCQE_eep_cut";
  char * mode2 = "(NEUTMode == 1)";
  char * mode_t2 = "NEUT_CCQE_mode_cut";
  char * mode3 = "(NEUTMode == 1) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96)";
  char * mode_t3 = "NEUT_CCQE_eep_no_q2_cut";

  TH1D * comp1 = new TH1D("temp1","",200,-1,0.6);

  t1->Draw("((EISLepRec_KineFSLep - EISLepTrue)/1000) >> temp1",mode1);

  plot_title_style(comp1,"QE Energy Reconstructed (E_{#nu} = 0.6GeV)","E^{QE}_{reco} - E_{true}","# of Events");

  comp1->Draw();
  TLine * line = new TLine(0,0,0,130);
  line->SetLineColor(kRed);
  line->SetLineStyle(10);
  line->Draw("same");
  c->SaveAs(Form("../GENIEescatcomps/pngdump/nu_0_6_Ar_ereco_%s.png",mode_t1));

  t1->Draw("((EISLepRec_KineFSLep - EISLepTrue)/1000) >> temp1",mode2);

  plot_title_style(comp1,"QE Energy Reconstructed (E_{#nu} = 0.6GeV)","E^{QE}_{reco} - E_{true}","# of Events");

  comp1->Draw();
  TLine * line = new TLine(0,0,0,27600);
  line->SetLineColor(kRed);
  line->SetLineStyle(10);
  line->Draw("same");
  c->SaveAs(Form("../GENIEescatcomps/pngdump/nu_0_6_Ar_ereco_%s.png",mode_t2));

  t1->Draw("((EISLepRec_KineFSLep - EISLepTrue)/1000) >> temp1",mode3);

  plot_title_style(comp1,"QE Energy Reconstructed (E_{#nu} = 0.6GeV)","E^{QE}_{reco} - E_{true}","# of Events");

  comp1->Draw();
  TLine * line = new TLine(0,0,0,18100);
  line->SetLineColor(kRed);
  line->SetLineStyle(10);
  line->Draw("same");
  c->SaveAs(Form("../GENIEescatcomps/pngdump/nu_0_6_Ar_ereco_%s.png",mode_t3));
}

////////////////////////////////////////////////////////////////////////////////////////////////

