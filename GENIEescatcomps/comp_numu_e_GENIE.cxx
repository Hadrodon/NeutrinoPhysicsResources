#include</mnt../inclusions.h>
#include</mnt../GENIEescatcomps/inclusions/e_func.h>
#include</mnt../GENIEescatcomps/inclusions/other_func.h>
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////

int comp_numu_e_GENIE();
TTree * get_tree(string InFile);
void make_validification_plots();
void get_tables();
TH1 * get_g(TTree * t,char * sel,char * cut,char * bounds);

/////////////////////////////////////////////////////////////////////////////////////////////////

int comp_numu_e_GENIE(){
  /*
  Compares the energy reconstruction of muon neutrino interactions at multiple energies and selections
  Currently holds all the comparisons done with GENIE files, to which a similar file holds the same comparisons but for NEUT files
  */
  int sing = 0;
  vector<string> couter = caged();
  hush_root();
  standardize_gStyle();

  vector<string> mode_sel; //Holds the selections for each event. Change this to change the kind of selection you want in the graphs made
  mode_sel.push_back("((NEUTMode == 1) || (NEUTMode == 2) || (NEUTMode == 11) || (NEUTMode == 12) || (NEUTMode == 13)) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96) && (W < 2000)");
  mode_sel.push_back("((NEUTMode == 1) || (NEUTMode == 2) || (NEUTMode == 11) || (NEUTMode == 12) || (NEUTMode == 13)) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96) && (FSProton_mom > 300) && (W < 2000)");
  mode_sel.push_back("(NEUTMode == 1) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96) && (FSProton_mom > 300) && (NTrueNeutrons == 0) && (W < 2000)");
  mode_sel.push_back("(NEUTMode == 2) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96) && (FSProton_mom > 300) && (W < 2000)");
  mode_sel.push_back("((NEUTMode == 11) || (NEUTMode == 12) || (NEUTMode == 13)) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96) && (FSProton_mom > 300) && (W < 2000)");
  mode_sel.push_back("((NEUTMode == 1) || (NEUTMode == 2) || (NEUTMode == 11) || (NEUTMode == 12) || (NEUTMode == 13)) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96) && (W < 2000)");
  mode_sel.push_back("(NEUTMode == 2) && (NTrueProtons == 2) && (NTrueNeutrons == 0) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96) && (FSProton_mom > 300) && NTruePions == 0) && (W < 2000)");
  mode_sel.push_back("(NEUTMode == 2) && (NTrueProtons == 1) && (NTrueNeutrons == 1) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96) && (FSProton_mom > 300) && (W < 2000)");
  mode_sel.push_back("((NEUTMode == 1) || (NEUTMode == 2) || (NEUTMode == 11) || (NEUTMode == 12) || (NEUTMode == 13)) && (FSLepton_ct < 0.96) && (FSProton_mom > 300) && (fabs(bjorken_x_rec - 1) < 0.2) && (W < 2000)");
  mode_sel.push_back("((NEUTMode == 1) || (NEUTMode == 2) || (NEUTMode == 11) || (NEUTMode == 12) || (NEUTMode == 13)) && (FSLepton_ct < 0.96) && (FSProton_mom > 300) && (Q2_true > 5E5) && (W < 2000)");

  vector<string> mode_title;
  mode_title.push_back("GENIE_CC0Pi_ee_cut");
  mode_title.push_back("GENIE_CC0Pi_eep_cut");
  mode_title.push_back("GENIE_CCQE_cut");
  mode_title.push_back("GENIE_2p2h_cut");
  mode_title.push_back("GENIE_RES_cut");
  mode_title.push_back("GENIE_CC0Pi_Prot_mom_cut");
  mode_title.push_back("GENIE_2p2h_pp_cut");
  mode_title.push_back("GENIE_2p2h_np_cut");
  mode_title.push_back("GENIE_bjorken_cut");
  mode_title.push_back("GENIE_q2_cut");

  if(mode_sel.size() != mode_title.size()){
    cout << "Check cut selection and title to make sure they match" << endl;
    return 0;
  }

  vector<string> tgt; //Defines the targets in the same order that they are inserted into the input files
  tgt.push_back("Ar");
  tgt.push_back("Ar");
  tgt.push_back("Ar");
  tgt.push_back("Ar");
  tgt.push_back("Ar");
  tgt.push_back("Ar");
  tgt.push_back("Ar");
  tgt.push_back("Ar");
  tgt.push_back("Fe");
  tgt.push_back("Fe");
  tgt.push_back("Fe");
  tgt.push_back("Fe");
  tgt.push_back("Fe");
  tgt.push_back("Fe");
  tgt.push_back("Fe");
  tgt.push_back("Fe");
  tgt.push_back("Ar");
  tgt.push_back("Ar");
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
  Enu.push_back("1.1");
  Enu.push_back("2.261");
  Enu.push_back("3.3");
  Enu.push_back("4.4");
  Enu.push_back("1.1");
  Enu.push_back("2.261");
  Enu.push_back("3.3");
  Enu.push_back("4.4");
  Enu.push_back("1.1");
  Enu.push_back("2.261");
  Enu.push_back("2.261");
  Enu.push_back("4.4");

  vector<string> inter; //Defines the interaction type in the same order that the input files are
  inter.push_back("nu");
  inter.push_back("nu");
  inter.push_back("nu");
  inter.push_back("nu");
  inter.push_back("scat");
  inter.push_back("scat");
  inter.push_back("scat");
  inter.push_back("scat");
  inter.push_back("nu");
  inter.push_back("nu");
  inter.push_back("nu");
  inter.push_back("nu");
  inter.push_back("scat");
  inter.push_back("scat");
  inter.push_back("scat");
  inter.push_back("scat");
  inter.push_back("nu");
  inter.push_back("scat");
  inter.push_back("nu");
  inter.push_back("nu");

  cout << "--- Claiming Trees ---\n";

  vector<TTree *> t; //Holds all the ROOT files we will work with, which files 1-8 are paired with 9-16 for comparisons (this can be changed below in a for loop)
  //1-4 Are for GENIE Ar Enu from 1.1GeV to 4.4GeV
  //5-8 Are for GENIE Ar Escat from 1.1GeV to 4.4GeV
  //9-12 Are for GENIE Fe Enu from 1.1GeV to 4.4GeV
  //13-16 Are for GENIE Fe Escat from 1.1GeV to 4.4GeV
  //17-18 Are for GibUU events
  //19-20 GENIE C12 events for 2.2GeV and 4.4GeV
  t.push_back(get_tree("../CLAS6gens/1.1_Ar_GENIE_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/2.2_Ar_GENIE_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/3.3_Ar_GENIE_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/4.4_Ar_GENIE_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/1.1_Ar_GENIE_escat_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/2.2_Ar_GENIE_escat_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/3.3_Ar_GENIE_escat_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/4.4_Ar_GENIE_escat_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/1.1_Fe_GENIE_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/2.2_Fe_GENIE_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/3.3_Fe_GENIE_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/4.4_Fe_GENIE_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/1.1_Fe_GENIE_escat_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/2.2_Fe_GENIE_escat_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/3.3_Fe_GENIE_escat_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/4.4_Fe_GENIE_escat_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(0x0); //Null here since it does not exist (which is what 0x0 represents)
  t.push_back(get_tree("../CLAS6gens/GiBUU.escat.Ar40.NoFSI.stdhep.CLAS6.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/2.2_C12_GENIE_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/4.4_C12_GENIE_CLAS6Acc.root","CLAS6AcceptSummary"));

  cout << "--- Trees Claimed ---\n";

  TCanvas * c = new TCanvas();
  TLegend * l1 = new TLegend(0.1,0.6,0.6,0.9);
  TLegend * l2 = new TLegend(0.5,0.7,0.9,0.9);

  for(int i = 0; i < 8; ++i){
//  for(int i = 1; i < 2; ++i){ //use for testing purposes
    if(t[i]==0x0){continue;} //skips files that havent been made yet or cant be found
    if(t[i+8]==0x0){continue;} //skips files that dont have a pair
    cout << "--- Looping through trees, cycle " << i+1 << "/8 ---\n";
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
      attr2.push_back(Enu[i+8]);
      attr2.push_back(tgt[i+8]);
      attr2.push_back(Form("%s%s%s",tgt[i+8].c_str(),mode_title[j].c_str(),inter[i].c_str()));
      attr2.push_back(inter[i+8]);
      e_comp(t[i],t[i+8],attr1,attr2,c,l1,l2);
      //other_comp(t[i],t[i+8],attr1,attr2,c);

      if( i == 1){
        cout << "--- Generating Carbon-12 Graphs ---\n";
        TH1D * carb2 = e_frac(t[18],c,attr1,Form("C2%s",attr1[1].c_str()));
        TH1D * carb4 = e_frac(t[19],c,attr1,Form("C4%s",attr1[1].c_str()));

        if(attr1[1] == "GENIE_CC0Pi_ee'p_cut"){ //outputs the fraction of events in certain energy bands
          cout << Form("\n|            |tgt:C12 |cut: %s|   Energy: 2.261   |\n",attr1[1].c_str());
          cout << "|   Energy   | 0.75-1 | 1-1.25 | 1.25-1.5 | 1.5-1.75 | 1.75-2 |\n";
          cout << Form("| Event Frac | %.4f | %.4f |  %.4f  |  %.4f  | %.4f |\n",carb2->Integral(27,36),carb2->Integral(37,45),carb2->Integral(46,53),carb2->Integral(54,62),carb2->Integral(63,72));
        
          cout << Form("\n|            |tgt:C12 |cut: %s|  Energy: 4.4  |\n",attr1[1].c_str());
          cout << "|   Energy   |  1.5-2 | 2-2.5  | 2.5-3  | 3-3.5  | 3.5-4  |\n";
          cout << Form("| Event Frac | %.4f | %.4f | %.4f | %.4f | %.4f |\n",carb4->Integral(27,36),carb4->Integral(37,45),carb4->Integral(46,53),carb4->Integral(54,62),carb4->Integral(63,72));

        }
        cout << "--- Comparing Energy Ranges ---\n";
        comp_energy_ranges(t[8],t[9],t[10],t[11],t[0],t[1],t[2],t[3],c,l1,attr2,attr1);

        cout << "--- Generating GiBUU Graphs(temporary) ---\n";
        if(j == 0)attr1[1] = "GiBUU_CC0Pi_cut";
        if(j == 1)attr1[1] = "GiBUU_CCQE_cut";
        if(j == 2)attr1[1] = "GiBUU_2p2h_cut";
        if(j == 3)attr1[1] = "GiBUU_RES_cut";
        if(j == 4)attr1[1] = "GiBUU_CC0Pi_no_Prot_mom_cut";

        TH1D * void1 = e_frac(t[17],c,attr1,Form("Gib1%s",attr1[1].c_str()));
        //TH1D * void2 = e_tail(t[17],c,attr1,Form("Gib2%s",attr1[1].c_str()));
        TH1D * void3 = e_kin(t[17],c,attr1,Form("Gib3%s",attr1[1].c_str()));
        sanity_check(t[17],c,attr1,Form("Gib3%s",attr1[1].c_str()));
        q2_bjorkenx_comp(t[17],c,attr1,Form("Gib4%s",attr1[1].c_str()));
        q0_q3_comp(t[17],c,attr1,Form("Gib5%s",attr1[1].c_str()));

        if(attr1[1] == "GiBUU_CC0Pi_cut"){ //outputs the fraction of events in certain energy bands
          cout << Form("\n|            | tgt:Ar | cut: %s |  Energy: 2.261 |\n",attr1[1].c_str());
          cout << "|   Energy   | 0.75-1 | 1-1.25 | 1.25-1.5 | 1.5-1.75 | 1.75-2 |\n";
          cout << Form("| Event Frac | %.4f | %.4f |  %.4f  |  %.4f  | %.4f |\n",void1->Integral(27,36),void1->Integral(37,45),void1->Integral(46,53),void1->Integral(54,62),void1->Integral(63,72));
        }

        vector<string> FSParticles;
        FSParticles.push_back("EFSProt");
        FSParticles.push_back("EFSLep");
        FSParticles.push_back("EFSNeutron");
        FSParticles.push_back("EFSPion");
        FSParticles.push_back("EFSOther");

        for(int k = 0; k < FSParticles.size(); ++k){
          e_cont(t[17],c,attr1,Form("Gib%s%s",FSParticles[k].c_str(),attr1[1].c_str()),FSParticles[k].c_str());
        }
      }
    }
  }

}

/////////////////////////////////////////////////////////////////////////////////////////////////

void get_tables(){
  /*
  outputs the fraction of events in certain energy bands for the files in t
  */

  vector<TTree *>t;
  t.push_back(get_tree("../CLAS6gens/2.2_Ar_GENIE_escat_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/2.2_Fe_GENIE_escat_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/GiBUU.escat.Ar40.NoFSI.stdhep.CLAS6.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/2.2_C12_GENIE_escat_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/4.4_Fe_GENIE_escat_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/4.4_C12_GENIE_escat_CLAS6Acc.root","CLAS6AcceptSummary"));

  char * sel = "EISLepRec_calo";
  char * cut = "((NEUTMode == 1) || (NEUTMode == 2) || (NEUTMode == 11) || (NEUTMode == 12) || (NEUTMode == 13)) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96) && (FSProton_mom > 300) && (W < 2000) && (Pperp > 200)";

  char * attr[6][2] = 
    {{"Ar","GENIE"},
    {"Fe","GENIE"},
    {"Ar","GiBUU"},
    {"C12","GENIE"},
    {"Fe","GENIE"},
    {"C12","GENIE"}};

  char * bounds2[6] = {"(EISLepRec_calo >= 750) && (EISLepRec_calo <= 1000)","(EISLepRec_calo >= 1000) && (EISLepRec_calo <= 1250)","(EISLepRec_calo >= 1250) && (EISLepRec_calo <= 1500)","(EISLepRec_calo >= 1500) && (EISLepRec_calo <= 1750)","(EISLepRec_calo >= 1750) && (EISLepRec_calo <= 2000)","(EISLepRec_calo >= 0) && (EISLepRec_calo <= 3000)"};
  char * bounds4[6] = {"(EISLepRec_calo >= 1500) && (EISLepRec_calo <= 2000)","(EISLepRec_calo >= 2000) && (EISLepRec_calo <= 2500)","(EISLepRec_calo >= 2500) && (EISLepRec_calo <= 3000)","(EISLepRec_calo >= 3000) && (EISLepRec_calo <= 3500)","(EISLepRec_calo >= 3500) && (EISLepRec_calo <= 4000)","(EISLepRec_calo >= 0) && (EISLepRec_calo <= 5000)"};

  cout << setprecision(3) << setw(8);
  for(int i = 0; i < 6; ++i){
    if(i < 4){
      cout << Form("|  inter:nu  |  tgt:%s  |     Energy: 2.2     |   Generator:%s   |\n",attr[i][0],attr[i][1]);
      cout << "|   Energy   | 0.75-1.0 | 1.0-1.25 | 1.25-1.5 | 1.5-1.75 | 1.75-2.0 | N_events |\n";
      cout << "|            | ";
      for(int j = 0; j < 6; ++j){
        cout << get_g(t[i],sel,cut,bounds2[j],Form("2%d%s%s",j,attr[i][0],attr[i][1]))->Integral() << " | ";
      }
      cout << "\n\n";
    }

    else{
      cout << Form("|  inter:nu  |  tgt:%s  |     Energy: 4.4     |   Generator:%s   |\n",attr[i][0],attr[i][1]);
      cout << "|   Energy   |   1.5-2  |   2-2.5  |   2.5-3  |   3-3.5  |   3.5-4  | N_events |\n";
      cout << "|            | ";
      for(int j = 0; j < 6; ++j){
        cout << get_g(t[i],sel,cut,bounds4[j],Form("%d%s%s",j,attr[i][0],attr[i][1]))->Integral() << " | ";
      }
      cout << "\n\n";
    }
  }

  
}

////////////////////////////////////////////////////////////////////////////////////////////////

TH1 * get_g(TTree * t,char * sel,char * cut,char * bounds,char * temp){
  /*
  Gets the graph given certain attributes
  */
  t->Draw(Form("%s >> %s",sel,temp),Form("%s && %s",cut,bounds));
  TH1 * g = (TH1 *)gDirectory->Get(temp);
  return g;
}

////////////////////////////////////////////////////////////////////////////////////////////////

void make_C_graphs(){
  /*
  Makes graphs for just the carbon events
  */
  TCanvas * c = new TCanvas();
  vector<TTree *> t;
  t.push_back(get_tree("../CLAS6gens/2.2_C12_GENIE_escat_CLAS6Acc.root","CLAS6AcceptSummary"));
  t.push_back(get_tree("../CLAS6gens/4.4_C12_GENIE_escat_CLAS6Acc.root","CLAS6AcceptSummary"));

  // Calorimetric cut with proton momentum cut
  char * calocut = "((NEUTMode == 1) || (NEUTMode == 2) || (NEUTMode == 11) || (NEUTMode == 12) || (NEUTMode == 13)) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96) && (FSProton_mom > 300) && (W < 2000) && (Pperp > 200)";

  // Calorimetric cut without proton momentum cut
  char * kincut = "((NEUTMode == 1) || (NEUTMode == 2) || (NEUTMode == 11) || (NEUTMode == 12) || (NEUTMode == 13)) && (Q2_true > 5E5) && (fabs(bjorken_x_rec - 1) < 0.2) && (FSLepton_ct < 0.96) && (W < 2000) && (Pperp > 200)";

  vector<string> attr;
  attr.push_back(calocut);
  attr.push_back("CC0Pi_eep_cut");
  attr.push_back("2.261");
  attr.push_back("C12");
  attr.push_back("foo1");
  attr.push_back("scat");

  TH1D * void1 = e_frac(t[0],c,attr,"tmp1");
  attr[0] = kincut;
  TH1D * void2 = e_kin(t[0],c,attr,"tmp2");
  attr[2] = "4.4";
  TH1D * void3 = e_kin(t[1],c,attr,"tmp3");
  attr[0] = calocut;
  TH1D * void4 = e_frac(t[1],c,attr,"tmp4");
}

////////////////////////////////////////////////////////////////////////////////////////////////

