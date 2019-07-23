#include<../inclusions.h>
#include<vector>
#include<iostream>
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////

void e_comp(TTree * t1, TTree * t2, vector<string> attr1, vector<string> attr2, TCanvas * c, TLegend * l1, TLegend * l2);
TH1D * e_frac(TTree * t, TCanvas * c, vector<string> attr, string tmp);
TH1D * e_kin(TTree * t, TCanvas * c, vector<string> attr, string tmp);
TH1D * e_tail(TTree * t, TCanvas * c, vector<string> attr, string tmp);
TH1D * e_lost(TTree * t, TCanvas * c, vector<string> attr, string tmp);
void e_cont(TTree * t, TCanvas * c, vector<string> attr, string tmp, string cont);
void e_cont_diff(TTree * t, TCanvas * c, vector<string> attr, string tmp, string cont);
TH1D * sanity_check(TTree * t, TCanvas * c, vector<string> attr, string tmp);
void comp_energy_ranges(TTree * f0, TTree * f1, TTree * f2, TTree * f3, TTree * a0, TTree * a1, TTree * a2, TTree * a3, TCanvas * c, TLegend * l, vector<string> attrf, vector<string> attra);
void e_lep_prot(TTree * t, TCanvas * c, TLegend *l, vector<string> attr);
void peak_diff(TCanvas *c, vector<string> attr, TLegend *l, TH1D * e_total, TH1D * e_reco);
void epion(TTree * t, TCanvas * c, vector<string> attr, string tmp);

/////////////////////////////////////////////////////////////////////////////////////////////////

void e_comp(TTree * t1, TTree * t2, vector<string> attr1, vector<string> attr2, TCanvas * c, TLegend * l1, TLegend * l2){
  /* Takes two file name in the current directory and will create 6 png's of the ereco comparison
   * to e_true using a few selections specified and outputs them in the current directory
   */

  string selection = attr1[0];
  string sel_title = attr1[1];
  string Enu = attr1[2];
  string inter = attr1[5];
  string Enu_t;

  if(Enu == "2.261") Enu_t = "2_2";
  else if(Enu == "3.3") Enu_t = "3_3";
  else if(Enu == "4.4") Enu_t = "4_4";
  else if(Enu == "1.1") Enu_t = "1_1";

  epion(t2,c,attr2,Form("epion%s1%s",attr2[4].c_str(),Enu.c_str()));
  epion(t1,c,attr1,Form("epion%s2%s",attr1[4].c_str(),Enu.c_str()));

  TH1D * Fe = e_frac(t2,c,attr2,Form("efrac%s1%s",attr2[4].c_str(),Enu.c_str()));
  TH1D * Ar = e_frac(t1,c,attr1,Form("efrac%s2%s",attr1[4].c_str(),Enu.c_str()));

  TH1D * Fe_kin = e_kin(t2,c,attr2,Form("ekin%s1%s",attr2[4].c_str(),Enu.c_str()));
  TH1D * Ar_kin = e_kin(t1,c,attr1,Form("ekin%s2%s",attr1[4].c_str(),Enu.c_str()));

  l1->AddEntry(Fe,"Iron");
  l1->AddEntry(Ar,"Argon");

  Ar->SetLineColor(2);
  Fe->SetLineColor(4);
  Ar->SetLineWidth(2);
  Fe->SetLineWidth(2);
  Ar->Draw();
  Fe->Draw("same");
  l1->Draw("same");
  c->SaveAs(Form("../GENIEescatcomps/pngdump/efrac/%s_%s_Fe_Ar_ereco_%s_comp.png",inter.c_str(),Enu_t.c_str(),sel_title.c_str()));
  l1->Clear();

  TH1D * Fe_tail = e_tail(t2,c,attr2,Form("tail%s1%s",attr2[4].c_str(),Enu.c_str()));
  TH1D * Ar_tail = e_tail(t1,c,attr1,Form("tail%s2%s",attr1[4].c_str(),Enu.c_str()));

  l2->AddEntry(Fe_tail,"Iron");
  l2->AddEntry(Ar_tail,"Argon");

  Ar_tail->SetLineColor(2);
  Fe_tail->SetLineColor(4);
  Ar_tail->SetLineWidth(2);
  Fe_tail->SetLineWidth(2);

  Fe_tail->Draw();
  Ar_tail->Draw("same");
  l2->Draw("same");
  c->SaveAs(Form("../GENIEescatcomps/pngdump/etail/%s_%s_Fe_Ar_etail_%s_comp.png",inter.c_str(),Enu_t.c_str(),sel_title.c_str()));
  l2->Clear();

  TH1D * Fe_tot = sanity_check(t2,c,attr2,Form("san%s1%s",attr2[4].c_str(),Enu.c_str()));
  TH1D * Ar_tot = sanity_check(t1,c,attr1,Form("san%s2%s",attr1[4].c_str(),Enu.c_str()));

  e_lep_prot(t1,c,l1,attr1);
  e_lep_prot(t2,c,l1,attr2);

  peak_diff(c, attr2, l1, Fe_tot, Fe);
  peak_diff(c, attr1, l1, Ar_tot, Ar);

  vector<string> FSParticles;
  FSParticles.push_back("EFSProt");
  FSParticles.push_back("EFSLep");
  FSParticles.push_back("EFSNeutron");
  FSParticles.push_back("EFSPion");
  FSParticles.push_back("EFSOther");

  for(int i = 0; i < FSParticles.size(); ++i){
    e_cont(t2,c,attr2,Form("%s4%s%s%s",attr2[4].c_str(),FSParticles[i].c_str(),Enu.c_str(),FSParticles[i]),FSParticles[i]);
    e_cont(t1,c,attr1,Form("%s5%s%s%s",attr1[4].c_str(),FSParticles[i].c_str(),Enu.c_str(),FSParticles[i]),FSParticles[i]);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////

TH1D * e_frac(TTree * t, TCanvas * c, vector<string> attr, string tmp){
  /* Takes a file name in the current directory and if given a selection, will create and return
  * a (E_reco - Etrue)/Etrue TH1D as a pointer with E_reco-Etrue/Etrue on the x-axis and the
  * # of events on the y-axis
  */

  string selection = attr[0];
  string sel_title = attr[1];
  string Enu = attr[2];
  string tgt = attr[3];
  string inter = attr[5];
  string Enu_t;
  double Enu_d = 5000;

  if(Enu == "2.261"){Enu_t = "2_2"; Enu_d = 2200;}
  else if(Enu == "3.3"){Enu_t = "3_3"; Enu_d = 3300;}
  else if(Enu == "4.4"){Enu_t = "4_4"; Enu_d = 4400;}
  else if(Enu == "1.1"){Enu_t = "1_1"; Enu_d = 1100;}
  TH1D * comp = new TH1D(tmp.c_str(),"",100,-0.4,0.2);
  t->Draw(Form("(EISLepRec_calo - EISLepTrue) / EISLepTrue >> %s",tmp.c_str()), selection.c_str()); //flag
  plot_title_style(comp,Form("Energy Reconstructed (Ereco - Etrue)/Etrue (E_{#nu} %s GeV)", Enu.c_str()),"Energy Reconstructed","# of Events");
  Double_t scale = 1/comp->Integral();
  comp->Scale(scale);
  comp->Draw();

  TPaveText *pt = new TPaveText(0.15,0.75,0.6,0.85,"brNDC");
//  pt->SetFillColor(18);
  pt->SetTextAlign(12);
  pt->AddText(Form("Fraction of events of .05 of True Enu: %.4f",comp->Integral(60,75)));
  pt->Draw("same");
  c->SaveAs(Form("../GENIEescatcomps/pngdump/efrac/%s_%s_%s_ereco_%s.png",inter.c_str(),Enu_t.c_str(),tgt.c_str(),sel_title.c_str()));
  c->Clear();
  return (TH1D*)comp->Clone();
}

/////////////////////////////////////////////////////////////////////////////////////////////////

TH1D * e_kin(TTree * t, TCanvas * c, vector<string> attr, string tmp){
  /* Takes a file name in the current directory and if given a selection, will create and return
  * a (E_reco - Etrue)/Etrue TH1D as a pointer with E_reco-Etrue/Etrue on the x-axis and the
  * # of events on the y-axis
  */

  string selection = attr[0];
  string sel_title = attr[1];
  string Enu = attr[2];
  string tgt = attr[3];
  string inter = attr[5];
  string Enu_t;
  double Enu_d = 5000;

  if(Enu == "2.261"){Enu_t = "2_2"; Enu_d = 2200;}
  else if(Enu == "3.3"){Enu_t = "3_3"; Enu_d = 3300;}
  else if(Enu == "4.4"){Enu_t = "4_4"; Enu_d = 4400;}
  else if(Enu == "1.1"){Enu_t = "1_1"; Enu_d = 1100;}
  TH1D * comp = new TH1D(tmp.c_str(),"",100,-1,0.6);
  t->Draw(Form("(EISLepRec_KineFSLep - EISLepTrue) / EISLepTrue >> %s",tmp.c_str()), selection.c_str()); //flag
  plot_title_style(comp,Form("Energy Reconstructed (Ekin - Etrue)/Etrue (E_{%s} %s GeV)",inter.c_str(),Enu.c_str()),"Energy Reconstructed","# of Events");
  Double_t scale = 1/comp->Integral();
  comp->Scale(scale);
  comp->Draw();

  TPaveText *pt = new TPaveText(0.15,0.75,0.6,0.85,"brNDC");
//  pt->SetFillColor(18);
  pt->SetTextAlign(12);
  pt->AddText(Form("Fraction of events of .05 of True Enu: %.4f",comp->Integral(60,67)/comp->Integral()));
  pt->Draw("same");
  c->SaveAs(Form("../GENIEescatcomps/pngdump/ekin/%s_%s_%s_ekin_%s.png",inter.c_str(),Enu_t.c_str(),tgt.c_str(),sel_title.c_str()));
  c->Clear();
  return (TH1D*)comp->Clone();
}

/////////////////////////////////////////////////////////////////////////////////////////////////

TH1D * e_tail(TTree * t, TCanvas * c, vector<string> attr, string tmp){
  /* Takes a file name in the current directory and if given a selection, will create and return
   * a E_reco - Etrue TH1D as a pointer
   */

  string selection = attr[0];
  string sel_title = attr[1];
  string Enu = attr[2];
  string tgt = attr[3];
  string inter = attr[5];
  string Enu_t;
  double Enu_d = 5000;

  if(Enu == "2.261"){Enu_t = "2_2"; Enu_d = 2200;}
  else if(Enu == "3.3"){Enu_t = "3_3"; Enu_d = 3300;}
  else if(Enu == "4.4"){Enu_t = "4_4"; Enu_d = 4400;}
  else if(Enu == "1.1"){Enu_t = "1_1"; Enu_d = 1100;}

  TH1D * comp = new TH1D(tmp.c_str(),"",100,0.0001,100.0001);
  t->Draw(Form("EISLepRec_calo - EISLepTrue >> %s",tmp.c_str()), selection.c_str());
  plot_title_style(comp,Form("Energy Reconstructed (Ereco - Etrue) (E_{%s} %s GeV)",inter.c_str(),Enu.c_str()),"Energy Reconstructed (MeV)","# of Events");
  Double_t scale = 1/1E6;
  comp->Scale(scale);
  comp->Draw();
  c->SaveAs(Form("../GENIEescatcomps/pngdump/etail/%s_%s_%s_etail_%s.png",inter.c_str(),Enu_t.c_str(),tgt.c_str(),sel_title.c_str()));
  c->Clear();
  return (TH1D*)comp->Clone();
}

/////////////////////////////////////////////////////////////////////////////////////////////////

TH1D * e_lost(TTree * t, TCanvas * c, vector<string> attr, string tmp){
  /* Takes a file name in the current directory and if given a selection, will create and return
   * a E_reco - Etrue TH1D as a pointer
   */

  string selection = attr[0];
  string sel_title = attr[1];
  string Enu = attr[2];
  string tgt = attr[3];
  string inter = attr[5];
  string Enu_t;
  double Enu_d = 5000;

  if(Enu == "2.261"){Enu_t = "2_2"; Enu_d = 2200;}
  else if(Enu == "3.3"){Enu_t = "3_3"; Enu_d = 3300;}
  else if(Enu == "4.4"){Enu_t = "4_4"; Enu_d = 4400;}
  else if(Enu == "1.1"){Enu_t = "1_1"; Enu_d = 1100;}

  TH1D * comp = new TH1D(tmp.c_str(),"",100,-(Enu_d + 100),400);
  t->Draw(Form("EISLepRec_calo - EISLepTrue >> %s",tmp.c_str()), selection.c_str());
  plot_title_style(comp,Form("Energy Reconstructed (Ereco - Etrue) (E_{%s} %s GeV)",inter.c_str(),Enu.c_str()),"Energy Reconstructed (MeV)","# of Events");
  Double_t scale = 1/1E6;
  comp->Scale(scale);
  comp->Draw();
  c->SaveAs(Form("../GENIEescatcomps/pngdump/elost/%s_%s_%s_elost_%s.png",inter.c_str(),Enu_t.c_str(),tgt.c_str(),sel_title.c_str()));
  c->Clear();
  return (TH1D*)comp->Clone();
}

/////////////////////////////////////////////////////////////////////////////////////////////////

TH1D * sanity_check(TTree * t, TCanvas * c, vector<string> attr, string tmp){
  /* Takes a file name in the current directory and if given a selection, will create and return
  * a (E_reco - Etrue)/Etrue TH1D as a pointer with E_reco-Etrue/Etrue on the x-axis and the
  * # of events on the y-axis
  */

  string selection = attr[0];
  string sel_title = attr[1];
  string Enu = attr[2];
  string tgt = attr[3];
  string inter = attr[5];
  string Enu_t;
  double Enu_d = 5000;
  if(Enu == "2.261"){Enu_t = "2_2"; Enu_d = 2261;}
  else if(Enu == "3.3"){Enu_t = "3_3"; Enu_d = 3300;}
  else if(Enu == "4.4"){Enu_t = "4_4"; Enu_d = 4400;}
  else if(Enu == "1.1"){Enu_t = "1_1"; Enu_d = 1100;}

  TH1D * comp = new TH1D(tmp.c_str(),"",100,-0.4,0.2);
  t->Draw(Form("(EISLepRec_calo + EFSNeutron + EFSPion - EISLepTrue) / EISLepTrue >> %s",tmp.c_str()), selection.c_str()); //flag
  comp->Scale(1/comp->Integral());
  comp->Draw("same");
  plot_title_style(comp,Form("Total Final Energy (Ereco + Neutron - Etrue)/Etrue (E_{%s} %s GeV)",inter.c_str(),Enu.c_str()),"Energy Reconstructed","# of Events");

  TFitResultPtr r = comp->Fit("gaus","SQN");
  if(r.Get() != NULL){
    double width = 2 * r->Parameter(2);
    TPaveText *pt = new TPaveText(0.2,0.7,0.5,0.9,"brNDC");
//    pt->SetFillColor(18);
    pt->SetTextAlign(12);
    pt->AddText(Form("Width of Gaussian fit: %f",width*Enu_d));
    pt->AddText(Form("At a mean of: %.4f",r->Parameter(1)));
    pt->Draw("same");
  }
  c->SaveAs(Form("../GENIEescatcomps/pngdump/etotal/%s_%s_%s_etotal_%s.png",inter.c_str(),Enu_t.c_str(),tgt.c_str(),sel_title.c_str()));
  c->Clear();
  return comp;
}

/////////////////////////////////////////////////////////////////////////////////////////////////

void comp_energy_ranges(TTree * f0, TTree * f1, TTree * f2, TTree * f3, TTree * a0, TTree * a1, TTree * a2, TTree * a3, TCanvas * c, TLegend * l, vector<string> attrf, vector<string> attra){
  /* Takes multiple files containing energies 1.1, 2.2, 3.3, and 4.4 GeV in that order and makes
   * a histogram using the e_frac method in order to compare their shapes by normalizing both the
   * x and y-axis. Then it saves the historgram as a png in the pngdump folder
   */

  string selection = attrf[0];
  string sel_title = attrf[1];
  string Enu = attrf[2];
  string tgtf = attrf[3];
  string tgta = attra[3];
  string inter = attrf[5];
  if((f0 == 0x0) || (f1 == 0x0) || (f2 == 0x0) || (f3 == 0x0) || (a0 == 0x0) || (a1 == 0x0) || (a2 == 0x0) || (a3 == 0x0))continue;
  TH1D * Ef1 = e_frac(f0,c,attrf,Form("comp0%s%s%s",tgtf.c_str(),sel_title.c_str(),Enu.c_str()));
  TH1D * Ef2 = e_frac(f1,c,attrf,Form("comp1%s%s%s",tgtf.c_str(),sel_title.c_str(),Enu.c_str()));
  TH1D * Ef3 = e_frac(f2,c,attrf,Form("comp2%s%s%s",tgtf.c_str(),sel_title.c_str(),Enu.c_str()));
  TH1D * Ef4 = e_frac(f3,c,attrf,Form("comp3%s%s%s",tgtf.c_str(),sel_title.c_str(),Enu.c_str()));
  TH1D * Ea1 = e_frac(a0,c,attra,Form("comp0%s%s%s",tgta.c_str(),sel_title.c_str(),Enu.c_str()));
  TH1D * Ea2 = e_frac(a1,c,attra,Form("comp1%s%s%s",tgta.c_str(),sel_title.c_str(),Enu.c_str()));
  TH1D * Ea3 = e_frac(a2,c,attra,Form("comp2%s%s%s",tgta.c_str(),sel_title.c_str(),Enu.c_str()));
  TH1D * Ea4 = e_frac(a3,c,attra,Form("comp3%s%s%s",tgta.c_str(),sel_title.c_str(),Enu.c_str()));


  Ef1->SetLineWidth(2);
  Ef2->SetLineWidth(2);
  Ef3->SetLineWidth(2);
  Ef4->SetLineWidth(2);
  
  Ef1->SetLineColor(1);
  Ef2->SetLineColor(2);
  Ef3->SetLineColor(4);
  Ef4->SetLineColor(8);

  Ef1->SetLineStyle(1);
  Ef2->SetLineStyle(10);
  Ef3->SetLineStyle(9);
  Ef4->SetLineStyle(6);

  l->AddEntry(Ef1,Form("%s (E_{nu} = 1.1GeV)",tgtf.c_str()));
  l->AddEntry(Ef2,Form("%s (E_{nu} = 2.261GeV)",tgtf.c_str()));
  l->AddEntry(Ef3,Form("%s (E_{nu} = 3.3GeV)",tgtf.c_str()));
  l->AddEntry(Ef4,Form("%s (E_{nu} = 4.4GeV)",tgtf.c_str()));

  Ef1->Draw();
  Ef2->Draw();
  Ef3->Draw("same");
  Ef4->Draw("same");
  l->Draw("same");
  c->SaveAs(Form("../GENIEescatcomps/pngdump/efrac/%s_%s_Enu_%s_comp.png",inter.c_str(),tgtf.c_str(),sel_title.c_str()));
  l->Clear();
  c->Clear();
  
  Ea1->SetLineWidth(2);
  Ea2->SetLineWidth(2);
  Ea3->SetLineWidth(2);
  Ea4->SetLineWidth(2);
  
  Ea1->SetLineColor(1);
  Ea2->SetLineColor(2);
  Ea3->SetLineColor(4);
  Ea4->SetLineColor(8);

  Ea1->SetLineStyle(1);
  Ea2->SetLineStyle(10);
  Ea3->SetLineStyle(9);
  Ea4->SetLineStyle(6);

  l->AddEntry(Ea1,Form("%s (E_{nu} = 1.1GeV)",tgta.c_str()));
  l->AddEntry(Ea2,Form("%s (E_{nu} = 2.261GeV)",tgta.c_str()));
  l->AddEntry(Ea3,Form("%s (E_{nu} = 3.3GeV)",tgta.c_str()));
  l->AddEntry(Ea4,Form("%s (E_{nu} = 4.4GeV)",tgta.c_str()));

  Ea1->Draw();
  Ea2->Draw();
  Ea3->Draw("same");
  Ea4->Draw("same");
  l->Draw("same");
  c->SaveAs(Form("../GENIEescatcomps/pngdump/efrac/%s_%s_Enu_%s_comp.png",inter.c_str(),tgtf.c_str(),sel_title.c_str()));
  l->Clear();
  c->Clear();

  Ef1->Add(Ea1,-1);
  Ef2->Add(Ea2,-1);
  Ef3->Add(Ea3,-1);
  Ef4->Add(Ea4,-1);
  Ef1->SetLineWidth(2);
  Ef2->SetLineWidth(2);
  Ef3->SetLineWidth(2);
  Ef4->SetLineWidth(2);
  
  Ef1->SetLineColor(1);
  Ef2->SetLineColor(2);
  Ef3->SetLineColor(4);
  Ef4->SetLineColor(8);

  Ef1->SetLineStyle(1);
  Ef2->SetLineStyle(10);
  Ef3->SetLineStyle(9);
  Ef4->SetLineStyle(6);

  l->AddEntry(Ef1,"Fe-Ar (E_{nu} = 1.1GeV)");
  l->AddEntry(Ef2,"Fe-Ar (E_{nu} = 2.261GeV)");
  l->AddEntry(Ef3,"Fe-Ar (E_{nu} = 3.3GeV)");
  l->AddEntry(Ef4,"Fe-Ar (E_{nu} = 4.4GeV)");

  Ef1->GetYaxis()->SetLimits(-0.1,0.1);
  Ef2->GetYaxis()->SetLimits(-0.1,0.1);
  Ef3->GetYaxis()->SetLimits(-0.1,0.1);
  Ef4->GetYaxis()->SetLimits(-0.1,0.1);

  Ef1->Draw();
  Ef2->Draw("same");
  Ef3->Draw("same");
  Ef4->Draw("same");
  l->Draw("same");
  c->SaveAs(Form("../GENIEescatcomps/pngdump/efrac/%s_Ar_Fe_diff_Enu_%s_comp.png",inter.c_str(),sel_title.c_str()));
  l->Clear();
  c->Clear();
 
}

/////////////////////////////////////////////////////////////////////////////////////////////////

void e_cont(TTree * t, TCanvas * c, vector<string> attr, string tmp, string cont){
 /* Takes a file name in the current directory and if given a selection, will create and return
   * a E_reco - Etrue TH1D as a pointer
   */

  string selection = attr[0];
  string sel_title = attr[1];
  string Enu = attr[2];
  string tgt = attr[3];
  string inter = attr[5];
  string Enu_t;
  double Enu_d = 1000;
  
  if(Enu == "2.261"){Enu_t = "2_2"; Enu_d = 2200;}
  else if(Enu == "3.3"){Enu_t = "3_3"; Enu_d = 3300;}
  else if(Enu == "4.4"){Enu_t = "4_4"; Enu_d = 4400;}
  else if(Enu == "1.1"){Enu_t = "1_1"; Enu_d = 1100;}

  TH2D * comp = new TH2D(tmp.c_str(),"",100,0.0001,1.2001,100,0.0001,1.0001);
  t->Draw(Form("%s/EISLepTrue:EISLepRec_calo/EISLepTrue >> %s",cont.c_str(),tmp.c_str()), selection.c_str());
  plot_title_style(comp,Form("Energy Reco Contribution of %s (E_{%s} %s GeV)",inter.c_str(),cont.c_str(),Enu.c_str()),"Energy Reconstruction (MeV)",Form("%s",cont.c_str()));
  comp->Scale(1/1E6);
  comp->Draw("colz");
  c->SaveAs(Form("../GENIEescatcomps/pngdump/econt/%s/%s_%s_%s_econt_%s_%s.png",cont.c_str(),inter.c_str(),Enu_t.c_str(),tgt.c_str(),cont.c_str(),sel_title.c_str()));
  c->Clear();

}

/////////////////////////////////////////////////////////////////////////////////////////////////

void e_lep_prot(TTree * t, TCanvas * c, TLegend * l, vector<string> attr){
  /* Takes a file name in the current directory and if given a selection, will create and return
   * a E_reco - Etrue TH1D as a pointer
   */

  int temp_col[] = {1,2,nOrange,nSkyBlue,nBlueGreen,nYellow,nBlue,nVermillion,nRedPurple,4,9,11,12,13,41,40,33,kRed-4,kYellow+4,kYellow-1,kCyan,kCyan+4,kGray};
  string selection = attr[0];
  string sel_title = attr[1];
  string Enu = attr[2];
  string tgt = attr[3];
  string inter = attr[5];
  string Enu_t;
  int Enu_d = 1000;
  
  if(Enu == "2.261"){Enu_t = "2_2"; Enu_d = 2261;}
  else if(Enu == "3.3"){Enu_t = "3_3"; Enu_d = 3300;}
  else if(Enu == "4.4"){Enu_t = "4_4"; Enu_d = 4400;}
  else if(Enu == "1.1"){Enu_t = "1_1"; Enu_d = 1100;}

  int e_low = 0;
  const int num = Enu_d/200;
  TH1D * comp[num];
  l->SetNColumns(2);
  for(int iter = 0; iter < num; ++iter){
    string tmp = Form("%s%s%s%d%s",inter.c_str(),tgt.c_str(),Enu.c_str(),e_low,sel_title.c_str());
    e_low += 200;
    comp[iter] = new TH1D(tmp.c_str(),"",100,0,1.2);
    t->Draw(Form("EFSLep/EISLepTrue >> %s",tmp.c_str()), Form("%s && EFSProt > %d && EFSProt < %d",selection.c_str(),e_low,e_low+200));
    comp[iter]->Scale(1/1E6);
  }
  comp[0]->SetLineWidth(2);
  comp[0]->SetLineColor(1);
  comp[0]->Draw();
  plot_title_style(comp[0],Form("Energy Reco Contribution of EFSLep (E_{%s} %s GeV)",inter.c_str(),Enu.c_str()),"Energy Reconstruction (MeV)","Fraction of Final State Energy");
  e_low = 0;
  l->AddEntry(comp[0],Form("EFSProt > %d",e_low));
  for(int i = 1; i < num; ++i){
    e_low+= 200;
    comp[i]->SetLineWidth(2);
    comp[i]->SetLineColor(temp_col[i]);
    comp[i]->Draw("same");
    l->AddEntry(comp[i],Form("EFSProt > %d",e_low));
  }
  l->Draw("same");
  c->SaveAs(Form("../GENIEescatcomps/pngdump/eleprot/%s_%s_%s_eleprot_EFSLep_%s.png",inter.c_str(),Enu_t.c_str(),tgt.c_str(),sel_title.c_str()));
  c->Clear();
  l->Clear();
  l->SetNColumns(1);
}

/////////////////////////////////////////////////////////////////////////////////////////////////

void peak_diff(TCanvas *c, vector<string> attr, TLegend *l, TH1D * e_total, TH1D * e_reco){

  string selection = attr[0];
  string sel_title = attr[1];
  string Enu = attr[2];
  string tgt = attr[3];
  string inter = attr[5];
  string Enu_t;
  int Enu_d = 1000;
  
  if(Enu == "2.261"){Enu_t = "2_2"; Enu_d = 2261;}
  else if(Enu == "3.3"){Enu_t = "3_3"; Enu_d = 3300;}
  else if(Enu == "4.4"){Enu_t = "4_4"; Enu_d = 4400;}
  else if(Enu == "1.1"){Enu_t = "1_1"; Enu_d = 1100;}

  TPaveText *pt = new TPaveText(0.2,0.3,0.6,0.5,"brNDC");
//  pt->SetFillColor(18);
  pt->SetTextAlign(12);

  TFitResultPtr r1 = e_total->Fit("gaus","SQN");
  TFitResultPtr r2 = e_reco->Fit("gaus","SQN");
  if(r1.Get() != NULL){
    if(r2.Get() != NULL){
      pt->AddText(Form("Difference between peaks: %.4f",r1->Parameter(1) - r2->Parameter(1)));
      pt->AddText(Form("Or in MeV: %.4f MeV",(r1->Parameter(1) - r2->Parameter(1)) * Enu_d));
      pt->Draw("same");
    }
  }

  e_total->SetLineColor(2);
  e_reco->SetLineColor(4);
  e_total->SetLineWidth(2);
  e_reco->SetLineWidth(2);

  plot_title_style(e_total,Form("Energy Reconstruction (E_{%s} %s GeV)",inter.c_str(),Enu.c_str()),"Energy Reconstruction","# Of Events");

  l->AddEntry(e_total,"Total Energy");
  l->AddEntry(e_reco,"Energy Reconstructed");

  e_reco->Draw();
  e_total->Draw("same");
  pt->Draw("same");
  l->Draw("same");
  c->SaveAs(Form("../GENIEescatcomps/pngdump/peakdiff/%s_%s_%s_peak_diff_%s.png",inter.c_str(),Enu_t.c_str(),tgt.c_str(),sel_title.c_str()));
  c->Clear();
  l->Clear();
}

/////////////////////////////////////////////////////////////////////////////////////////////////

void epion(TTree * t, TCanvas * c, vector<string> attr, string tmp){

  string selection = attr[0];
  string sel_title = attr[1];
  string Enu = attr[2];
  string tgt = attr[3];
  string inter = attr[5];
  string Enu_t;
  double Enu_d = 5000;

  if(Enu == "2.261"){Enu_t = "2_2"; Enu_d = 2200;}
  else if(Enu == "3.3"){Enu_t = "3_3"; Enu_d = 3300;}
  else if(Enu == "4.4"){Enu_t = "4_4"; Enu_d = 4400;}
  else if(Enu == "1.1"){Enu_t = "1_1"; Enu_d = 1100;}
  TH1D * comp = new TH1D(tmp.c_str(),"",100,-0.9999,0.2501);
  t->Draw(Form("(EISLepRec_calo) >> %s",tmp.c_str()), Form("(EFSPion < 1400) && %s",selection.c_str())); //flag
  plot_title_style(comp,Form("Energy Reconstructed with E_{lep} < 1400MeV(E_{#nu} %s GeV)", Enu.c_str()),"Energy Reconstructed (MeV)","# of Events");
  Double_t scale = 1/1E6;
  comp->Scale(scale);
  comp->Draw();
  c->SaveAs(Form("../GENIEescatcomps/pngdump/efrac/%s_%s_%s_elep_%s.png",inter.c_str(),Enu_t.c_str(),tgt.c_str(),sel_title.c_str()));
  c->Clear();
}

/////////////////////////////////////////////////////////////////////////////////////////////////

#endif
