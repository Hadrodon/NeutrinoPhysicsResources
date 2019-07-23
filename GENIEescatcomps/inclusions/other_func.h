#ifndef OTHER_FUNC_H
#define OTHER_FUNC_H
#include<../inclusions.h>
#include<vector>
#include<iostream>
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////

void other_comp(TTree * t1, TTree * t2, vector<string> attr1, vector<string> attr2, TCanvas * c, vector<string> extra1, vector<string> extra2);
void q2_bjorkenx_comp(TTree * t, TCanvas * c, vector<string> attr, string tmp);
void q0_q3_comp(TTree * t, TCanvas * c, vector<string> attr, string tmp);
void q2_q0_comp(TTree * t, TCanvas * c, vector<string> attr, string tmp);
void view_constraint_prog(TTree * t, TCanvas * c, TLegend * l, vector<string> attr, vector<string> sel, vector<string> sel_ttl);
void view_pperp(TTree * t, TCanvas * c, TLegend * l, vector<string> attr, vector<string> sel, vector<string> sel_ttl);

/////////////////////////////////////////////////////////////////////////////////////////////////

void other_comp(TTree * t1, TTree * t2, vector<string> attr1, vector<string> attr2, TCanvas * c, vector<string> extra1, vector<string> extra2){
  TLegend * l = new TLegend(0.1,0.7,0.5,0.9);
  q2_bjorkenx_comp(t2,c,attr2,Form("1%s%s",attr2[4].c_str(),attr2[2].c_str()));
  q2_bjorkenx_comp(t1,c,attr1,Form("2%s%s",attr1[4].c_str(),attr1[2].c_str()));
  q0_q3_comp(t2,c,attr2,Form("temp%s%s",attr2[4].c_str(),attr2[2].c_str()));
  q0_q3_comp(t1,c,attr1,Form("temp%s%s",attr1[4].c_str(),attr1[2].c_str()));
  view_constraint_prog(t1,c,l,attr1,extra1,extra2);
  view_constraint_prog(t2,c,l,attr2,extra1,extra2);
}

/////////////////////////////////////////////////////////////////////////////////////////////////

void q2_bjorkenx_comp(TTree * t, TCanvas * c, vector<string> attr, string tmp){

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

  TH2D * comp = new TH2D(tmp.c_str(),"",100,0,4500,80,0.75,1.25);
  t->Draw(Form("bjorken_x_rec:Q2_true/1E3 >> %s",tmp.c_str()), selection.c_str()); //flag
  plot_title_style(comp,Form("bjorken_x vx Q^{2} (E_{%s} %s GeV)",inter.c_str(),Enu.c_str()),"Q^{2}","bjorken x");
  comp->Draw("colz");
  c->SaveAs(Form("../GENIEescatcomps/pngdump/bjorken_q2/%s_%s_%s_bjorken_q2_%s.png",inter.c_str(),Enu_t.c_str(),tgt.c_str(),sel_title.c_str()));
  c->Clear();

}

/////////////////////////////////////////////////////////////////////////////////////////////////

void q0_q3_comp(TTree * t, TCanvas * c, vector<string> attr, string tmp){

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

  TH2D * comp = new TH2D(tmp.c_str(),"",100,0,4500,80,0,4500);
  t->Draw(Form("q0_true:q3_true >> %s",tmp.c_str()), selection.c_str()); //flag
  plot_title_style(comp,Form("q_{0} v q_{3} (E_{%s} %s GeV)",inter.c_str(),Enu.c_str()),"q_{3}","q_{0}");
  comp->Draw("colz");
  c->SaveAs(Form("../GENIEescatcomps/pngdump/q0_q3/%s_%s_%s_q0_q3_%s.png",inter.c_str(),Enu_t.c_str(),tgt.c_str(),sel_title.c_str()));
  c->Clear();

}

/////////////////////////////////////////////////////////////////////////////////////////////////

void q2_q0_comp(TTree * t, TCanvas * c, vector<string> attr, string tmp){
  /*
  Makes a Q2 vs Q0 comparison given a tree and some sort of acceptance
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

  TH2D * comp = new TH2D(tmp.c_str(),"",100,0,4500,80,0,4500);
  t->Draw(Form("Q2_true/1E3:q0_true >> %s",tmp.c_str()), selection.c_str()); //flag
  plot_title_style(comp,Form("Q^{2} v q_{0} (E_{%s} %s GeV)",inter.c_str(),Enu.c_str()),"q_{0}","Q^{2}");
  comp->Draw("colz");
  c->SaveAs(Form("../GENIEescatcomps/pngdump/q2_q0/%s_%s_%s_q2_q0_%s.png",inter.c_str(),Enu_t.c_str(),tgt.c_str(),sel_title.c_str()));
  c->Clear();

}

/////////////////////////////////////////////////////////////////////////////////////////////////

void view_constraint_prog(TTree * t, TCanvas * c, TLegend * l, vector<string> attr, vector<string> sel, vector<string> sel_ttl){
  /*
  Shows multiple constraints fo a single graph given a vector of selections
  */

  string selection = attr[0];
  string sel_title = attr[1];
  string Enu = attr[2];
  string tgt = attr[3];
  string inter = attr[5];
  TH1D * temp[22];
  string Enu_t;
  double Enu_d;
  if(Enu == "2.261"){Enu_t = "2_2"; Enu_d = 2261;}
  else if(Enu == "3.3"){Enu_t = "3_3"; Enu_d = 3300;}
  else if(Enu == "4.4"){Enu_t = "4_4"; Enu_d = 4400;}
  else if(Enu == "1.1"){Enu_t = "1_1"; Enu_d = 1100;}

  for(int i = 0; i < sel.size(); ++i){
    attr[0] = sel[i];
    attr[1] = sel_ttl[i];
    temp[i] = e_frac(t,c,attr,Form("%s%s%s",Enu.c_str(),sel_ttl[i].c_str(),attr[4].c_str()));
  }
  plot_title_style(temp[0],"Ereco Constraint Progression","Energy Reconstructed","# of Events");
  temp[0]->SetLineWidth(2);
  temp[0]->SetLineColor(1);
  temp[0]->Draw();
  l->AddEntry(temp[0],sel_ttl[0].c_str());
  for(int i = 1; i < sel.size(); ++i){
    temp[i]->SetLineWidth(2);
    temp[i]->SetLineColor(i+1);
    temp[i]->Draw("same");
    l->AddEntry(temp[i],sel_ttl[i].c_str());
  }
  l->Draw("same");
  c->SaveAs(Form("../GENIEescatcomps/pngdump/efrac/%s_%s_%s_econt_mode_prog.png",inter.c_str(),Enu_t.c_str(),tgt.c_str()));
  c->Clear();
  l->Clear();
}

/////////////////////////////////////////////////////////////////////////////////////////////////

void view_ekin_pperp(TTree * tGENIE, TTree * tNEUT, TCanvas * c, TLegend * l, vector<string> attr, vector<string> sel, vector<string> sel_ttl, TFile * dat){
  /*
  Makes a comparison of Kinetic energy reconstructed vs the perpendicular momentum given a some required values
  */

  string selection = attr[0];
  string sel_title = attr[1];
  string Enu = attr[2];
  string tgt = attr[3];
  string inter = attr[5];
  TH1 * temp[22];
  string Enu_t;
  double Enu_d;
  if(Enu == "2.261"){Enu_t = "2_2"; Enu_d = 2261;}
  else if(Enu == "4.4"){Enu_t = "4_4"; Enu_d = 4400;}

  TPaveText* pt = new TPaveText(0.7,0.8,0.9,0.9);
//  pt->SetFillColor(18);
  pt->SetTextAlign(12);
  pt->AddText("^{12}C");
  cout << "-- Creating Graphs --\n";
  for(int i = 0; i < 2; ++i){
    attr[0] = sel[i];
    attr[1] = sel_ttl[i];
    temp[i] = e_temp2(tNEUT,c,attr,Form("N%s%s%s",Enu.c_str(),sel_ttl[i].c_str(),attr[4].c_str()));
    temp[i+2] = e_temp2(tGENIE,c,attr,Form("G%s%s%s",Enu.c_str(),sel_ttl[i].c_str(),attr[4].c_str()));
  }
  cout << "-- Setting Graph Attributes --\n";
  for(int i = 0; i < 2; ++i){
    temp[i]->SetLineWidth(2);
    temp[i+2]->SetLineWidth(2);
    temp[i]->SetLineColor(i+1);
    temp[i+2]->SetLineColor(i+1);
  }
  cout << "-- Drawing NEUT graphs --\n";
  temp[0]->Draw();l->AddEntry(temp[0],sel_ttl[0].c_str());
  temp[1]->Draw("same");l->AddEntry(temp[1],sel_ttl[1].c_str());
  l->Draw("same");
  c->SaveAs(Form("../GENIEescatcomps/test/NEUT_%s_%s_%s_ekin_pperp_prog.png",inter.c_str(),Enu_t.c_str(),tgt.c_str()));
  l->Clear();
  
  cout << "-- Drawing GENIE graphs --\n";
  temp[2]->Draw();l->AddEntry(temp[2],sel_ttl[0].c_str());
  temp[3]->Draw("same");l->AddEntry(temp[3],sel_ttl[1].c_str());
  l->Draw("same");
  c->SaveAs(Form("../GENIEescatcomps/test/GENIE_%s_%s_%s_ekin_pperp_prog.png",inter.c_str(),Enu_t.c_str(),tgt.c_str()));
  l->Clear();

  cout << "-- Adding GENIE to NEUT graphs --\n";
  temp[0]->SetLineColor(nBlue);
  temp[1]->SetLineColor(nOrange);
  temp[2]->SetLineColor(nYellow);
  temp[3]->SetLineColor(nVermillion);
  temp[3]->Scale(temp[1]->Integral() / temp[3]->Integral());
  temp[1]->Draw();
  l->AddEntry(temp[0],Form("NEUT %s",sel_ttl[0].c_str()));
  l->AddEntry(temp[1],Form("NEUT %s",sel_ttl[1].c_str()));
  l->AddEntry(temp[2],Form("GENIE %s",sel_ttl[0].c_str()));
  l->AddEntry(temp[3],Form("GENIE %s",sel_ttl[1].c_str()));
  temp[3]->SetLineStyle(2);
  temp[0]->Draw("same");
  temp[2]->Draw("same");
  temp[3]->Draw("same");
  l->Draw("same");
  c->SaveAs(Form("../GENIEescatcomps/test/GENIE_NEUT_%s_%s_%s_ekin_pperp_prog.png",inter.c_str(),Enu_t.c_str(),tgt.c_str()));

  TH1F * CLASdat;
  dat->GetObject("h_Erec_subtruct_piplpimi_noprot_frac_feed",CLASdat);
  CLASdat->Scale(temp[1]->Integral() / CLASdat->Integral());
  l->AddEntry(CLASdat,"CLAS data");
  CLASdat->Draw("same");
  c->SaveAs(Form("../GENIEescatcomps/test/GENIE_NEUT_CLAS_%s_%s_%s_ekin_pperp_prog.png",inter.c_str(),Enu_t.c_str(),tgt.c_str()));

  l->Clear();
}

/////////////////////////////////////////////////////////////////////////////////////////////////

#endif
