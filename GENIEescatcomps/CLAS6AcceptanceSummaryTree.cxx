#include "Measurement1D.h"
#include "MeasurementVariableBox1D.h"
#include "PhysConst.h"
#include "CLAS6-EG2_Accepter.cxx"

static double const MProt = PhysConst::mass_proton * 1E3;
 
double Q2(TLorentzVector EPTransfer) { return -EPTransfer.Mag2(); }

double bjorkenX(TLorentzVector EPTransfer) {
  return Q2(EPTransfer) / (2 * MProt * EPTransfer.E());
}
 
class CLAS6AcceptanceSummaryTree : public Measurement1D {

  /// TTree used to dump out event-level variables
  TTree *eventVariables;
 
  /// Branch-backing data members.

  double FSLepton_mom, FSProton_mom, FSLepton_ct, FSProton_ct, FSLepton_phi, EFSProt, EFSLep, EFSPion, EFSNeutron, EFSOther, Pperp, EFSLepTrue;

  int NTrueFSChargedLep, NTrueProtons, NTruePions, NTrueNeutrons,
  NTrueOther, NEUTMode;
 
  int FSCLep_PDG, OtherPDG;

  double bjorken_x_rec, Q2_true, q0_true, q3_true, XSec_weight;
 
  double EISLepTrue, EISLepRec_calo, EISLepRec_KineFSLep, EISLepRec_KineFSLep_afro; //new

  double Eb, MottXSec, MottXSecRecip, W; //new

  CLASAccepter *clasAccepter;
 
  public:
    CLAS6AcceptanceSummaryTree(nuiskey samplekey) {
      clasAccepter = nullptr;

      fSettings = LoadSampleSettings(samplekey);
      fSettings.SetXTitle("");
      fSettings.SetYTitle("");
      fSettings.SetAllowedTypes("FIX/DIAG", "FIX,FREE,SHAPE/DIAG");
      fSettings.SetEnuRange(0.0, 100.0);
      fSettings.SetTitle("CLAS6AcceptanceSummaryTree");
      FinaliseSampleSettings();
     
      // Scale to xsec, may not be necc.
      fScaleFactor = (GetEventHistogram()->Integral("width") *1E-38) / (fNEvents + 0.) * 2. / 1.;
     
      /// Leave because important for NUSIANCE to not cry, but wont affect this
      /// sample
      fDataHist = new TH1D(("empty_data"), ("empty-data"), 100, 0, 10);
      SetupDefaultHist();
      SetCovarFromDiagonal();
      FinaliseMeasurement();

      Eb = 25; //new
      if (samplekey.Has("bindingenergy")){Eb = samplekey.GetD("bindingenergy");} //new

      // for (nuiskey &ch : samplekey.GetListOfChildNodes()) {
        // if (ch.GetElementName() == "CLASAccepter") {
          // std::cout << "[INFO]: Setting up CLASAccepter." << std::endl;
          // clasAccepter = new CLASAccepter();
          // clasAccepter->Setup(ch);
          // break;
        // }
      // }
    
      // Get output file so that we know this tree will live there
      Config::Get().out->cd();
      eventVariables = new TTree("CLAS6AcceptSummary", "");
      eventVariables->Branch("Pperp",&Pperp,"Pperp/D");
      eventVariables->Branch("FSLepton_mom", &FSLepton_mom, "FSLepton_mom/D");
      eventVariables->Branch("FSProton_mom", &FSProton_mom, "FSProton_mom/D");
      eventVariables->Branch("FSLepton_ct", &FSLepton_ct, "FSLepton_ct/D");
      eventVariables->Branch("FSProton_ct", &FSProton_ct, "FSProton_ct/D");
      eventVariables->Branch("FSLepton_phi", &FSLepton_phi, "FSLepton_phi/D");
      eventVariables->Branch("NTrueFSChargedLep", &NTrueFSChargedLep,"NTrueFSChargedLep/I");
      eventVariables->Branch("NTrueProtons", &NTrueProtons, "NTrueProtons/I");
      eventVariables->Branch("NTruePions", &NTruePions, "NTruePions/I");
      eventVariables->Branch("NTrueNeutrons", &NTrueNeutrons, "NTrueNeutrons/I");
       
      eventVariables->Branch("NTrueOther", &NTrueOther, "NTrueOther/I");
      eventVariables->Branch("OtherPDG", &OtherPDG, "OtherPDG/I");
      eventVariables->Branch("FSCLep_PDG", &FSCLep_PDG, "FSCLep_PDG/I");
      eventVariables->Branch("NEUTMode", &NEUTMode, "NEUTMode/I");
      eventVariables->Branch("bjorken_x_rec", &bjorken_x_rec, "bjorken_x_rec/D");
      eventVariables->Branch("Q2_true", &Q2_true, "Q2_true/D");
      eventVariables->Branch("q0_true", &q0_true, "q0_true/D");
      eventVariables->Branch("q3_true", &q3_true, "q3_true/D");
      eventVariables->Branch("EISLepTrue", &EISLepTrue, "EISLepTrue/D");
      eventVariables->Branch("EFSLepTrue", &EFSLepTrue, "EFSLepTrue/D");
      eventVariables->Branch("EISLepRec_calo", &EISLepRec_calo,"EISLepRec_calo/D");
      eventVariables->Branch("EISLepRec_KineFSLep", &EISLepRec_KineFSLep,"EISLepRec_KineFSLep/D");
      eventVariables->Branch("EISLepRec_KineFSLep_afro", &EISLepRec_KineFSLep_afro,"EISLepRec_KineFSLep_afro/D"); //new
      eventVariables->Branch("EFSProt", &EFSProt, "EFSProt/D");
      eventVariables->Branch("EFSLep", &EFSLep, "EFSLep/D");
      eventVariables->Branch("EFSPion", &EFSPion, "EFSPion/D");
      eventVariables->Branch("EFSNeutron", &EFSNeutron, "EFSNeutron/D");
      eventVariables->Branch("EFSOther", &EFSOther, "EFSOther/D");
      eventVariables->Branch("XSec_weight",&XSec_weight,"XSec_weight/D");
      eventVariables->Branch("MottXSec", &MottXSec, "MottXSec/D"); //new
      eventVariables->Branch("MottXSecRecip", &MottXSecRecip, "MottXSecRecip/D"); //new
      eventVariables->Branch("W",&W,"W/D"); //new
    }

    virtual ~CLAS6AcceptanceSummaryTree(){};
 
    //! Grab info from nvect
    void FillEventVariables(FitEvent *nvect) {
      FitParticle *hmfs_lep = nvect->GetHMFSLeptons();
      FitParticle *hmfs_prot = nvect->GetHMFSProton();
      if (!hmfs_prot || !hmfs_lep) {
        return;
      }
      
      NEUTMode = nvect->Mode;
      FitParticle *islep = nvect->GetHMISAnyLeptons();
      EISLepTrue = islep->E();
      EFSLepTrue = hmfs_lep->E(); //new
      
      TLorentzVector FourMomentumTransfer = (islep->P4() - hmfs_lep->P4());

      W = FitUtils::Wrec(islep->P4(),hmfs_lep->P4());

      FSLepton_mom = hmfs_lep->p();
      FSProton_mom = hmfs_prot->p();
      FSLepton_ct = hmfs_lep->P3().CosTheta();
      FSLepton_phi = hmfs_lep->P3().Phi();
      FSProton_ct = hmfs_prot->P3().CosTheta();
      Q2_true = Q2(FourMomentumTransfer);
      q0_true = FourMomentumTransfer.E();
      q3_true = FourMomentumTransfer.Vect().Mag();
      XSec_weight = nvect->InputWeight;
      bjorken_x_rec = bjorkenX(FourMomentumTransfer);

      constexpr double FSC2 = pow(0.0072973525664, 2); //new
      double A = FSC2 * (FSLepton_ct + 1.0);  //new
      double B = 2.0 * pow(EFSLepTrue, 2.0) * pow((1.0 - FSLepton_ct), 2.0); //new

      MottXSec = A/B;
      MottXSecRecip = 1.0 / MottXSec;
      FSCLep_PDG = hmfs_lep->PDG();
      EFSLep = hmfs_lep->E();
      EISLepRec_KineFSLep = FitUtils::EnuQErec(hmfs_lep->p() * 1E-3, hmfs_lep->P3().CosTheta(), Eb, (islep->PDG() > 0)) * 1E3; //changed
      EISLepRec_KineFSLep_afro = ((2.0 * MProt * Eb) + (2.0 * MProt * EFSLepTrue) - pow(PhysConst::mass_muon * 1E3, 2.0)) / (2.0 * (MProt - EFSLepTrue + FSLepton_mom * FSLepton_ct)); //new

      TVector3 FSLepProtMom(0, 0, 0); //new

      NTrueNeutrons = 0;
      NTrueProtons = 0;
      NTruePions = 0;
      NTrueFSChargedLep = 0;
      NTrueOther = 0;
      EISLepRec_calo = 0;
      EFSLep = 0;
      EFSProt = 0;
      EFSPion = 0;
      EFSNeutron = 0;
      EFSOther = 0;
      OtherPDG = 0;

      // Loop over particle stack
      for (size_t p_it = 0; p_it != nvect->NParticles(); ++p_it) {
        FitParticle *p = nvect->GetParticle(p_it);
        if (p->Status() != kFinalState) { // Only care about FS particles
          continue;
        }
        switch (abs(p->PDG())) {
          case 11:
          case 13: {
            EISLepRec_calo += p->E();
            FSLepProtMom += p->P3(); //new
            EFSLep += p->E();
            NTrueFSChargedLep++;
            break;
          }
          case 2212: {
            EISLepRec_calo += p->KE();
            EFSProt += p->KE();
            FSLepProtMom += p->P3();
            NTrueProtons++;
            break;
          }
          case 211:
          case 111: {
            //EISLepRec_calo += p->E(); removed
            //EISLepRec_KineFSLep += p->E(); removed
            EFSPion += p->E();
            NTruePions++;
            break;
          }
          case 2112: {
            EFSNeutron += p->KE();
            NTrueNeutrons++;
            break;
          }
          default: {
            //EISLepRec_calo += p->E(); removed
            //EISLepRec_KineFSLep += p->E(); removed
            OtherPDG = p->PDG();
            EFSOther += p->E();
            NTrueOther++;
          }
        }
      }
      Pperp = TVector3(FSLepProtMom[0], FSLepProtMom[1], 0).Mag(); //new
      eventVariables->Fill();
    };
     
  //! Define this samples signal
  bool isSignal(FitEvent *nvect) {
    FitParticle *hmfs_lep = nvect->GetHMFSLeptons();
    FitParticle *hmfs_prot = nvect->GetHMFSProton();
    
    if (!hmfs_prot || !hmfs_lep) {
      return false;
    }
    return true;
  }
 
  void FillHistograms() {}

  double GetLikelihood(void) { return 0; }
 
  void Write(std::string drawOpt) {
    Config::Get().out->cd();
    eventVariables->Write();
  }
};

