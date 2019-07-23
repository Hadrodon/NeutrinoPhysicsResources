/* 
Magdelena Allen | Summer 2018 SULI Program
Final version of macro for analysis. Generates distributions for kinematic variables for GENIE MC by reading in GST files, computing variables, applying weighting, and taking cuts. Adapted from Kevin Ewart's code (2017). 

Instructions for use: 
	- All configuration occurs on lines 102-129.
	- Change file name and path (both for reading in and saving) on lines 104 and 105. Note that .root/gst is added for you. (The save name is automatically generated from the file name, in the form "Histograms_MSA6_*_variable_weightingoption.root". To override this, change the save path/name on lines 567 and 585.)
	- Set final lepton mass and bound nuclear mass on lines 108 and 110 depending on your choice of scattering mode and target. 
	- Set weighting for events on line 112. Use 0 for no weighting and 4 to apply an electron-to-nu cross section correction. 
	- Apply cuts on lines 114-123 and set cutoffs on lines 125-129. Turn on the CC cut for nu runs and turn it off for electron events (or no events will be selected). The Q^2, W, x, and final particle cuts exist to isolate the QE channel in an experimentally reproducible way.
	- Descriptions of the output variables can be found on lines 62-72. 
*/


#include "TH1.h"
#include "TH3.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TVector3.h"

#include <iostream>
#include <sstream>

using namespace std;


// ---------- Function to find charge of a particle (Ewart) ----------

int findCharge(int pdg) {
	int apdg = abs(pdg);
	if(apdg<100) return 0; //skips bosons and other weird stuff, which are handled elsewhere
	else if(apdg>1000000000) return 0; //skip nuclei
	int q[3];
	q[0] = apdg%10000/1000;
    q[1] = apdg%1000/100;
    q[2] = apdg%100/10;
    int charge = 0;
    for(int i = 0;i<3;i++){
        if(q[i]%2==0 && q[i]!=0){
            charge+=2;
        	}
        else if(q[i]%2==1){
            charge+=1;
        	}
    	}
    switch(charge){
        case 2: return 0; //dd' quark antiquark or similar
        case 3: return pdg/apdg; //ud or ddd, charge determined by sign of pdg
        case 4: return 0; //uu' or ud'd'
        //Baryons
        case 5: return pdg/apdg; //uud'
        case 6: return 2*pdg/apdg; //uuu
        default: throw "Error! Invalid charge!";
    }
}


// ---------- Variables to go into histograms ----------

enum Vars// : int
{vE,                    //incoming neutrino beam energy
EQE, EQEres,            //quasielastic reconstruction, resolution
cal, calRes,            //calorimetric reconstuction, resolution
q0, 					//q0 (energy transfer = E_i,l - E_f,l)
q3,                 	//q3 (momentum transfer = sqrt(Q^2 + q0^2) with Q the four-momentum transfer)
lE, lP, lpT, lpz, lTheta, lCosTh, lPhi, //primary lepton energy, momentum, transverse momentum, longitudinal momentum, theta angle, cos(theta), phi angle
pP, pT, ppT, pTheta, pCosTh, pPhi,	    //primary proton momentum, kinetic energy, transverse momentum, theta angle, cos(theta), phi angle
pNum,                   //number protons produced
piP, nP, knP,			//pion and neutron momenta (GENIE), kinematic neutron momentum
DPT, DPHIT, DALPHAT,    //imbalance momentum, phi, and alpha (see X.-G. Lu et al. 2018)
Q2_ind, W_ind, bjx_ind, nu_ind,		    //distributions of kinematic variables
VARS_NUM}; 				//this value tracks number of items in enum and must be kept at the end

char const* varNames[] = 
{"vE", "EQE", "EQEres", "cal", "calRes", "q0", "q3", "lE", "lP", "lpT", "lpz", "lTheta", "lCosTh", "lPhi", "pP", "pT", "ppT", "pTheta", "pCosTh", "pPhi", "pNum", "piP", "nP", "knP", "dpt", "dphit", "dalphat", "Q2", "W", "bjx", "nu"};


// ---------- Cuts available perform ----------

enum Cuts// : int
{QE_fp_cut, //select for QE by requiring a single final-product particle as a proton (GENIE)
Q2_cut,		//Q2 > Q2_cutoff
p_cut, 		//events that have exactly one proton with momentum > 0.3 GeV/c and no other charged particles with momentum > p_cutoff
QEL_cut, 	//events that are quasi-elastic (GENIE)
CC_cut, 	//events that are charged current (if applicable)
pion_cut, 	//no produced pions with momentum > pion_cutoff (GeV)
W_cut, 		//invariant mass W < W_cutoff (GeV)
bjx_cut, 	//Bjorken x factor abs(x-1) < bjx_cutoff
RES_cut,	//events that are resonant (GENIE)
CUTS_NUM};

char const* cutNames[] = 
{"QE_fp_cut", "Q2_cut", "p_cut", "QEL_cut", "CC_cut", "pion_cut", "W_cut", "bjx_cut"};

bool cut_switches[CUTS_NUM] = {false};  //These flags will switch on which cuts to make on events
bool cut_flags[CUTS_NUM] = {false}; 	//These flags will be turned true if an event passes the cuts


void simHists_MSA6_multip(){

	// --------- CONFIGURE FILE INFORMATION AND PROCESSING OPTIONS ---------
	//Select file options: 
	const char* NAME = "nu_12C_hA_MC_gst"; //source file name and save name basis
	string PATH = "../genie_mc/nu_hA_MC/"; //path to source file and save directory
	//const string READNAME = string(NAME)+".root/gst"; //makes filename to access GST
        const string READNAME = "/pnfs/minerva/persistent/users/betan009/DUNE/genie_ND_carbon_T2K_gst.root/gst"; 
	//Select scattering mode: final lepton mass
	const Double_t m_l = 0.105658;  //0.000510999 electron mass (amu); 0.105658 muon mass (amu);
	//Select target: nuclear mass (used in kinematic neutron momentum calculation)
	const Double_t m_a = 11.17486; //bound nuclear mass (GeV); 11.17486 for 12C, 52.08977 for 56Fe, 37.55098 for 40Ar
	//Select weighting options
	int weighting = 0; //0=none, 1=CLAS, 2=CLAS and inverse Mottxsection, 3=inverse Mottxsection, 4 = e- to nu mode cross section correction (ratio of neutrino Mottxsec equivalent and Mottxsec - best correction for e-/nu comparision). NOTE: 1 & 2 accounting for CLAS acceptance map are NOT yet available! 
	//TFile* file_acceptance = new TFile("Map/eg2_maps_12C_E_4_461.root"); //weighting map (not available)
	//Cut options: true == cut is taken
	cut_switches[CC_cut] = true;	 //charged-current events only: turn on for nu, turn off for e- events (because nc==0 && cc==0 for them)! 
	cut_switches[Q2_cut] = true; 	 //Q2 > Q2_cutoff
	cut_switches[p_cut] = true; 	 //events that have exactly one proton with momentum > 0.3 GeV/c and no other charged particles with momentum > p_cutoff
	cut_switches[W_cut] = true;      //invariant mass W < W_cutoff (GeV)
	cut_switches[pion_cut] = true;	 //no produced pions with momentum > pion_cutoff (GeV)
	cut_switches[bjx_cut] = true;	 //Bjorken x factor abs(x-1) < bjx_cutoff
	cut_switches[QE_fp_cut] = false; //GENIE-given QEL channel by taking only single proton final particles---keep off
	cut_switches[QEL_cut] = false;   //GENIE-given QEL channel---keep off. Can turn on for CCQE run comparison (if so manually change save name to include _CCQE)
	cut_switches[RES_cut] = false;	 //only includes resonant channel
	
	const double Q2_cutoff = 0.5; //GeV. 0.5 GeV for 2.2 GeV analysis, 1.0 GeV for 4.4 GeV analysis
	const double p_cutoff = 0.3;  //GeV. final-state proton momentum cutoff
	const double pion_cutoff = 0; //GeV. 0 to exclude pions, 150 MeV for JLab pion momentum threshold
	const double W_cutoff = 2.0;  //GeV
	const double bjx_cutoff = 0.2;

	// --------- PHYSICAL CONSTANTS (from NIST unless otherwise specified) --------- 
	const Double_t Eb = 0.02; 		//binding energy (from Afro)
    const Double_t m_p =  0.938272; //mass of proton
    const Double_t m_n =  0.939565; //mass of neutron
    const Double_t m_pic = 0.13958; //mass of charged pion https://urldefense.proofpoint.com/v2/url?u=https-3A__doi.org_10.1103_PhysRev.163.1451&d=DwIGAg&c=gRgGjJ3BkIsb5y6s49QqsA&r=Qd3NpOP9HPHXo6y5ljG34X1gMsPSfzG6cZrTqOH1Bu4&m=XAsD64ome9rSXO9exFhQqN3TEosZnVc2PlePdP7zZYg&s=RfHD4kq1k3ID-TabZItdj1LxdtSIOfIapIuj0jY69-g&e=
    const Double_t m_pi0 = 0.13498; //(GeV/c^2) mass of neutral pion (wikipedia, couldn't find better)
    const Double_t fine_struc_const = 1./137.035999139;

	// --------- CONFIGURE 1D HISTOGRAMS (parameters same shape and order as Vars) --------- 
    //Bin numbers
    const int BINS[] =
        {150,
        80,80,
        100,100,
        80,80,
        80,40, 40, 40, 20, 80, 60,
        60,80, 40, 80,80, 60,
        9,
		20, 30, 30,
        20,24,9,
		20,20,20,20}; 
    //Histogram minimums and maximums
    const double HMIN[] = 
        {0.0,
        0.0, -0.3,
        0.0, -0.3,
        0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -90.,
        0.3, 0.0, 0.0, 0.0, -1.0, -90.,
        0,
		0., 0., 0.,
        0.0, 0.0, 0.0,
		0.0,0.0,0.75,0.0};
    const double HMAX[] = 
        {6.0,
        6.0, 0.3,
        5.0, 0.3,
        5.0, 5.0,
        5.0, 5, 2.0, 4.25, 180., 1., 90.,
        3.5, 5.0 , 2.0, 180., 1., 90.,
        9,
		0.2, 3., 3.,
        .7, 90., 180.,
		4.4,2.0,1.25,2.5};
	//Histogram titles
    const char* TITLES[] = 
        {"Incoming Lepton Energy",
        "Quasi-Elastic Reconstruction (E_QE)", "E_QE Resolution",
        "Calorimetric Reconstruction", "Calorimetric Reco Resolution",
        "q0", "q3",
        "Primary Lepton Energy", "Primary Lepton Momentum", "Primary Lepton Transverse Momentum", "Primary Lepton Longitudinal Momentum", "Primary Lepton Theta", "Primary Lepton cos(Theta)", "Primary lepton Phi",
        "Proton Momentum", "Proton Kinetic Energy", "Primary Proton Transverse Momentum", "Proton Angle", "Proton cos(Theta)", "Proton Phi",
        "Protons Produced",
		"Pion Momentum", "GENIE neutron momentum", "Kinematic neutron momentum",
        "Delta Momentum Magnitude", "Delta Phi", "Delta Alpha",
		"Momentum Transfer Q^2", "Invariant Mass W", "Bjorken x", "Kinematic Nu"};
	//Make histograms
	TH1D *h[VARS_NUM];
	for(int i=0; i<VARS_NUM; i++) {
		string name = string("h") + varNames[i]; //+ "_" + cutNames;
		h[i] = new TH1D(name.c_str(), TITLES[i], BINS[i], HMIN[i], HMAX[i]);
	}

	// --------------- CONFIGURE 2D HISTOGRAMS -----------------
	gStyle->SetStripDecimals(false);
    TGaxis::SetMaxDigits(3);
    //gStyle->SetOptStat(0);
    gStyle->SetPadLeftMargin(0.15);	
	gStyle->SetPadBottomMargin(0.15);
	TH2D * h2D = new TH2D("h2D", "", 200, 0, 3, 200, 0, 3);
	h2D->SetStats(false);

	// --------- SETUP GST --------- 
	TChain *tree = new TChain(); //for whole-directory processing
    tree->Add((READNAME).c_str(),0); //get tree
	Long64_t nentries = tree->GetEntries(); //For troubleshooting, can set to 10 events
	//Find maximum number of final hadronic products produced in TChain (GST) of events (i.e., max value of nf branch)
	tree->SetBranchStatus("*", 0); //Turn off all the branches
    tree->SetBranchStatus("nf",1); //Turn on nf branch (number final-state particles)...
	Int_t nf;
	int maxNF = 0; //counter
    tree->SetBranchAddress("nf",&nf);
	// Loop through entries to find max
    for(Long64_t ij=0;ij<nentries;ij++){
        tree->GetEntry(ij);   
        if(nf>maxNF){
            maxNF = nf;
        }
    }
    cout<<"maxNF: "<<maxNF<<"\n";
	//Initialize variables for GST
    Double_t pxl; //final lepton x momentum (GeV)
    Double_t pyl; //final lepton y momentum (GeV)
    Double_t pzl; //final lepton z momentum (GeV)
    Double_t El;  //final lepton energy (GeV)
	Int_t nfp;	  //number final p and pbar
	Int_t nfn;	  //number final n and nbar
	Int_t nfpi0;  //number final pi0
	Int_t nfpip;  //number final pi+
	Int_t nfpim;  //number final pi-
    Double_t Ev;  //initial lepton (beam) energy (GeV)
    Double_t Q2;  //4-momentum transfer (GeV)
    Double_t x;   //Bjorken x (calculated from event)
    Double_t W;   //invariant mass (calculated from event) (GeV)
    Int_t neu;    //neutrino pdg
	Bool_t qel;   //Quasi-elastic channel?
	Bool_t mec;	  //MEC channel?
	Bool_t res;   //Resonance channel? 
	Bool_t dis;   //DIS channel?
	Bool_t cc;	  //Charged current?
	Bool_t nc;    //Neutral current? 
    //Initialize arrays for GST final-particle information storage
    /*Double_t Ef[maxNF];
    Int_t pdgf[maxNF];
    Double_t pxf[maxNF], pyf[maxNF], pzf[maxNF];
    */

    Double_t *Ef=new Double_t[maxNF];
    Int_t *pdgf=new Int_t[maxNF];
    Double_t *pxf= new Double_t[maxNF];
    Double_t *pyf= new Double_t[maxNF];
    Double_t *pzf= new Double_t[maxNF];
    Double_t vtxx, vtxy, vtxz;
	//Turn (back) on all used branches
    tree->SetBranchStatus("El", 1);
    tree->SetBranchStatus("Ef", 1);
    tree->SetBranchStatus("pdgf",1);
    tree->SetBranchStatus("pxl",1);
    tree->SetBranchStatus("pyl",1);
    tree->SetBranchStatus("pzl",1);
    tree->SetBranchStatus("qel",1);
    tree->SetBranchStatus("mec",1);
    tree->SetBranchStatus("res",1);
    tree->SetBranchStatus("dis",1);
    tree->SetBranchStatus("Ev",1);
    tree->SetBranchStatus("Q2",1);
    tree->SetBranchStatus("pxf",1);
    tree->SetBranchStatus("pyf",1);
    tree->SetBranchStatus("pzf",1);
    tree->SetBranchStatus("cc",1);
	tree->SetBranchStatus("nc",1);
    tree->SetBranchStatus("x",1);
    tree->SetBranchStatus("W",1);
    tree->SetBranchStatus("neu",1);
	tree->SetBranchStatus("nfp",1);
	tree->SetBranchStatus("nfn",1);
	tree->SetBranchStatus("nfpi0",1);
	tree->SetBranchStatus("nfpim",1);
	tree->SetBranchStatus("nfpip",1);
	//Set branch addresses
    tree->SetBranchAddress("El",&El);
    tree->SetBranchAddress("Ef",&Ef);     //energy of [kth] final state particle in hadronic system (GeV)
    tree->SetBranchAddress("pdgf",&pdgf); //PDG code of [kth] final state particle in hadronic system
    tree->SetBranchAddress("pxl",&pxl);
    tree->SetBranchAddress("pyl",&pyl);
    tree->SetBranchAddress("pzl",&pzl);
	tree->SetBranchAddress("qel",&qel);
    tree->SetBranchAddress("mec",&mec);
    tree->SetBranchAddress("res",&res);
    tree->SetBranchAddress("dis",&dis);
    tree->SetBranchAddress("Ev",&Ev);
    tree->SetBranchAddress("Q2",&Q2);
    tree->SetBranchAddress("pxf",&pxf); //Px of [kth] final state particle in hadronic system (GeV)
    tree->SetBranchAddress("pyf",&pyf); //Px of [kth] final state particle in hadronic system (GeV)
    tree->SetBranchAddress("pzf",&pzf); //Px of [kth] final state particle in hadronic system (GeV)
    tree->SetBranchAddress("cc",&cc);
	tree->SetBranchAddress("nc",&nc);
    tree->SetBranchAddress("x",&x);
    tree->SetBranchAddress("W",&W);
    tree->SetBranchAddress("neu",&neu);
	tree->SetBranchAddress("nfp",&nfp);
	tree->SetBranchAddress("nfn",&nfn);
	tree->SetBranchAddress("nfpi0",&nfpi0);
	tree->SetBranchAddress("nfpip",&nfpip);
	tree->SetBranchAddress("nfpim",&nfpim);

	// --------- ITERATE OVER ALL EVENTS --------- 

    for(Long64_t ill=0;ill<nentries;ill++) {
		//Progress bar
       /* if(i%(nentries/100)==0){
            cout<<"\r";
            cout<<"Filling histograms... "<<(i/(nentries/100))<<"%";
            cout.flush();
        }
        */
        //Fill pointers
        tree->GetEntry(ill);

		//MAKE CUTS: determine if the ith event passes available cuts

		cut_flags[QE_fp_cut] = (nf==1 && pdgf[0]==2212);
		cut_flags[Q2_cut] = (Q2 > Q2_cutoff);
		cut_flags[QEL_cut] = (qel);
		cut_flags[RES_cut] = (res);
		cut_flags[CC_cut] = (cc);
		cut_flags[W_cut] = (W < W_cutoff);
		cut_flags[bjx_cut] = (TMath::Abs(x-1) < bjx_cutoff);

		//Iterate over all final-state particles: p_cut, pion_cut

		cut_flags[p_cut] = true;    //passes unless proved otherwise
		cut_flags[pion_cut] = true; //passes unless proved otherwise
		
		int p_count_cut = 0;  //count number of final-state protons in event with protonMomentum > p_cutoff
		int p_count = 0;      //count for number of final-state protons in event
		int pi_count = 0;     //count number of final pions (+,-,0)
		int n_count = 0;      //count number of final neutrons
		//Attributes of leading (most energetic) final-state proton
		int ppi = -1; 		  //leading proton index in the nfp-length p_count array of p+pbars
		int protonIndex = -1; //leading proton index in nf-length GST arrays of final particles
		double protonTotE = -1;
		double protonKE = -1;
		double protonMomentum = -1;
		double protonTheta = -1;
		double protonCosTheta = -1;
		double protonPhi = -1;
		double protonTransverseMomentum = -1;
		//Attributes for p+ production histograms
		/*double pMomenta[nfp];
		int pIndices[nfp];
		double pKEs[nfp];
		double pThetas[nfp];
		double pCosThetas[nfp];
		double pPhis[nfp];
		double pTransverseMomenta[nfp];
                 */
                Double_t *pMomenta= new Double_t[nfp];
                Double_t *pKEs= new Double_t[nfp];
                Double_t *pThetas= new Double_t[nfp];
                Double_t *pCosThetas= new Double_t[nfp];
                Double_t *pPhis= new Double_t[nfp];
                Double_t *pTransverseMomenta= new Double_t[nfp];
                Int_t *pIndices= new Int_t[nfp];
		//Attributes for n0, pi+/-/0 histograms
		int n_pis = nfpi0 + nfpip + nfpim; //total number of pions produced
		//double piPs[n_pis];  //final-state pion momenta
		//double nPs[nfn];	 //final-state neutron momenta
                Double_t *piPs= new Double_t[n_pis];  //final-state pion momenta
                Double_t *nPs= new Double_t[nfn];         //final-state neutron momenta

		for(Int_t k=0; k<nf; k++) { //loop through k final particles
	
			//Calculations for kth final-state particle:
			int pdg = pdgf[k];	
			double kMomentum = TMath::Sqrt(TMath::Power(pxf[k],2)+TMath::Power(pyf[k],2)+TMath::Power(pzf[k],2));
            double kTheta = TMath::ATan2(TMath::Sqrt(TMath::Power(pxf[k],2)+TMath::Power(pyf[k],2)),pzf[k])*180/TMath::Pi(); //deg
            double kCosTheta = pzf[k]/kMomentum; //TMath::Cos(kTheta*TMath::Pi()/180) checks out
            double kPhi = TMath::ATan2(pyf[k],pxf[k]);
			double kTransverseMomentum = TMath::Sqrt(TMath::Power(pxf[k],2)+TMath::Power(pyf[k],2));

			//Final-state particle cuts:
			//charged, non-proton particle with momentum exceeding p_cutoff
			if((abs(pdg)==0 || abs(pdg) == 11 || abs(pdg)==13 || abs(pdg)==15 || abs(pdg)==24 || abs(pdg)==34 || abs(pdg)==37 || (findCharge(abs(pdg))!=0 && pdg!=2212)) && kMomentum>p_cutoff) {
				cut_flags[p_cut] = false; 
			}
			//pions
			if(abs(pdg)==111 || abs(pdg)==211) { //pi+, pi-, and pi0
				if (kMomentum > pion_cutoff) { //pion cut
					cut_flags[pion_cut] = false;
				}
				piPs[pi_count] = kMomentum;
				pi_count++;
			}
			//neutrons
			if(abs(pdg) == 2112) { //neutrons and antineutrons
				nPs[n_count] = kMomentum;
				n_count++;
			}
			//protons
			if(pdg == 2212){
				if(kMomentum > p_cutoff) { //assume primary proton
					if (kMomentum > protonMomentum) { //get most primary (most energetic) proton
						protonIndex = k;  //index of leading proton in GENIE's k final state particles
						protonKE = Ef[k]-m_p;
						protonMomentum = kMomentum;
						protonTheta = kTheta;
						protonCosTheta = kCosTheta;
						protonPhi = kPhi;
						protonTransverseMomentum = kTransverseMomentum;
						ppi = p_count;  //primary proton index in nfp final p/pbars
					}
					p_count_cut++; //count all passing protons
				}
				pMomenta[p_count] = kMomentum;
				pIndices[p_count] = k;
				pKEs[p_count] = Ef[k]-m_p;
				pThetas[p_count] = kTheta;
				pCosThetas[p_count] = kCosTheta;
				pPhis[p_count] = kPhi;
				pTransverseMomenta[p_count] = kTransverseMomentum;
				p_count++;
			}else if(pdg == -2212) { //pbars counted in nfp, but we don't want them in the distributions
				pMomenta[p_count] = -99; 
				p_count++;
			}
		} //end loop through final-state particles

		if(p_count_cut != 1) cut_flags[p_cut] = false; //p_cut (EXACTLY one proton with momentum > p_cutoff)

		// --------- CALCULATIONS FOR HISTOGRAMS --------- 
		//final lepton information
        double leptonMomentum = TMath::Sqrt(TMath::Power(pxl,2)+TMath::Power(pyl,2)+TMath::Power(pzl,2));
        double leptonTransverseMomentum = TMath::Sqrt(TMath::Power(pxl,2)+TMath::Power(pyl,2));
        double leptonTheta = TMath::ATan2(TMath::Sqrt(TMath::Power(pxl,2)+TMath::Power(pyl,2)),pzl)*180/TMath::Pi();
        double leptonCosTheta = pzl/leptonMomentum; //TMath::Cos(leptonTheta*TMath::Pi()/180) checks out
        double leptonPhi = TMath::ATan2(pyl,pxl);
		//(leading proton calculations already done in final-state particle loop)
		//kinematics
        double E_QE = (pow(m_n,2)-pow(m_p-Eb,2)-pow(m_l,2)+2*(m_p-Eb)*El)/(2*(m_p-Eb-El+leptonMomentum*leptonCosTheta));
        double E_QEres = (Ev-E_QE)/Ev;
        double q0f = Ev-El;
        double q3f = TMath::Sqrt(TMath::Power(Ev-El,2)+Q2);
		double calf = El + protonKE + Eb;
        double calResf = (Ev-calf)/Ev;
		double nu = Ev - El; //kinematic variable nu = energy lost by incoming lepton (in lab frame)
		//imbalance vars (for each energetic proton)
        TVector3 ptmuon = *(new TVector3(pxl,pyl,0)); //xy planar projection to get transverse momentum vector
		/*double dpt[p_count];
		double dphit[p_count];
		double dalphat[p_count];
               */
               Double_t *dpt= new Double_t[p_count];
                Double_t *dphit= new Double_t[p_count];
                Double_t *dalphat= new Double_t[p_count];
		for(int fi=0; fi<p_count; fi++) {
			if(pMomenta[fi] > 0) { //screen out pbars with -99 flag
				TVector3 ptproton = *(new TVector3(pxf[pIndices[fi]],pyf[pIndices[fi]],0));
				TVector3 tmpd = ptmuon + ptproton;
				dpt[fi] = tmpd.Mag(); //GeV
				dphit[fi] = TMath::ACos(ptmuon.Dot(ptproton)*(-1)/ptmuon.Mag()/ptproton.Mag())*TMath::RadToDeg(); //Deg
				dalphat[fi] = TMath::ACos(tmpd.Dot(ptmuon)*(-1)/(tmpd.Mag()*ptmuon.Mag()))*TMath::RadToDeg();     //Deg
			}else {
				dpt[fi] = -99;
				dphit[fi] = -99;
				dalphat[fi] = -99;
			}
		}
		//kinematic neutron momentum
		double m_astar = m_a - m_n + Eb; //GeV
		TVector3 protonMomentumVector = *(new TVector3(pxf[protonIndex],pyf[protonIndex],pzf[protonIndex]));  //full vector
		TVector3 leptonMomentumVector = *(new TVector3(pxl,pyl,pzl));
		double protonE = TMath::Sqrt(m_p*m_p + protonMomentumVector.Mag2());
		double leptonE = TMath::Sqrt(m_l*m_l + leptonMomentumVector.Mag2());
		double factor = m_a - protonE - leptonE + pzf[protonIndex] + pzl;
		double nPT = dpt[ppi];  //dpt of leading proton
		double nPL = -1*(m_astar*m_astar + nPT*nPT - factor*factor)/(2.*factor);
		double kinematicNeutronMomentum = TMath::Sqrt(nPL*nPL + nPT*nPT);

		// --------- WEIGHTING ---------  
		// Constants don't matter because these will be renormalized by integral during plotting (ComparePlotter.C)!
        double Mott_cross_sec = (TMath::Power(fine_struc_const,2)*(leptonCosTheta+1))/TMath::Power((2*TMath::Power(El,2)*(1-leptonCosTheta)),2); //Original Ewart macro
		//double Mott_cross_sec = (TMath::Power(fine_struc_const,2)*(leptonCosTheta+1))/(2*El*El*TMath::Power((1-leptonCosTheta),2)); //Thompson version? BUGGY
        // Acceptance weighting: not currently functional
        //double weight[2] = {1., 1.}; //default no weighting (e.g. weighting == 0)
        //double e_acc_ratio, p_acc_ratio;
		//if(weighting == 1 || weighting == 2) { //CLAS weighting
		//	e_acc_ratio = acceptance_c(leptonMomentum, leptonCosTheta, leptonPhi, 11, file_acceptance);
        //    p_acc_ratio = acceptance_c(protonMomentum, protonCosTheta, protonPhi, 2212, file_acceptance);
		//	weight[0] *= e_acc_ratio; //electron only
        //    weight[1] *= e_acc_ratio*p_acc_ratio;
		//	if(fabs(weight[0])!=weight[0] || fabs(weight[1])!=weight[1]) continue; // I [Afro] need this line because there are cases where I get infinity
		//}
		//if(weighting == 2 or weighting == 3) { //Mottxsection weighting
        //    //Define weights
        //    weight[0] /= Mott_cross_sec; //electron only
        //    weight[1] /= Mott_cross_sec; //proton and electron
        //}
		double weight = 1.;
		if(weighting == 3) { //weight only by Mottxsection
			weight /= Mott_cross_sec;
		}
		if(weighting == 4) { //weight by e- to nu mode cross section correction (Adi)
			weight *= Q2*Q2;
		}

		// --------- FILL HISTOGRAMS --------- 
		//Determine if event passed all implemented cuts
		bool pass = true;
		for(int k = 0; k < CUTS_NUM; k++) {
			if(cut_switches[k]==true && cut_flags[k]==false) {
				pass = false;
			}
		}
		//Fill histograms with passing events
		if(pass == true) {
			h[vE]->Fill(Ev,weight);
		    h[EQE]->Fill(E_QE,weight);
		    h[EQEres]->Fill(E_QEres,weight);
		    h[q0]->Fill(q0f,weight);
		    h[q3]->Fill(q3f,weight);
		    h[lE]->Fill(El,weight);
		    h[lP]->Fill(leptonMomentum,weight);
		    h[lTheta]->Fill(leptonTheta,weight);
		    h[lCosTh]->Fill(leptonCosTheta,weight);
			h[lPhi]->Fill(leptonPhi,weight);
		    h[lpT]->Fill(leptonTransverseMomentum,weight);
		    h[lpz]->Fill(pzl,weight);
	        h[cal]->Fill(calf,weight);
	        h[calRes]->Fill(calResf,weight);
			h[Q2_ind]->Fill(Q2,weight);
			h[W_ind]->Fill(W,weight);
			h[bjx_ind]->Fill(x,weight);
			h[nu_ind]->Fill(nu,weight);
			h[pNum]->Fill(p_count,weight); //counts all produced protons
			//Vars for leading proton
			//h[pP]->Fill(protonMomentum,weight);
			//h[pT]->Fill(protonKE,weight);
			//h[pTheta]->Fill(protonTheta,weight);
			//h[pCosTh]->Fill(protonCosTheta,weight);
			//h[pPhi]->Fill(protonPhi,weight);
			//h[ppT]->Fill(protonTransverseMomentum,weight);
			//h[DPT]->Fill(dpt[ppi],weight);
	    	//h[DPHIT]->Fill(dphit[ppi],weight);
	    	//h[DALPHAT]->Fill(dalphat[ppi],weight);
			h[knP]->Fill(kinematicNeutronMomentum,weight);
			//Distributions for all final-state particles:
			for (int fi=0; fi<pi_count; fi++) {
				if (piPs[fi] > 0) {  //screens out any errors set to 0 or -99
					h[piP]->Fill(piPs[fi],weight);
				}
			}
			for (int fi=0; fi<n_count; fi++) {
				if (nPs[fi] > 0) {  //screens out any errors set to 0 or -99
					h[nP]->Fill(nPs[fi],weight);
				}
			}
			//* If we want to include vars for ALL protons produced: 
			for (int fi=0; fi<p_count; fi++) {
				if (pMomenta[fi] > 0) { //screen out pbars with the -99 flag
					h[pP]->Fill(pMomenta[fi],weight);
					h[pT]->Fill(pKEs[fi],weight);
					h[pTheta]->Fill(pThetas[fi],weight);
					h[pCosTh]->Fill(pCosThetas[fi],weight);
					h[pPhi]->Fill(pPhis[fi],weight);
					h[ppT]->Fill(pTransverseMomenta[fi],weight);
				}
			}
			for (int fi=0; fi<p_count; fi++) {
				if(dpt[fi] > 0) { //screen out pbars with the -99 flag
					h[DPT]->Fill(dpt[fi],weight);
			    	h[DPHIT]->Fill(dphit[fi],weight);
			    	h[DALPHAT]->Fill(dalphat[fi],weight);
				}
			}
			//*/
			h2D->Fill(q0f, q3f, weight); //fill 2D histogram
		}

	} //end main for-loop through events

	// --------- SAVE FILES --------- 
	cout<<"\rFilling histograms... 100%\n"; //finish progress bar
	//set file name (notes weighting)
    string saveName = string(NAME);
    std::stringstream ss;
    ss << weighting;
    std::string ws = ss.str();
    saveName+=string("_weighted")+ws;
    string fileName = PATH+"Histograms_MSA6_"+saveName.c_str()+".root";
    //save 1D histograms
    TFile *fwrite = new TFile(fileName.c_str(),"recreate");
    for(int ii=0;ii<VARS_NUM;ii++){
        if(h[ii]->GetEntries() != 0) h[ii]->Write();
    }
	//draw 2D histograms
	h2D->SetXTitle("q0 (GeV)");
	h2D->GetXaxis()->SetTitleSize(0.06);
	h2D->GetXaxis()->CenterTitle();
	h2D->GetXaxis()->SetLabelSize(0.04);
	h2D->GetYaxis()->SetTitle("q3 (GeV)");
	h2D->GetYaxis()->SetTitleSize(0.06);
	h2D->GetYaxis()->CenterTitle();
	h2D->GetYaxis()->SetLabelSize(0.04);
	//h2D->Write();  //writes to the root file but in black and white
	TCanvas * ch2D = new TCanvas("ch2D", "ch2D", 1000, 750);
    h2D->Draw("ch2D,colz");
	ch2D->Print((PATH+"Histograms_MSA6_"+saveName+"_q0vsq3.png").c_str());
	ch2D->Close();

    cout<<"Done! Saved to "<<fileName<<"\n";

}

