
#include <cmath>
#include <cstdlib>
#include <ostream>
#include <fstream>
#include <string>
#include <vector>

// ROOT include files.
#include "TTree.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TObject.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TFile.h"
#include "TH2.h"

using namespace std;


int main(int argc,char **argv){
  TApplication *ta=new TApplication("ta",&argc,argv);

  //===========================================================
  //Extracting the root files for analysis
  //===========================================================
  // K-40 Lines t1/2=1.251e9
  TFile k40("/Users/akindele1/Software/PMT/VetoStudy/ROOT_Files/twelve_std_glass/K40Bon.root");
  TTree *K40Data = (TTree*)k40.Get("data");

  //Th-232 Lines t1/2=1.405e10
  TFile ac228("/Users/akindele1/Software/PMT/VetoStudy/ROOT_Files/twelve_std_glass/Ac228Bon.root");
  TTree *Ac228Data = (TTree*)ac228.Get("data");

  TFile pb212("/Users/akindele1/Software/PMT/VetoStudy/ROOT_Files/twelve_std_glass/pb212Bon.root");
  TTree *Pb212Data = (TTree*)pb212.Get("data");

  TFile bi212("/Users/akindele1/Software/PMT/VetoStudy/ROOT_Files/twelve_std_glass/Bi212Bon.root");
  TTree *Bi212Data = (TTree*)bi212.Get("data");

  TFile tl208("/Users/akindele1/Software/PMT/VetoStudy/ROOT_Files/twelve_std_glass/TL208Bon.root");
  TTree *Tl208Data = (TTree*)tl208.Get("data"); //note branching ratio is 35.94%


  //U-238 Lines t1/2=4.468e9
  TFile bi210("/Users/akindele1/Software/PMT/VetoStudy/ROOT_Files/twelve_std_glass/Bi210Bon.root");
  TTree *Bi210Data = (TTree*)bi210.Get("data");

  TFile bi214("/Users/akindele1/Software/PMT/VetoStudy/ROOT_Files/twelve_std_glass/Bi214Bon.root");
  TTree *Bi214Data = (TTree*)bi214.Get("data");

  TFile pa234("/Users/akindele1/Software/PMT/VetoStudy/ROOT_Files/twelve_std_glass/Pa234Bon.root");
  TTree *Pa234Data = (TTree*)pa234.Get("data");

  TFile pb214("/Users/akindele1/Software/PMT/VetoStudy/ROOT_Files/twelve_std_glass/pb214Bon.root");
  TTree *Pb214Data = (TTree*)pb214.Get("data");//note branching ratio is 99.98%

  TFile tl210("/Users/akindele1/Software/PMT/VetoStudy/ROOT_Files/twelve_std_glass/Tl210Bon.root");
  TTree *Tl210Data = (TTree*)tl210.Get("data");//note branching ratio is 0.02%
  //======================================================
  //Making the Histograms of Interest
  //======================================================
  TH2D* inner_veto = new TH2D("In_Veto", "Inner vs Veto PMTs;Inner PMTs; Veto PMTs",201,-0.5,200.5,101,-0.5,100.5);
  TH2D* inner_n9   = new TH2D("In_n9", "Inner vs n9; Inner PMTs; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* veto_n9    = new TH2D("Ve_n9","Veto vs n9; Vetos; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* veto_v_veto= new TH2D("Vetotal","Veto Total vs n9; Vetos; n9",101,-0.5,100.5,101,-0.5,100.5);
  TH2D* n9_radius  = new TH2D("n9_radius","n9 versus Radius; Radius (cm); n9 ",200,-0.5,10000.5,11,-0.5,10.5);
  TH2D* n9_z       = new TH2D("n9_z","n9 versus Radius; Radius (cm); n9 ",200,-0.5,10000.5,11,-0.5,10.5);
  TH2D* r_z        = new TH2D("r_z","Spacial Event Location; Radius(cm)",200,-0.5,10000.5,200,-0.5,10000.5);
  TH1D* zm         = new TH1D("Z","Location of Events; Z(cm);Number of Events",200,-0.5,10000.5);
  TH1D* rm         = new TH1D("Radius","Location of Events; R(cm);Number of Events",200,-0.5,10000.5);

  // Setting up the Activtiy for the glass calculations
  double  t_u238  = 4.468e9;
  double  t_th232 = 1.405e10;
  double  t_k40   = 1.251e9;

  double abund_u238  = 0.992745;
  double abund_th232 = 1.0;
  double abund_k40   = 0.00117;

  double mass_pmt = 2000; //grams
  double k40Mass  = 40.0;
  double u238Mass = 238.0;
  double th232Mass= 232.0;

  double Na = 6.022e23;
  double time = 365.25*24.0*3600.0;
  double events = 1.0e5;
  double num_pmts = 197.0;

  //Determining the Actuvity of the Standard Glass
  double k40_ppm  = 260.0*120e-6;
  double u238_ppm = 0.341;
  double th232_ppm= 1.33;

  double A_u238 = 0.69315/(t_u238*time)*(Na*mass_pmt*u238_ppm*1.0e-6)/u238Mass*num_pmts;
  double A_th232= 0.69315/(t_th232*time)*(Na*mass_pmt*th232_ppm*1.0e-6)/th232Mass*num_pmts;
  double A_k40= 0.69315/(t_k40*time)*(Na*mass_pmt*k40_ppm*1.0e-6)/k40Mass*num_pmts;

  printf("STANDARD Glass: Activity of U238: %f, Activity of Th232: %f, Activity of K40: %f \n",A_u238,A_th232,A_k40);

  // ==========================================================
  //       Potassium 40 Radioactivtiy
  //===========================================================
  int inner_hit, veto_hit, veto_plus_dr_hit;
  double n9, x, y, z;

  K40Data->SetBranchAddress("inner_hit",&inner_hit);
  K40Data->SetBranchAddress("veto_hit",&veto_hit);
  K40Data->SetBranchAddress("n9",&n9);
  K40Data->SetBranchAddress("x",&x);
  K40Data->SetBranchAddress("y",&y);
  K40Data->SetBranchAddress("z",&z);
  K40Data->SetBranchAddress("veto_plus_dr_hit",&veto_plus_dr_hit);

  TH2D* k40_inner_veto = new TH2D("k40_In_Veto", "Inner vs Veto PMTs;Inner PMTs; Veto PMTs",201,-0.5,200.5,101,-0.5,100.5);
  TH2D* k40_inner_n9   = new TH2D("k40_In_n9", "Inner vs n9; Inner PMTs; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* k40_veto_n9    = new TH2D("k40_Ve_n9","Veto vs n9; Vetos; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* k40_veto_v_veto= new TH2D("k40_Vetotal","Veto vs Veto Total; Vetos; Veto Total",101,-0.5,100.5,101,-0.5,100.5);
  TH2D* k40_n9_radius  = new TH2D("k40_n9_radius","n9 versus Radius; Radius (cm); n9",200,-0.5,10000.5,11,-0.5,10.5);
  TH2D* k40_n9_z       = new TH2D("k40_n9_z","n9 versus Radius; Radius (cm); n9 ",200,-0.5,10000.5,11,-0.5,10.5);
  TH2D* k40_r_z        = new TH2D("k40_r_z","Spacial Event Location; Radius(cm)",200,-0.5,10000.5,200,-0.5,10000.5);
  TH1D* k40_zm         = new TH1D("k40_z","Location of Events; Z(cm);Number of Events",200,-0.5,10000.5);
  TH1D* k40_rm         = new TH1D("k40_Radius","Location of Events; R(cm);Number of Events",200,-0.5,10000.5);

  for (int i=0; i<K40Data->GetEntries(); i++){

    K40Data->GetEntry(i);

    k40_inner_veto->Fill(inner_hit,veto_hit);
    k40_inner_n9->Fill(inner_hit,n9);
    k40_veto_n9->Fill(veto_hit,n9);
    k40_n9_radius->Fill(sqrt(x*x+y*y),n9);
    k40_veto_v_veto->Fill(veto_hit, veto_plus_dr_hit);
    k40_n9_z->Fill(abs(z),n9);
    k40_r_z->Fill(sqrt(x*x+y*y),abs(z));
    k40_zm->Fill(abs(z));
    k40_rm->Fill(sqrt(x*x+y*y));

  }
  // ==========================================================
  //       AC 228 Radioactivtiy
  //===========================================================
  Ac228Data->SetBranchAddress("inner_hit",&inner_hit);
  Ac228Data->SetBranchAddress("veto_hit",&veto_hit);
  Ac228Data->SetBranchAddress("n9",&n9);
  Ac228Data->SetBranchAddress("x",&x);
  Ac228Data->SetBranchAddress("y",&y);
  Ac228Data->SetBranchAddress("z",&z);
  Ac228Data->SetBranchAddress("veto_plus_dr_hit",&veto_plus_dr_hit);

  TH2D* ac228_inner_veto = new TH2D("ac228_In_Veto", "Inner vs Veto PMTs;Inner PMTs; Veto PMTs",201,-0.5,200.5,101,-0.5,100.5);
  TH2D* ac228_inner_n9   = new TH2D("ac228_In_n9", "Inner vs n9; Inner PMTs; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* ac228_veto_n9    = new TH2D("ac228_Ve_n9","Veto vs n9; Vetos; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* ac228_veto_v_veto= new TH2D("ac228_Vetotal","Veto vs Veto Total; Vetos; Veto Total",101,-0.5,100.5,101,-0.5,100.5);
  TH2D* ac228_n9_radius  = new TH2D("ac228_n9_radius","n9 versus Radius; Radius (cm); n9",200,-0.5,10000.5,11,-0.5,10.5);
  TH1D* ac228_rm         = new TH1D("ac228_Radius","Location of Events; R(cm);Number of Events",200,-0.5,10000.5);
  TH2D* ac228_n9_z       = new TH2D("ac228_n9_z","n9 versus Radius; Radius (cm); n9 ",200,-0.5,10000.5,11,-0.5,10.5);
  TH2D* ac228_r_z        = new TH2D("ac228_r_z","Spacial Event Location; Radius(cm)",200,-0.5,10000.5,200,-0.5,10000.5);
  TH1D* ac228_zm         = new TH1D("ac228_z","Location of Events; Z(cm);Number of Events",200,-0.5,10000.5);

  for (int i=0; i<Ac228Data->GetEntries(); i++){

    Ac228Data->GetEntry(i);

    ac228_inner_veto->Fill(inner_hit,veto_hit);
    ac228_inner_n9->Fill(inner_hit,n9);
    ac228_veto_n9->Fill(veto_hit,n9);
    ac228_n9_radius->Fill(sqrt(x*x+y*y),n9);
    ac228_veto_v_veto->Fill(veto_hit, veto_plus_dr_hit);
    ac228_rm->Fill(sqrt(x*x+y*y));
    ac228_n9_z->Fill(abs(z),n9);
    ac228_r_z->Fill(sqrt(x*x+y*y),abs(z));
    ac228_zm->Fill(abs(z));

  }

  // ==========================================================
  //      Pb 212 Radioactivtiy
  //===========================================================
  Pb212Data->SetBranchAddress("inner_hit",&inner_hit);
  Pb212Data->SetBranchAddress("veto_hit",&veto_hit);
  Pb212Data->SetBranchAddress("n9",&n9);
  Pb212Data->SetBranchAddress("x",&x);
  Pb212Data->SetBranchAddress("y",&y);
  Pb212Data->SetBranchAddress("z",&z);
  Pb212Data->SetBranchAddress("veto_plus_dr_hit",&veto_plus_dr_hit);

  TH2D* pb212_inner_veto = new TH2D("pb212_In_Veto", "Inner vs Veto PMTs;Inner PMTs; Veto PMTs",201,-0.5,200.5,101,-0.5,100.5);
  TH2D* pb212_inner_n9   = new TH2D("pb212_In_n9", "Inner vs n9; Inner PMTs; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* pb212_veto_n9    = new TH2D("pb212_Ve_n9","Veto vs n9; Vetos; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* pb212_veto_v_veto= new TH2D("pb212_Vetotal","Veto vs Veto Total; Vetos; Veto Total",101,-0.5,100.5,101,-0.5,100.5);
  TH2D* pb212_n9_radius  = new TH2D("pb212_n9_radius","n9 versus Radius; Radius (cm); n9",200,-0.5,10000.5,11,-0.5,10.5);
  TH1D* pb212_rm         = new TH1D("pb212_Radius","Location of Events; R(cm);Number of Events",200,-0.5,10000.5);
  TH2D* pb212_n9_z       = new TH2D("pb212_n9_z","n9 versus Radius; Radius (cm); n9 ",200,-0.5,10000.5,11,-0.5,10.5);
  TH2D* pb212_r_z        = new TH2D("pb212_r_z","Spacial Event Location; Radius(cm)",200,-0.5,10000.5,200,-0.5,10000.5);
  TH1D* pb212_zm         = new TH1D("pb212_z","Location of Events; Z(cm);Number of Events",200,-0.5,10000.5);

  for (int i=0; i<Pb212Data->GetEntries(); i++){

    Pb212Data->GetEntry(i);

    pb212_inner_veto->Fill(inner_hit,veto_hit);
    pb212_inner_n9->Fill(inner_hit,n9);
    pb212_veto_n9->Fill(veto_hit,n9);
    pb212_n9_radius->Fill(sqrt(x*x+y*y),n9);
    pb212_veto_v_veto->Fill(veto_hit, veto_plus_dr_hit);
    pb212_rm->Fill(sqrt(x*x+y*y));
    pb212_n9_z->Fill(abs(z),n9);
    pb212_r_z->Fill(sqrt(x*x+y*y),abs(z));
    pb212_zm->Fill(abs(z));

  }

  // ==========================================================
  //      Bi 212 Radioactivtiy
  //===========================================================
  Bi212Data->SetBranchAddress("inner_hit",&inner_hit);
  Bi212Data->SetBranchAddress("veto_hit",&veto_hit);
  Bi212Data->SetBranchAddress("n9",&n9);
  Bi212Data->SetBranchAddress("x",&x);
  Bi212Data->SetBranchAddress("y",&y);
  Bi212Data->SetBranchAddress("veto_plus_dr_hit",&veto_plus_dr_hit);

  TH2D* bi212_inner_veto = new TH2D("bi212_In_Veto", "Inner vs Veto PMTs;Inner PMTs; Veto PMTs",201,-0.5,200.5,101,-0.5,100.5);
  TH2D* bi212_inner_n9   = new TH2D("bi212_In_n9", "Inner vs n9; Inner PMTs; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* bi212_veto_n9    = new TH2D("bi212_Ve_n9","Veto vs n9; Vetos; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* bi212_veto_v_veto= new TH2D("bi212_Vetotal","Veto vs Veto Total; Vetos; Veto Total",101,-0.5,100.5,101,-0.5,100.5);
  TH2D* bi212_n9_radius  = new TH2D("bi212_n9_radius","n9 versus Radius; Radius (cm); n9",200,-0.5,10000.5,11,-0.5,10.5);
  TH1D* bi212_rm         = new TH1D("bi212_Radius","Location of Events; R(cm);Number of Events",200,-0.5,10000.5);
  TH2D* bi212_n9_z       = new TH2D("bi212_n9_z","n9 versus Radius; Radius (cm); n9 ",200,-0.5,10000.5,11,-0.5,10.5);
  TH2D* bi212_r_z        = new TH2D("bi212_r_z","Spacial Event Location; Radius(cm)",200,-0.5,10000.5,200,-0.5,10000.5);
  TH1D* bi212_zm         = new TH1D("bi212_z","Location of Events; Z(cm);Number of Events",200,-0.5,10000.5);

  for (int i=0; i<Bi212Data->GetEntries(); i++){

    Bi212Data->GetEntry(i);

    bi212_inner_veto->Fill(inner_hit,veto_hit);
    bi212_inner_n9->Fill(inner_hit,n9);
    bi212_veto_n9->Fill(veto_hit,n9);
    bi212_n9_radius->Fill(sqrt(x*x+y*y),n9);
    bi212_veto_v_veto->Fill(veto_hit, veto_plus_dr_hit);
    bi212_rm->Fill(sqrt(x*x+y*y));
    bi212_n9_z->Fill(abs(z),n9);
    bi212_r_z->Fill(sqrt(x*x+y*y),abs(z));
    bi212_zm->Fill(abs(z));

  }

  // ==========================================================
  //      Tl 208 Radioactivtiy
  //===========================================================
  Tl208Data->SetBranchAddress("inner_hit",&inner_hit);
  Tl208Data->SetBranchAddress("veto_hit",&veto_hit);
  Tl208Data->SetBranchAddress("n9",&n9);
  Tl208Data->SetBranchAddress("x",&x);
  Tl208Data->SetBranchAddress("y",&y);
  Tl208Data->SetBranchAddress("z",&z);
  Tl208Data->SetBranchAddress("veto_plus_dr_hit",&veto_plus_dr_hit);

  TH2D* tl208_inner_veto = new TH2D("tl208_In_Veto", "Inner vs Veto PMTs;Inner PMTs; Veto PMTs",201,-0.5,200.5,101,-0.5,100.5);
  TH2D* tl208_inner_n9   = new TH2D("tl208_In_n9", "Inner vs n9; Inner PMTs; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* tl208_veto_n9    = new TH2D("tl208_Ve_n9","Veto vs n9; Vetos; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* tl208_veto_v_veto= new TH2D("tl208_Vetotal","Veto vs Veto Total; Vetos; Veto Total",101,-0.5,100.5,101,-0.5,100.5);
  TH2D* tl208_n9_radius  = new TH2D("tl208_n9_radius","n9 versus Radius; Radius (cm); n9",200,-0.5,10000.5,11,-0.5,10.5);
  TH1D* tl208_rm         = new TH1D("tl208_Radius","Location of Events; R(cm);Number of Events",200,-0.5,10000.5);
  TH2D* tl208_n9_z       = new TH2D("tl208_n9_z","n9 versus Radius; Radius (cm); n9 ",200,-0.5,10000.5,11,-0.5,10.5);
  TH2D* tl208_r_z        = new TH2D("tl208_r_z","Spacial Event Location; Radius(cm)",200,-0.5,10000.5,200,-0.5,10000.5);
  TH1D* tl208_zm         = new TH1D("tl208_z","Location of Events; Z(cm);Number of Events",200,-0.5,10000.5);

  for (int i=0; i<Tl208Data->GetEntries(); i++){

    Tl208Data->GetEntry(i);

    tl208_inner_veto->Fill(inner_hit,veto_hit);
    tl208_inner_n9->Fill(inner_hit,n9);
    tl208_veto_n9->Fill(veto_hit,n9);
    tl208_n9_radius->Fill(sqrt(x*x+y*y),n9);
    tl208_veto_v_veto->Fill(veto_hit, veto_plus_dr_hit);
    tl208_rm->Fill(sqrt(x*x+y*y));
    tl208_n9_z->Fill(abs(z),n9);
    tl208_r_z->Fill(sqrt(x*x+y*y),abs(z));
    tl208_zm->Fill(abs(z));
  }

  // ==========================================================
  //      Bi 210 Radioactivtiy
  //===========================================================
  Bi210Data->SetBranchAddress("inner_hit",&inner_hit);
  Bi210Data->SetBranchAddress("veto_hit",&veto_hit);
  Bi210Data->SetBranchAddress("n9",&n9);
  Bi210Data->SetBranchAddress("x",&x);
  Bi210Data->SetBranchAddress("y",&y);
  Bi210Data->SetBranchAddress("z",&z);
  Bi210Data->SetBranchAddress("veto_plus_dr_hit",&veto_plus_dr_hit);

  TH2D* bi210_inner_veto = new TH2D("bi210_In_Veto", "Inner vs Veto PMTs;Inner PMTs; Veto PMTs",201,-0.5,200.5,101,-0.5,100.5);
  TH2D* bi210_inner_n9   = new TH2D("bi210_In_n9", "Inner vs n9; Inner PMTs; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* bi210_veto_n9    = new TH2D("bi210_Ve_n9","Veto vs n9; Vetos; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* bi210_veto_v_veto= new TH2D("bi210_Vetotal","Veto vs Veto Total; Vetos; Veto Total",101,-0.5,100.5,101,-0.5,100.5);
  TH2D* bi210_n9_radius  = new TH2D("bi210_n9_radius","n9 versus Radius; Radius (cm); n9",200,-0.5,10000.5,11,-0.5,10.5);
  TH1D* bi210_rm         = new TH1D("bi210_Radius","Location of Events; R(cm);Number of Events",200,-0.5,10000.5);
  TH2D* bi210_n9_z       = new TH2D("bi210_n9_z","n9 versus Radius; Radius (cm); n9 ",200,-0.5,10000.5,11,-0.5,10.5);
  TH2D* bi210_r_z        = new TH2D("bi210_r_z","Spacial Event Location; Radius(cm)",200,-0.5,10000.5,200,-0.5,10000.5);
  TH1D* bi210_zm         = new TH1D("bi210_z","Location of Events; Z(cm);Number of Events",200,-0.5,10000.5);

  for (int i=0; i<Bi210Data->GetEntries(); i++){

    Bi210Data->GetEntry(i);

    bi210_inner_veto->Fill(inner_hit,veto_hit);
    bi210_inner_n9->Fill(inner_hit,n9);
    bi210_veto_n9->Fill(veto_hit,n9);
    bi210_n9_radius->Fill(sqrt(x*x+y*y),n9);
    bi210_veto_v_veto->Fill(veto_hit, veto_plus_dr_hit);
    bi210_rm->Fill(sqrt(x*x+y*y));
    bi210_n9_z->Fill(abs(z),n9);
    bi210_r_z->Fill(sqrt(x*x+y*y),abs(z));
    bi210_zm->Fill(abs(z));

  }

  // ==========================================================
  //      Bi 214 Radioactivtiy
  //===========================================================
  Bi214Data->SetBranchAddress("inner_hit",&inner_hit);
  Bi214Data->SetBranchAddress("veto_hit",&veto_hit);
  Bi214Data->SetBranchAddress("n9",&n9);
  Bi214Data->SetBranchAddress("x",&x);
  Bi214Data->SetBranchAddress("y",&y);
  Bi214Data->SetBranchAddress("z",&z);
  Bi214Data->SetBranchAddress("veto_plus_dr_hit",&veto_plus_dr_hit);

  TH2D* bi214_inner_veto = new TH2D("bi214_In_Veto", "Inner vs Veto PMTs;Inner PMTs; Veto PMTs",201,-0.5,200.5,101,-0.5,100.5);
  TH2D* bi214_inner_n9   = new TH2D("bi214_In_n9", "Inner vs n9; Inner PMTs; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* bi214_veto_n9    = new TH2D("bi214_Ve_n9","Veto vs n9; Vetos; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* bi214_veto_v_veto= new TH2D("bi214_Vetotal","Veto vs Veto Total; Vetos; Veto Total",101,-0.5,100.5,101,-0.5,100.5);
  TH2D* bi214_n9_radius  = new TH2D("bi214_n9_radius","n9 versus Radius; Radius (cm); n9",200,-0.5,10000.5,11,-0.5,10.5);
  TH1D* bi214_rm         = new TH1D("bi214_Radius","Location of Events; R(cm);Number of Events",200,-0.5,10000.5);
  TH2D* bi214_n9_z       = new TH2D("bi214_n9_z","n9 versus Radius; Radius (cm); n9 ",200,-0.5,10000.5,11,-0.5,10.5);
  TH2D* bi214_r_z        = new TH2D("bi214_r_z","Spacial Event Location; Radius(cm)",200,-0.5,10000.5,200,-0.5,10000.5);
  TH1D* bi214_zm         = new TH1D("bi214_z","Location of Events; Z(cm);Number of Events",200,-0.5,10000.5);

  for (int i=0; i<Bi214Data->GetEntries(); i++){

    Bi214Data->GetEntry(i);

    bi214_inner_veto->Fill(inner_hit,veto_hit);
    bi214_inner_n9->Fill(inner_hit,n9);
    bi214_veto_n9->Fill(veto_hit,n9);
    bi214_n9_radius->Fill(sqrt(x*x+y*y),n9);
    bi214_veto_v_veto->Fill(veto_hit, veto_plus_dr_hit);
    bi214_rm->Fill(sqrt(x*x+y*y));
    bi214_n9_z->Fill(abs(z),n9);
    bi214_r_z->Fill(sqrt(x*x+y*y),abs(z));
    bi214_zm->Fill(abs(z));
  }

  // ==========================================================
  //      Pa 234 Radioactivtiy
  //===========================================================
  Pa234Data->SetBranchAddress("inner_hit",&inner_hit);
  Pa234Data->SetBranchAddress("veto_hit",&veto_hit);
  Pa234Data->SetBranchAddress("n9",&n9);
  Pa234Data->SetBranchAddress("x",&x);
  Pa234Data->SetBranchAddress("y",&y);
  Pa234Data->SetBranchAddress("z",&z);
  Pa234Data->SetBranchAddress("veto_plus_dr_hit",&veto_plus_dr_hit);

  TH2D* pa234_inner_veto = new TH2D("pa234_In_Veto", "Inner vs Veto PMTs;Inner PMTs; Veto PMTs",201,-0.5,200.5,101,-0.5,100.5);
  TH2D* pa234_inner_n9   = new TH2D("pa234_In_n9", "Inner vs n9; Inner PMTs; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* pa234_veto_n9    = new TH2D("pa234_Ve_n9","Veto vs n9; Vetos; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* pa234_veto_v_veto= new TH2D("pa234_Vetotal","Veto vs Veto Total; Vetos; Veto Total",101,-0.5,100.5,101,-0.5,100.5);
  TH2D* pa234_n9_radius  = new TH2D("pa234_n9_radius","n9 versus Radius; Radius (cm); n9",200,-0.5,10000.5,11,-0.5,10.5);
  TH1D* pa234_rm         = new TH1D("pa234_Radius","Location of Events; R(cm);Number of Events",200,-0.5,10000.5);
  TH2D* pa234_n9_z       = new TH2D("pa234_n9_z","n9 versus Radius; Radius (cm); n9 ",200,-0.5,10000.5,11,-0.5,10.5);
  TH2D* pa234_r_z        = new TH2D("pa234_r_z","Spacial Event Location; Radius(cm)",200,-0.5,10000.5,200,-0.5,10000.5);
  TH1D* pa234_zm         = new TH1D("pa234_z","Location of Events; Z(cm);Number of Events",200,-0.5,10000.5);

  for (int i=0; i<Pa234Data->GetEntries(); i++){

    Pa234Data->GetEntry(i);

    pa234_inner_veto->Fill(inner_hit,veto_hit);
    pa234_inner_n9->Fill(inner_hit,n9);
    pa234_veto_n9->Fill(veto_hit,n9);
    pa234_n9_radius->Fill(sqrt(x*x+y*y),n9);
    pa234_veto_v_veto->Fill(veto_hit, veto_plus_dr_hit);
    pa234_rm->Fill(sqrt(x*x+y*y));
    pa234_n9_z->Fill(abs(z),n9);
    pa234_r_z->Fill(sqrt(x*x+y*y),abs(z));
    pa234_zm->Fill(abs(z));

  }

  // ==========================================================
  //      Pb 214 Radioactivtiy
  //===========================================================
  Pb214Data->SetBranchAddress("inner_hit",&inner_hit);
  Pb214Data->SetBranchAddress("veto_hit",&veto_hit);
  Pb214Data->SetBranchAddress("n9",&n9);
  Pb214Data->SetBranchAddress("x",&x);
  Pb214Data->SetBranchAddress("y",&y);
  Pb214Data->SetBranchAddress("z",&z);
  Pb214Data->SetBranchAddress("veto_plus_dr_hit",&veto_plus_dr_hit);

  TH2D* pb214_inner_veto = new TH2D("pb214_In_Veto", "Inner vs Veto PMTs;Inner PMTs; Veto PMTs",201,-0.5,200.5,101,-0.5,100.5);
  TH2D* pb214_inner_n9   = new TH2D("pb214_In_n9", "Inner vs n9; Inner PMTs; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* pb214_veto_n9    = new TH2D("pb214_Ve_n9","Veto vs n9; Vetos; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* pb214_veto_v_veto= new TH2D("pb214_Vetotal","Veto vs Veto Total; Vetos; Veto Total",101,-0.5,100.5,101,-0.5,100.5);
  TH2D* pb214_n9_radius  = new TH2D("pb214_n9_radius","n9 versus Radius; Radius (cm); n9",200,-0.5,10000.5,11,-0.5,10.5);
  TH1D* pb214_rm         = new TH1D("pb214_Radius","Location of Events; R(cm);Number of Events",200,-0.5,10000.5);
  TH2D* pb214_n9_z       = new TH2D("pb214_n9_z","n9 versus Radius; Radius (cm); n9 ",200,-0.5,10000.5,11,-0.5,10.5);
  TH2D* pb214_r_z        = new TH2D("pb214_r_z","Spacial Event Location; Radius(cm)",200,-0.5,10000.5,200,-0.5,10000.5);
  TH1D* pb214_zm         = new TH1D("pb214_z","Location of Events; Z(cm);Number of Events",200,-0.5,10000.5);

  for (int i=0; i<Pb214Data->GetEntries(); i++){

    Pb214Data->GetEntry(i);

    pb214_inner_veto->Fill(inner_hit,veto_hit);
    pb214_inner_n9->Fill(inner_hit,n9);
    pb214_veto_n9->Fill(veto_hit,n9);
    pb214_n9_radius->Fill(sqrt(x*x+y*y),n9);
    pb214_veto_v_veto->Fill(veto_hit, veto_plus_dr_hit);
    pb214_rm->Fill(sqrt(x*x+y*y));
    pb214_n9_z->Fill(abs(z),n9);
    pb214_r_z->Fill(sqrt(x*x+y*y),abs(z));
    pb214_zm->Fill(abs(z));
  }

  // ==========================================================
  //      Tl 210 Radioactivtiy
  //===========================================================
  Tl210Data->SetBranchAddress("inner_hit",&inner_hit);
  Tl210Data->SetBranchAddress("veto_hit",&veto_hit);
  Tl210Data->SetBranchAddress("n9",&n9);
  Tl210Data->SetBranchAddress("x",&x);
  Tl210Data->SetBranchAddress("y",&y);
  Tl210Data->SetBranchAddress("z",&z);
  Tl210Data->SetBranchAddress("veto_plus_dr_hit",&veto_plus_dr_hit);

  TH2D* tl210_inner_veto = new TH2D("tl210_In_Veto", "Inner vs Veto PMTs;Inner PMTs; Veto PMTs",201,-0.5,200.5,101,-0.5,100.5);
  TH2D* tl210_inner_n9   = new TH2D("tl210_In_n9", "Inner vs n9; Inner PMTs; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* tl210_veto_n9    = new TH2D("tl210_Ve_n9","Veto vs n9; Vetos; n9",201,-0.5,200.5,11,-0.5,10.5);
  TH2D* tl210_veto_v_veto= new TH2D("tl210_Vetotal","Veto vs Veto Total; Vetos; Veto Total",101,-0.5,100.5,101,-0.5,100.5);
  TH2D* tl210_n9_radius  = new TH2D("tl210_n9_radius","n9 versus Radius; Radius (cm); n9",200,-0.5,10000.5,11,-0.5,10.5);
  TH1D* tl210_rm         = new TH1D("tl210_Radius","Location of Events; R(cm);Number of Events",200,-0.5,10000.5);
  TH2D* tl210_n9_z       = new TH2D("tl210_n9_z","n9 versus Radius; Radius (cm); n9 ",200,-0.5,10000.5,11,-0.5,10.5);
  TH2D* tl210_r_z        = new TH2D("tl210_r_z","Spacial Event Location; Radius(cm)",200,-0.5,10000.5,200,-0.5,10000.5);
  TH1D* tl210_zm         = new TH1D("tl210_z","Location of Events; Z(cm);Number of Events",200,-0.5,10000.5);

  for (int i=0; i<Tl210Data->GetEntries(); i++){

    Tl210Data->GetEntry(i);

    tl210_inner_veto->Fill(inner_hit,veto_hit);
    tl210_inner_n9->Fill(inner_hit,n9);
    tl210_veto_n9->Fill(veto_hit,n9);
    tl210_n9_radius->Fill(sqrt(x*x+y*y),n9);
    tl210_veto_v_veto->Fill(veto_hit, veto_plus_dr_hit);
    tl210_rm->Fill(sqrt(x*x+y*y));
    tl210_n9_z->Fill(abs(z),n9);
    tl210_r_z->Fill(sqrt(x*x+y*y),abs(z));
    tl210_zm->Fill(abs(z));
  }

  //============================================================
  //     Scale by the Activity and Add
  //============================================================
  k40_inner_veto->Scale(A_k40/events);
  ac228_inner_veto->Scale(A_th232/events);pb212_inner_veto->Scale(A_th232/events);bi212_inner_veto->Scale(A_th232/events);tl208_inner_veto->Scale(A_th232*.3594/events);
  bi210_inner_veto->Scale(A_u238/events);bi214_inner_veto->Scale(A_u238/events);pa234_inner_veto->Scale(A_u238/events);pb214_inner_veto->Scale(A_u238*.9998/events);tl210_inner_veto->Scale(A_u238*.0002/events);

  inner_veto->Add(k40_inner_veto,inner_veto);inner_veto->Add(ac228_inner_veto,inner_veto);
  inner_veto->Add(pb212_inner_veto,inner_veto);inner_veto->Add(bi212_inner_veto,inner_veto);
  inner_veto->Add(tl208_inner_veto,inner_veto);inner_veto->Add(bi210_inner_veto,inner_veto);
  inner_veto->Add(bi214_inner_veto,inner_veto);inner_veto->Add(pa234_inner_veto,inner_veto);
  inner_veto->Add(pb214_inner_veto,inner_veto);inner_veto->Add(tl210_inner_veto,inner_veto);
  //==============================================================
  k40_inner_n9->Scale(A_k40/events);
  ac228_inner_n9->Scale(A_th232/events);pb212_inner_n9->Scale(A_th232/events);bi212_inner_n9->Scale(A_th232/events);tl208_inner_n9->Scale(A_th232*.3594/events);
  bi210_inner_n9->Scale(A_u238/events);bi214_inner_n9->Scale(A_u238/events);pa234_inner_n9->Scale(A_u238/events);pb214_inner_n9->Scale(A_u238*.9998/events);tl210_inner_n9->Scale(A_u238*.0002/events);

  inner_n9->Add(k40_inner_n9,inner_n9);inner_n9->Add(ac228_inner_n9,inner_n9);
  inner_n9->Add(pb212_inner_n9,inner_n9);inner_n9->Add(bi212_inner_n9,inner_n9);
  inner_n9->Add(tl208_inner_n9,inner_n9);inner_n9->Add(bi210_inner_n9,inner_n9);
  inner_n9->Add(bi214_inner_n9,inner_n9);inner_n9->Add(pa234_inner_n9,inner_n9);
  inner_n9->Add(pb214_inner_n9,inner_n9);inner_n9->Add(tl210_inner_n9,inner_n9);
  //==============================================================
  k40_veto_n9->Scale(A_k40/events);
  ac228_veto_n9->Scale(A_th232/events);pb212_veto_n9->Scale(A_th232/events);bi212_veto_n9->Scale(A_th232/events);tl208_veto_n9->Scale(A_th232*.3594/events);
  bi210_veto_n9->Scale(A_u238/events);bi214_veto_n9->Scale(A_u238/events);pa234_veto_n9->Scale(A_u238/events);pb214_veto_n9->Scale(A_u238*.9998/events);tl210_veto_n9->Scale(A_u238*.0002/events);

  veto_n9->Add(k40_veto_n9,veto_n9);veto_n9->Add(ac228_veto_n9,veto_n9);
  veto_n9->Add(pb212_veto_n9,veto_n9);veto_n9->Add(bi212_veto_n9,veto_n9);
  veto_n9->Add(tl208_veto_n9,veto_n9);veto_n9->Add(bi210_veto_n9,veto_n9);
  veto_n9->Add(bi214_veto_n9,veto_n9);veto_n9->Add(pa234_veto_n9,veto_n9);
  veto_n9->Add(pb214_veto_n9,veto_n9);veto_n9->Add(tl210_veto_n9,veto_n9);

  //==============================================================
  k40_n9_radius->Scale(A_k40/events);
  ac228_n9_radius->Scale(A_th232/events);pb212_n9_radius->Scale(A_th232/events);bi212_n9_radius->Scale(A_th232/events);tl208_n9_radius->Scale(A_th232*.3594/events);
  bi210_n9_radius->Scale(A_u238/events);bi214_n9_radius->Scale(A_u238/events);pa234_n9_radius->Scale(A_u238/events);pb214_n9_radius->Scale(A_u238*.9998/events);tl210_n9_radius->Scale(A_u238*.0002/events);

  n9_radius->Add(k40_n9_radius,n9_radius);n9_radius->Add(ac228_n9_radius,n9_radius);
  n9_radius->Add(pb212_n9_radius,n9_radius);n9_radius->Add(bi212_n9_radius,n9_radius);
  n9_radius->Add(tl208_n9_radius,n9_radius);n9_radius->Add(bi210_n9_radius,n9_radius);
  n9_radius->Add(bi214_n9_radius,n9_radius);n9_radius->Add(pa234_n9_radius,n9_radius);
  n9_radius->Add(pb214_n9_radius,n9_radius);n9_radius->Add(tl210_n9_radius,n9_radius);

  //==============================================================
  k40_veto_v_veto->Scale(A_k40/events);
  ac228_veto_v_veto->Scale(A_th232/events);pb212_veto_v_veto->Scale(A_th232/events);bi212_veto_v_veto->Scale(A_th232/events);tl208_veto_v_veto->Scale(A_th232*.3594/events);
  bi210_veto_v_veto->Scale(A_u238/events);bi214_veto_v_veto->Scale(A_u238/events);pa234_veto_v_veto->Scale(A_u238/events);pb214_veto_v_veto->Scale(A_u238*.9998/events);tl210_veto_v_veto->Scale(A_u238*.0002/events);

  veto_v_veto->Add(k40_veto_v_veto,veto_v_veto);veto_v_veto->Add(ac228_veto_v_veto,veto_v_veto);
  veto_v_veto->Add(pb212_veto_v_veto,veto_v_veto);veto_v_veto->Add(bi212_veto_v_veto,veto_v_veto);
  veto_v_veto->Add(tl208_veto_v_veto,veto_v_veto);veto_v_veto->Add(bi210_veto_v_veto,veto_v_veto);
  veto_v_veto->Add(bi214_veto_v_veto,veto_v_veto);veto_v_veto->Add(pa234_veto_v_veto,veto_v_veto);
  veto_v_veto->Add(pb214_veto_v_veto,veto_v_veto);veto_v_veto->Add(tl210_veto_v_veto,veto_v_veto);

  //==============================================================
  k40_rm->Scale(A_k40/events);
  ac228_rm->Scale(A_th232/events);pb212_rm->Scale(A_th232/events);bi212_rm->Scale(A_th232/events);tl208_rm->Scale(A_th232*.3594/events);
  bi210_rm->Scale(A_u238/events);bi214_rm->Scale(A_u238/events);pa234_rm->Scale(A_u238/events);pb214_rm->Scale(A_u238*.9998/events);tl210_rm->Scale(A_u238*.0002/events);

  rm->Add(k40_rm,rm);rm->Add(ac228_rm,rm);
  rm->Add(pb212_rm,rm);rm->Add(bi212_rm,rm);
  rm->Add(tl208_rm,rm);rm->Add(bi210_rm,rm);
  rm->Add(bi214_rm,rm);rm->Add(pa234_rm,rm);
  rm->Add(pb214_rm,rm);rm->Add(tl210_rm,rm);
  //===============================================================
  k40_zm->Scale(A_k40/events);
  ac228_zm->Scale(A_th232/events);pb212_zm->Scale(A_th232/events);bi212_zm->Scale(A_th232/events);tl208_zm->Scale(A_th232*.3594/events);
  bi210_zm->Scale(A_u238/events);bi214_zm->Scale(A_u238/events);pa234_zm->Scale(A_u238/events);pb214_zm->Scale(A_u238*.9998/events);tl210_zm->Scale(A_u238*.0002/events);

  zm->Add(k40_zm,zm);zm->Add(ac228_zm,zm);
  zm->Add(pb212_zm,zm);zm->Add(bi212_zm,zm);
  zm->Add(tl208_zm,zm);zm->Add(bi210_zm,zm);
  zm->Add(bi214_zm,zm);zm->Add(pa234_zm,zm);
  zm->Add(pb214_zm,zm);zm->Add(tl210_zm,zm);
  //==============================================================
  k40_n9_z->Scale(A_k40/events);
  ac228_n9_z->Scale(A_th232/events);pb212_n9_z->Scale(A_th232/events);bi212_n9_z->Scale(A_th232/events);tl208_n9_z->Scale(A_th232*.3594/events);
  bi210_n9_z->Scale(A_u238/events);bi214_n9_z->Scale(A_u238/events);pa234_n9_z->Scale(A_u238/events);pb214_n9_z->Scale(A_u238*.9998/events);tl210_n9_z->Scale(A_u238*.0002/events);

  n9_z->Add(k40_n9_z,n9_z);n9_z->Add(ac228_n9_z,n9_z);
  n9_z->Add(pb212_n9_z,n9_z);n9_z->Add(bi212_n9_z,n9_z);
  n9_z->Add(tl208_n9_z,n9_z);n9_z->Add(bi210_n9_z,n9_z);
  n9_z->Add(bi214_n9_z,n9_z);n9_z->Add(pa234_n9_z,n9_z);
  n9_z->Add(pb214_n9_z,n9_z);n9_z->Add(tl210_n9_z,n9_z);
  //==============================================================
  k40_r_z->Scale(A_k40/events);
  ac228_r_z->Scale(A_th232/events);pb212_r_z->Scale(A_th232/events);bi212_r_z->Scale(A_th232/events);tl208_r_z->Scale(A_th232*.3594/events);
  bi210_r_z->Scale(A_u238/events);bi214_r_z->Scale(A_u238/events);pa234_r_z->Scale(A_u238/events);pb214_r_z->Scale(A_u238*.9998/events);tl210_r_z->Scale(A_u238*.0002/events);

  r_z->Add(k40_r_z,r_z);r_z->Add(ac228_r_z,r_z);
  r_z->Add(pb212_r_z,r_z);r_z->Add(bi212_r_z,r_z);
  r_z->Add(tl208_r_z,r_z);r_z->Add(bi210_r_z,r_z);
  r_z->Add(bi214_r_z,r_z);r_z->Add(pa234_r_z,r_z);
  r_z->Add(pb214_r_z,r_z);r_z->Add(tl210_r_z,r_z);

  TCanvas *tc0= new TCanvas("Canvas0","ROOT Canvas",1);
  tc0->Connect("TCanvas","Closed()","TApplication",gApplication, "Terminate()");
  tc0->SetGrid();

  //tl208_inner_veto->Draw("Same");
  inner_veto->Draw("Same");

  TCanvas *tc1= new TCanvas("Canvas1","ROOT Canvas",1);
  tc1->Connect("TCanvas","Closed()","TApplication",gApplication, "Terminate()");
  tc1->SetGrid();

  inner_n9->Draw("Same");

  TCanvas *tc2= new TCanvas("Canvas2","ROOT Canvas",1);
  tc2->Connect("TCanvas","Closed()","TApplication",gApplication, "Terminate()");
  tc2->SetGrid();

  veto_n9->Draw("Same");

  TCanvas *tc3= new TCanvas("Canvas3","ROOT Canvas",1);
  tc3->Connect("TCanvas","Closed()","TApplication",gApplication, "Terminate()");
  tc3->SetGrid();

  n9_radius->Draw("Same");

  TCanvas *tc4= new TCanvas("Canvas4","ROOT Canvas",1);
  tc4->Connect("TCanvas","Closed()","TApplication",gApplication, "Terminate()");
  tc4->SetGrid();

  veto_v_veto->Draw("Same");

  TCanvas *tc5= new TCanvas("Canvas5","ROOT Canvas",1);
  tc5->Connect("TCanvas","Closed()","TApplication",gApplication, "Terminate()");
  tc5->SetGrid();

  rm->SetLineColor(kBlack);k40_rm->SetLineColor(kGreen);ac228_rm->SetLineColor(kBlue);pb212_rm->SetLineColor(kCyan);bi212_rm->SetLineColor(kOrange);tl208_rm->SetLineColor(kMagenta);bi210_rm->SetLineColor(kRed);bi214_rm->SetLineColor(kTeal);pa234_rm->SetLineColor(kViolet);pb214_rm->SetLineColor(kAzure-8);tl210_rm->SetLineColor(kYellow);
  rm->SetLineWidth(2);k40_rm->SetLineWidth(2);ac228_rm->SetLineWidth(2);pb212_rm->SetLineWidth(2);bi212_rm->SetLineWidth(2);tl208_rm->SetLineWidth(2);bi210_rm->SetLineWidth(2);bi214_rm->SetLineWidth(2);pa234_rm->SetLineWidth(2);pb214_rm->SetLineWidth(2);tl210_rm->SetLineWidth(2);
  rm->SetLineStyle(21);k40_rm->SetLineStyle(1);ac228_rm->SetLineStyle(1);pb212_rm->SetLineStyle(1);bi212_rm->SetLineStyle(1);tl208_rm->SetLineStyle(1);bi210_rm->SetLineStyle(1);bi214_rm->SetLineStyle(1);pa234_rm->SetLineStyle(1);pb214_rm->SetLineStyle(1);tl210_rm->SetLineStyle(1);

  rm->Draw();k40_rm->Draw("Same");ac228_rm->Draw("Same");pb212_rm->Draw("Same");bi212_rm->Draw("Same");tl208_rm->Draw("Same");bi210_rm->Draw("Same");bi214_rm->Draw("Same");pa234_rm->Draw("Same");pb214_rm->Draw("Same");tl210_rm->Draw("Same");

  TLegend *l5 = new TLegend (0.1,0.7,0.48,0.9);
  l5->AddEntry(rm,"Total","l");l5->AddEntry(k40_rm,"K-40","l");l5->AddEntry(ac228_rm,"Ac-228","l");l5->AddEntry(pb212_rm,"Pb-212","l");l5->AddEntry(bi212_rm,"Bi-212","l");l5->AddEntry(tl208_rm,"Tl-208","l");l5->AddEntry(bi210_rm,"Bi-210","l");l5->AddEntry(bi214_rm,"Bi-214","l");l5->AddEntry(pa234_rm,"Pa-234","l");l5->AddEntry(pb214_rm,"Pb-214","l");l5->AddEntry(tl210_rm,"Tl-210","l");
  l5->Draw();
  // Update canvas.
	tc5->Update();
	tc5->Paint();
	tc5->Draw();
  tc5->Modified();

  ta->Run();
  return 0;
}
