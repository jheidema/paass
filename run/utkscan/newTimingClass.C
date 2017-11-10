#define newTimingClass_cxx
#include "newTimingClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void newTimingClass::Loop(Long64_t nentries,const Char_t *filename)
{
//   In a ROOT session, you can do:
//      root> .L newTimingClass.C
//      root> newTimingClass t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   TFile *outputFile;
   if(!filename)
     outputFile=new TFile("Output.root","RECREATE");
   else
     outputFile=new TFile(filename,"RECREATE");
   TTree *outTree=new TTree("T","Output Tree");
   outTree->Branch("ROOT_phase[4]",&ROOT_phase,"phase[4]/D");
   outTree->Branch("ROOT_beta[4]",&ROOT_beta,"beta[4]/D");
   outTree->Branch("ROOT_gamma[4]",&ROOT_gamma,"gamma[4]/D");
   outTree->Branch("ROOT_ChiSq[4]",&ROOT_chisq,"ChiSq[4]/D");
   outTree->Branch("ROOT_Status[4]",&fStatus,"Stat[4]/I");
   outTree->Branch("ROOT_QDC[4]",&ROOT_qdc,"QDC[4]/D");
   outTree->Branch("ROOT_time[4]",&ROOT_time,"time[4]/D");
   outTree->Branch("ROOT_max[4]",&ROOT_max,"Max[4]/I");
   
   if(nentries==-1)
     nentries = fChain->GetEntriesFast();
   
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

     if(jentry%10000==0)
	cout<<"."<<flush;
     //Fit(jentry,kTRUE);
     Fit(jentry,kFALSE); //No fix beta and gamma
     //cout<<jentry<<"Filling... "<<endl;
     outTree->Fill();
     
     // if (Cut(ientry) < 0) continue;
   }
   outTree->Write();
   cout<<endl;
   
}
