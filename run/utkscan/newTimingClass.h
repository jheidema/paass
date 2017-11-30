//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 25 10:46:56 2017 by ROOT version 6.06/08
// from TTree timing/
// found on file: SmallPlastic_4A7C_ArrayJ_20170922_003_newProc.root
//////////////////////////////////////////////////////////

#ifndef newTimingClass_h
#define newTimingClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TTree.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <iostream>


// Header file for the classes stored in the TTree if any.
#include "vector"
#include <vector>
#include <algorithm>
#include <utility>
#include <iterator> 

#include "SiPMTimingFunction.cpp"
using namespace std;
class newTimingClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        start1_qdc;
   Double_t        start1_time;
   Double_t        start1_snr;
   Double_t        start1_wtime;
   Double_t        start1_phase;
   Double_t        start1_abase;
   Double_t        start1_sbase;
   UChar_t         start1_id;
   Double_t        stop1_qdc;
   Double_t        stop1_time;
   Double_t        stop1_snr;
   Double_t        stop1_wtime;
   Double_t        stop1_phase;
   Double_t        stop1_abase;
   Double_t        stop1_sbase;
   UChar_t         stop1_id;
   Double_t        start2_qdc;
   Double_t        start2_time;
   Double_t        start2_snr;
   Double_t        start2_wtime;
   Double_t        start2_phase;
   Double_t        start2_abase;
   Double_t        start2_sbase;
   UChar_t         start2_id;
   Double_t        stop2_qdc;
   Double_t        stop2_time;
   Double_t        stop2_snr;
   Double_t        stop2_wtime;
   Double_t        stop2_phase;
   Double_t        stop2_abase;
   Double_t        stop2_sbase;
   UChar_t         stop2_id;
   vector<unsigned int> *trace_start1;
   vector<unsigned int> *trace_start2;
   vector<unsigned int> *trace_stop1;
   vector<unsigned int> *trace_stop2;
   Double_t        StartTimeStamp[2];
   Double_t        StopTimeStamp[2];
   Int_t           StartMaximum[2];
   Int_t           StopTimeMaximum[2];
   Double_t        StartChiSq;
   Double_t        StopChiSq;

   //My Stuff
   Double_t beta[4];
   Double_t gamma[4];
   
   Double_t ROOT_phase[4];
   Double_t ROOT_time[4];
   Double_t ROOT_beta[4];
   Double_t ROOT_gamma[4];
   Double_t ROOT_chisq[4];
   Double_t ROOT_qdc[4];
   Double_t ROOT_max[4];

//   vector<unsigned int> *normtrace1;
//   vector<unsigned int> *normtrace2;
     
   TGraph *fTraces[4];
   TF1 *fFits[4];
   TGraph *fResiduals[4];
   TGraph *fXCorrelation[2];
   Int_t fStatus[4]={-1,-1,-1,-1};

   const Int_t factor=1;
   

   // List of branches
   TBranch        *b_start1;   //!
   TBranch        *b_stop1;   //!
   TBranch        *b_start2;   //!
   TBranch        *b_stop2;   //!
   TBranch        *b_trace_start1;   //!
   TBranch        *b_trace_start2;   //!
   TBranch        *b_trace_stop1;   //!
   TBranch        *b_trace_stop2;   //!
   TBranch        *b_StartTimestamp;   //!
   TBranch        *b_StopTimestamp;   //!
   TBranch        *b_StartMax;   //!
   TBranch        *b_StopMax;   //!
   TBranch        *b_StartChi;   //!
   TBranch        *b_StopChi;   //!

   newTimingClass(TTree *tree=0);
   virtual ~newTimingClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void Loop(Long64_t nentries =-1,const Char_t *filename=NULL);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1,Bool_t fix=kFALSE);

   //My methods
   virtual void     Plot(Long64_t entry = -1,Bool_t draw=kFALSE);
   virtual void     Fit(Long64_t entry = -1, Bool_t fix=kFALSE);
   virtual void     Residuals(Long64_t entry = -1, Bool_t fix=kFALSE);
   virtual void     Trap_filter(Long64_t entry = -1,UInt_t length=4,UInt_t gap=0);

   void SetBeta(Int_t n,Double_t val){beta[n]=val;}
   void SetGamma(Int_t n,Double_t val){gamma[n]=val;}

   Int_t CrossCorrelation(vector<unsigned int> *trace1,vector<unsigned int> *trace2, TGraph *fGraph);
//   void Normalize(vector<unsigned int> *trace1,vector<unsigned int> *trace2,vector<unsigned int> *normtrace1,vector<unsigned int> *normtrace2);
};

#endif

#ifdef newTimingClass_cxx
newTimingClass::newTimingClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("SmallPlastic_4A7C_ArrayJ_20170922_003_newProc.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("SmallPlastic_4A7C_ArrayJ_20170922_003_newProc.root");
      }
      f->GetObject("timing",tree);

   }
   Init(tree);

   Char_t fname[256];
  SiPMTimingFunction sipmfunc;

   for(Int_t j=0;j<4;j++){
   beta[j]=0.08;
   gamma[j]=0.26;
   sprintf(fname,"f_%d",j);
   fTraces[j]=new TGraph();
   fResiduals[j]=new TGraph();
   fFits[j]=new TF1(fname,sipmfunc,0,1,6);
   if(j<2)
   fXCorrelation[j]=new TGraph();
   }

}

newTimingClass::~newTimingClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
  for(Int_t j=0;j<4;j++){
    delete fResiduals[j];
    delete fTraces[j];
    delete fFits[j];
   if(j<2)
   delete fXCorrelation[j];

  }
}

Int_t newTimingClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t newTimingClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void newTimingClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trace_start1 = 0;
   trace_start2 = 0;
   trace_stop1 = 0;
   trace_stop2 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("start1", &start1_qdc, &b_start1);
   fChain->SetBranchAddress("stop1", &stop1_qdc, &b_stop1);
   fChain->SetBranchAddress("start2", &start2_qdc, &b_start2);
   fChain->SetBranchAddress("stop2", &stop2_qdc, &b_stop2);
   fChain->SetBranchAddress("trace_start1", &trace_start1, &b_trace_start1);
   fChain->SetBranchAddress("trace_start2", &trace_start2, &b_trace_start2);
   fChain->SetBranchAddress("trace_stop1", &trace_stop1, &b_trace_stop1);
   fChain->SetBranchAddress("trace_stop2", &trace_stop2, &b_trace_stop2);
   fChain->SetBranchAddress("StartTimeStamp[2]", StartTimeStamp, &b_StartTimestamp);
   fChain->SetBranchAddress("StopTimeStamp[2]", StopTimeStamp, &b_StopTimestamp);
   fChain->SetBranchAddress("StartMaximum[2]", StartMaximum, &b_StartMax);
   fChain->SetBranchAddress("StopTimeMaximum[2]", StopTimeMaximum, &b_StopMax);
   fChain->SetBranchAddress("StartChiSq", &StartChiSq, &b_StartChi);
   fChain->SetBranchAddress("StopChiSq", &StopChiSq, &b_StopChi);
   Notify();
}

Bool_t newTimingClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void newTimingClass::Show(Long64_t entry,Bool_t fix)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
   Fit(entry,fix);

   TCanvas *c;
   if(!gPad){
     c=new TCanvas();
     c->Divide(2,2);
   }
   else{
     c=gPad->GetCanvas();
     c->Clear();
     c->Divide(2,2);
   }   
   for(Int_t i=0;i<4;i++){
     c->cd(i+1);   
     if(fTraces[i]->GetN()>0){
       fTraces[i]->Draw("AL*");
       cout<<"ROOT_Status["<<i<<"] = "<<fStatus[i]<<endl;
       cout<<"ROOT_phase["<<i<<"] = "<<ROOT_phase[i]<<endl;
       cout<<"ROOT_beta["<<i<<"] = "<<ROOT_beta[i]<<endl;
       cout<<"ROOT_gamma["<<i<<"] = "<<ROOT_gamma[i]<<endl;
       cout<<"ROOT_chisq["<<i<<"] = "<<ROOT_chisq[i]<<endl;
       cout<<"ROOT_time["<<i<<"] = "<<ROOT_time[i]<<endl;
       cout<<"ROOT_max["<<i<<"] = "<<ROOT_max[i]<<endl;
       cout<<"----------------------------------------"<<endl;
       fFits[i]->Draw("same");
     }
     else
       continue;
   }


}
Int_t newTimingClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


void newTimingClass::Plot(Long64_t entry,Bool_t draw){

  GetEntry(entry);


  const UInt_t size0=trace_start1->size();
  const UInt_t size1=trace_start2->size();
  const UInt_t size2=trace_stop1->size();
  const UInt_t size3=trace_stop2->size();

  Int_t phase[2]={0,0};

  if(size0>0&&size1>0)
    phase[0]=CrossCorrelation(trace_start1,trace_start2,fXCorrelation[0]);
  if(size2>0&&size3>0)
    phase[1]=CrossCorrelation(trace_stop1,trace_stop2,fXCorrelation[1]);
  cout << "Phase: " << phase[0] << " " << phase[1] << endl;
  /* cout<<"Size "<<size0<<endl;  */
  /* cout<<"Size "<<size1<<endl;  */
  /* cout<<"Size "<<size2<<endl;  */
  /* cout<<"Size "<<size3<<endl;  */
  
  for(Int_t j=0;j<4;j++){
    fTraces[j]->Set(0); 
  }

  if(size0!=0){
    for(UInt_t j=0;j<size0;j++){
      fTraces[0]->SetPoint(j,factor*j,trace_start1->at(j));
    }
  }
  if(size1!=0){
    for(UInt_t j=0;j<size1;j++){
      fTraces[1]->SetPoint(j,factor*(j+phase[0]),trace_start2->at(j));
    }
  }
  if(size2!=0){
    for(UInt_t j=0;j<size2;j++){
      fTraces[2]->SetPoint(j,factor*j,trace_stop1->at(j));
    }
  }
  if(size3!=0){
    for(UInt_t j=0;j<size3;j++){
      fTraces[3]->SetPoint(j,factor*(j+phase[1]),trace_stop2->at(j));
    }
  }

  if(size0==0&&size1==0&&size0==0&&size3==0)
    return ;
  else{
    if(draw){
      TCanvas *c;
      if(!gPad){
	c=new TCanvas();
	c->Divide(2,2); 
      }
      else{
	//  c=(TCanvas*)gPad;
	c=gPad->GetCanvas();
	c->Clear();
	c->Divide(2,2);     
      }
      for(Int_t i=0;i<4;i++){
	c->cd(i+1);   
	if(fTraces[i]){
	  fTraces[i]->Draw("AL*");
	}
	else
	  continue;
      }
    }

    fXCorrelation[0]->Draw("*AL");
    fXCorrelation[1]->SetLineColor(kBlue);
    fXCorrelation[1]->SetMarkerColor(kBlue);
    fXCorrelation[1]->Draw("*L");
    return;
  }
  
}

void  newTimingClass::Fit(Long64_t entry, Bool_t fix){
  
  Plot(entry,false);
  TF1 *fit[4];


  
  //Double_t phase[4]={start1_phase,start2_phase,stop1_phase,stop2_phase};
  Double_t phase[4]={start1_phase/4.,start2_phase/4.,stop1_phase/4.,stop2_phase/4.};
  Double_t qdc[4]={start1_qdc,start2_qdc,stop1_qdc,stop2_qdc};
  Double_t baseline[4]={start1_abase,start2_abase,stop1_abase,stop2_abase};
  Double_t timestamp[4]={StartTimeStamp[0],StartTimeStamp[1],StopTimeStamp[0],StopTimeStamp[1]};
  Int_t maximum[4]={StartMaximum[0],StartMaximum[1],StopTimeMaximum[0],StopTimeMaximum[1]};
  Double_t delta[4]={3.5,3.5,3.5,3.5};
  //Old Boards
  //Double_t fit_gamma[4]={0.2641,0.2675,0.121,0.1009};
  //Double_t fit_beta[4]={0.0973,0.1281,0.01008,0.1424};
  //New Boards
  //Double_t fit_beta[4]={0.1162,0.117,0.1178,0.119};
  //Double_t fit_gamma[4]={0.2504,0.2611,0.2611,0.2611};
  //UtkScan
  Double_t fit_beta[4]={0.12,0.12,0.12,0.12};
  Double_t fit_gamma[4]={0.245,0.245,0.245,0.245};
  vector<pair<Int_t,Int_t>>fit_limits;
  fit_limits.push_back(make_pair(10,10));
  fit_limits.push_back(make_pair(10,10));
  fit_limits.push_back(make_pair(10,10));
  fit_limits.push_back(make_pair(10,10));
  //Double_t delta[4]={4,4,4,4};
  Char_t fname[256];
  for(Int_t m=0;m<4;m++){
    ROOT_phase[m]=-100;
    ROOT_beta[m]=-100;
    ROOT_gamma[m]=-100;
    beta[m]=fit_beta[m];
    gamma[m]= fit_gamma[m];
    ROOT_qdc[m]=qdc[m];
    vector <UInt_t>::iterator it;
    UInt_t max_position=0;
    switch(m){
    case 0:
      if(trace_start1->size()!=0){
	it=max_element(trace_start1->begin(),trace_start1->end());
	max_position=distance(trace_start1->begin(),it);
      }
      break;
    case 1:
      if(trace_start2->size()!=0){
	it=max_element(trace_start2->begin(),trace_start2->end());
	max_position=distance(trace_start2->begin(),it);
      }
      break;
    case 2:
      if(trace_stop1->size()!=0){
	it=max_element(trace_stop1->begin(),trace_stop1->end());
	max_position=distance(trace_stop1->begin(),it);
      }
      break;
    case 3:
      if(trace_stop2->size()!=0){
	it=max_element(trace_stop2->begin(),trace_stop2->end());
	max_position=distance(trace_stop2->begin(),it);
      }
      break;
    default:
      break;
    }
    // cout<<"Maximum is "<<max_position*4<<endl;
    
    if(fTraces[m]->GetN()>0){
      std::pair <Double_t,Double_t> range((max_position-fit_limits[m].first)*factor,(max_position+fit_limits[m].second)*factor);
      //fit[m] =new TF1(fname,sipmfunc,range.first,range.second,6);
      fFits[m]->SetRange(range.first,range.second);
      fFits[m]->SetParameters(range.first,qdc[m]*0.5,beta[m],gamma[m],baseline[m],delta[m]);
      fFits[m]->SetParNames("#phi","#alpha","#beta","#gamma","Baseline","#delta");
      fFits[m]->FixParameter(4,baseline[m]);
      if(fix){
	fFits[m]->FixParameter(2,beta[m]);
	fFits[m]->FixParameter(3,gamma[m]);
      }      
      fFits[m]->FixParameter(5,delta[m]);
      TFitResultPtr status=fTraces[m]->Fit(fFits[m],"RNQSW");
      fStatus[m] =status->Status(); 
      if(status->IsValid()){
	//status->Print();
	ROOT_phase[m]=fFits[m]->GetParameter(0)-range.first+max_position;
	ROOT_beta[m]=fFits[m]->GetParameter(2);
	ROOT_gamma[m]=fFits[m]->GetParameter(3);
	ROOT_chisq[m]=fFits[m]->GetChisquare();
	//ROOT_time[m]=fFits[m]->GetParameter(0)+timestamp[m];
	ROOT_time[m]=ROOT_phase[m]*4+timestamp[m];
	ROOT_max[m]=maximum[m];
      }
      
      //cout<<"PAASS phase "<<phase[m]<<"\t ROOT phase "<<ROOT_phase[m]<<endl;
      
    }
  }
}


void  newTimingClass::Residuals(Long64_t entry,Bool_t fix){
 
  Fit(entry,fix);

  TCanvas *c;
  if(!gPad){
    c=new TCanvas();
    c->Divide(2,2);
  }
  else{
    c=gPad->GetCanvas();
    c->Clear();
    c->Divide(2,2);
  }   
  for(Int_t i=0;i<4;i++){
    fResiduals[i]->Set(0);
    c->cd(i+1);   
    if(fTraces[i]->GetN()>0){
      //fFits[i]->SetRange(0,fTraces[i]->GetN()*4);
      Double_t *yval=fTraces[i]->GetY();
      Double_t *xval=fTraces[i]->GetX();
      Int_t n=0;
      for(Int_t m=0;m<fTraces[i]->GetN();m++){
	//cout<<xval[m]<<" "<<yval[m]<<endl;
	if(fFits[i]->Eval(factor*m)!=0){
	  fResiduals[i]->SetPoint(n,xval[m],yval[m]-fFits[i]->Eval(factor*m));
	  n++;
	}
      }
      //fTraces[i]->Draw("AL*");
      fResiduals[i]->Draw("AL*");
      //fFits[i]->Draw("same");
    }
    else
      continue;
  }
}

void  newTimingClass::Trap_filter(Long64_t entry,UInt_t length,UInt_t gap){
  
  GetEntry(entry);
  if(gPad)
    gPad->GetCanvas()->Clear();
  
  Int_t *filter=new Int_t[2*length+gap];
  for(UInt_t j=0;j<length;j++){
    filter[j]=1;
    if(gap>0)
      filter[j+gap]=0;
    filter[j+gap+length]=-1;
  }
  
  cout<<"size "<<trace_start1->size()<<endl;
  
  if(trace_start1->size()!=0){
    
    //Int_t *filter_output=new Int_t[trace_start1->size()];
    Double_t *filter_output=new Double_t[trace_start1->size()];
    TGraph *g=new TGraph();
    g->SetLineColor(2);
    g->SetLineWidth(2);
    TGraph *g2=new TGraph();
    g2->SetLineColor(4);
    g2->SetLineWidth(2);
    Int_t n_points=0;
    for(UInt_t n=2*length+gap;n<trace_start1->size();n++){
      Int_t sum=0;
      //cout<<"N "<<n<<endl; 
      for(UInt_t j=0;j<2*length+gap;j++){
	//cout<<filter[j]<<endl;
	if(n>j)
	  sum+=filter[j]*trace_start1->at(n-j);     
	//cout<<n<<" "<<j<<" "<<" "<<filter[j]<<" "<<sum<<endl;
      }
      filter_output[n]=(Double_t)sum/length;
      //cout<<n<<" "<<sum<<" "<<filter_output[n]<<endl;
      g->SetPoint(n_points,n,filter_output[n]);
      g2->SetPoint(n_points,n,trace_start1->at(n));
      n_points++;
    }
    g->Draw("AL");
    g2->Draw("L");
  }
  return;
}


Int_t newTimingClass::CrossCorrelation(vector<unsigned int> *trace1,vector<unsigned int> *trace2, TGraph *fGraph){

  // cout<<"HERE"<< trace1->size()<<" "<<trace2->size()<<endl;

  Int_t offset= trace1->size();
  //  cout << "Trace Length: " << offset << endl;
  Int_t max=-9999;
  Int_t max_bin=-9999;

  Int_t minbin[2]= {0};
  UInt_t basemin[2] = {trace1->at(0),trace2->at(0)};
  basemin[0] = *std::min_element(std::begin(*trace1),std::end(*trace1));
  basemin[1] = *std::min_element(std::begin(*trace2),std::end(*trace2));
   
  cout << "Minimums: t1(" << basemin[0] << ") t2(" << basemin[1] << ")" << endl;  
	///Subtracting the minimum for better normalization

  UInt_t traceint1 = 0, traceint2 = 0;
  vector<double> norm1(200);
  vector<double> norm2(200);
  
  for (int in = 0; in<offset; in++){
	norm1[in] = trace1->at(in)-basemin[0]; if(in>0) traceint1 += 0.5*(norm1[in]+norm1[in-1]);
	norm2[in] = trace2->at(in)-basemin[1]; if(in>0) traceint2 += 0.5*(norm2[in]+norm2[in-1]);
   }

   cout << "Normalization Constants: t1=" << traceint1 << " t2=" << traceint2 << endl;


 for (int is=0; is<offset; is++){
	norm1[is] = trace1->at(is)*1.0/traceint1;
        norm2[is] = trace2->at(is)*1.0/traceint2;
   } 
  

  for(Int_t n=-20;n<20;n++){
    //      cout<<n<<" "<<offset<<endl;  
    Int_t conv=0;

    for(Int_t m=0;m<offset;m++){
      
      if(m>=abs(n)&&(m<200-abs(n))&&m<200){
	conv+=trace1->at(m)*trace2->at(m+n)/1000000.;
	//	cout<<"trace "<<m<<" "<<m+n<<" "<<trace2->at(m+n)<<" "<<n<<endl;
      }
    }
    if(max<conv){
      max=conv;
      max_bin = n;
    }
    fGraph->SetPoint(n+20,n,conv);
    
  }
 
  return max_bin;
}

//void Normalize(vector<unsigned int> *trace1,vector<unsigned int> *trace2){
//  /// First find the minimum subtract baseline
//
//  Int_t size1 = trace1->size();
//  Int_t size2 = trace2->size();
//
//  Int_t minbin[2]= {0};
//  Int_t basemin[2] = {trace1->at(0),trace2->at(0)};
//  basemin[0] = *std::min_element(trace1,trace1+size1);
//  basemin[1] = *std::min_element(trace2,trace2+size2);
//   
//
//
////  for (int ij = 1; ij<size1 || ij<size2; ij++){
////	if (basemin[0]>trace1->at(ij)){ basemin[0]=trace1->at(ij); minbin[0] = ij;}
////	if (basemin[1]>trace2->at(ij)){ basemin[1]=trace2->at(ij); minbin[1] = ij;}
////    }
//  cout << "Minimums: t1(" << basemin[0] << ") t2(" << basemin[1] << ")" << endl;  
//	///Subtracting the minimum for better normalization
//
//  UInt_t traceint1 = 0, traceint2 = 0;
//
//  for (int in = 0; in<size1 || in<size2; in++){
//	trace1[in] = trace1->at(in)-basemin[0]; if(in>0) traceint1 += 0.5*(trace1[in]+trace1[in-1]);
//	trace2[in] = trace2->at(in)-basemin[1]; if(in>0) traceint2 += 0.5*(trace2[in]+trace2[in-1]);
//   }
//
//   cout << "Normalization Constants: t1=" << traceint1 << " t2=" << traceint2 << endl;
//
// for (int is=0; is<size1 || is<size2; is++){
//	trace1[is] = trace1[is]*1./traceint1;
//        trace2[is] = trace2[is]*1./traceint2;
//   } 
//
//}

#endif // #ifdef newTimingClass_cxx
