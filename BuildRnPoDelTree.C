#include <iostream>
#include <list>
#include <vector>

#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TString.h"
#include "TSystem.h"

using namespace std;

//half life of Po215 in s
const double kPo215HalfLife = 1.78e-3;

//earliest time (s) to look for correlated partner 
const double kStart = 10e-6; 

//number of half lives in near window
const int nHalf = 3;

//latest time (ns) to look for correlated partner
const double kEnd = kStart + nHalf * kPo215HalfLife; 

//number of half lives between end of near window and start of far
const int nHalfSep = 10;

//start time of far window for finding accidental correlated partners
const double kFarStart = kEnd + nHalfSep * kPo215HalfLife; 

//how much more statistics to accumulate in far window than near
const int kScale = 4;

//latest time to look for correlated partner
const double kFarEnd = kFarStart + kScale * (kEnd - kStart); 

//cut on total charge for delayed alphas
const double kQlow = 0.01, kQhigh = 0.05;

struct Pulse_t{
  Double_t t;
  Double_t Qtotal;
  Double_t Qtail;
  Double_t average;
  Double_t ped;
  Double_t PSD;
  uint16_t min;
  uint16_t minIdx;
};

struct RnPo_t{
  int size; //if prompt filled size=1, prompt and delayed filled>=2,
  int mult_p;//multiplicity of prompt candidate pulses
  int mult_f;//multiplicity of far pulses
  vector<Pulse_t> prompt;
  Pulse_t delayed;
  vector<Pulse_t> far;
};

//Builds tree of RnPo candidates with each entry a different single delayed 
//candidate each with any number of associated prompt and far events  

int BuildRnPoDelTree(TString fname, TString dir = "/home/jonesdc76/PROSPECT/data/Ac227/") 
{
  Pulse_t nullevt;
  nullevt.t = -99999;
  nullevt.Qtotal = -99999;
  nullevt.Qtail = -99999;
  nullevt.average = -99999;
  nullevt.ped = -99999;
  nullevt.PSD = -99999;
  nullevt.min = 65535;
  nullevt.minIdx = 65535;
   //Open the Root tree
  //////////////////////////
  TString fn = dir + fname;
  TFile *file = TFile::Open(fn.Data());
  if(file==0){
    cout<<"Cannot open file "<<fn.Data()<<". Exiting."<<endl;
    return -1;
  }
  TTree *tree = (TTree*)file->Get("tree");

  if(tree == 0){
    cout<<"Cannot find requested ROOT tree in "<<fn.Data()<<". Exiting."<<endl;
    return -1;
  }


  //parse time stamp
  //////////////////
  TString timestamp = fname;
  timestamp.Remove(timestamp.Last('.'),5);
  timestamp.Remove(0, timestamp.First('_')+3);
  cout<<"timestamp: "<<timestamp.Data()<<endl;
  if(!timestamp.IsDigit()){
    cout<<"Problem parsing timestamp: "<<timestamp.Data()<<" is not a digit\n";
    FILE *flog = fopen("/home/jonesdc/prospect/time.log", "a+");
    fprintf(flog, "%s", fname.Data());
    fclose(flog);
  }
  double ts = (double)timestamp.Atoll();


  //set up branch address pointers
  ///////////////////////////////////////
  Pulse_t pulse, prev_pulse; 
  uint16_t good, nSamples;
  tree->SetBranchAddress("abs_time", &pulse.t);
  tree->SetBranchAddress("psdCh0", &pulse.PSD);
  tree->SetBranchAddress("minCh0", &pulse.min);
  tree->SetBranchAddress("minIdxCh0", &pulse.minIdx);
  tree->SetBranchAddress("QtotalCh0", &pulse.Qtotal);
  tree->SetBranchAddress("QtailCh0", &pulse.Qtail);
  tree->SetBranchAddress("pedCh0", &pulse.ped);
  tree->SetBranchAddress("averageCh0", &pulse.average);
  tree->SetBranchAddress("goodCh0", &good);

  //Create new ROOT tree to hold BiPo Events
  //////////////////////////////////////////
  fn = dir + "RnPoDel_" + fname;
  TFile *newfile = TFile::Open(fn.Data(), "RECREATE");
  if(newfile == 0){
    cout<<"Cannot open file "<<fn.Data()<<". Exiting."<<endl;
    return -1;
  }

  Int_t nFar = 0, nPrompt = 0;
  Pulse_t Delayed;
  RnPo_t rnpo;
  vector<Pulse_t>vPrompt, vFar;
  TTree *RnPoTree = new TTree("RnPo", "Tree containing Rn-Po candidate events");
  RnPoTree->Branch("nFar", &nFar,"nFar/I");
  RnPoTree->Branch("nPrompt", &nPrompt,"nPrompt/I");
  RnPoTree->Branch("far", &vFar);
  RnPoTree->Branch("prompt", &vPrompt);
  RnPoTree->Branch("delayed", &Delayed);

  //function used to divide between between recoil and ionization bands
  TF1 *f = new TF1("f","[0]+[1]*pow(0.5,(x-[3])*[2])",0,0.1);
  f->SetParameters(0.31,0.12,500,0.012);

  bool pass = false;
  int n = 0;
  Long64_t N = tree->GetEntries();
  for(int i=0;i<N;++i){
    //    cout<<i<<" "<<N<<endl;
    vFar.clear();
    vPrompt.clear();
    rnpo.prompt.clear();
    rnpo.far.clear();
    rnpo.size = 0;
    rnpo.mult_p = 0;
    rnpo.mult_f = 0;

    tree->GetEntry(i);
    //    cout<<"time: "<<pulse.t<<" psd: "<<pulse.PSD<<" min: "<<pulse.min<<" minIdx: "<<pulse.minIdx<<" Qtot: "<<pulse.Qtotal<<" Qtail: "<<pulse.Qtail<<" ped: "<<pulse.ped<<" avg: "<<pulse.average<<" good: "<<good<<endl;

    //Identify delayed candidates
    if(!(pulse.Qtotal > kQlow //lower charge cut
       && pulse.Qtotal < kQhigh //upper charge cut
	 && pulse.PSD > f->Eval(pulse.Qtotal)))//in recoil PSD band
      continue;

    rnpo.delayed = pulse;
    ++rnpo.size;

    //Ignore if multiple delayed candidates close enough that 
    //one might be counted twice once as prompt and once as delayed.
    // n = 1;
    // while(i+n < N){
    //   tree->GetEntry(i+n);
    //   ++n;
    //   //cout<<pulse.t<<endl;
    //   //stop if outside time window
    //   if(pulse.t - rnpo.delayed.t > (kEnd - kStart))break;

    //   //multiplicity cut to ensure no candidate is used twice
    //   //pass criteria?
    //   if(pulse.Qtotal > kQlow //lower charge cut
    // 	 && pulse.Qtotal < kQhigh //upper charge cut
    // 	 && pulse.PSD > f->Eval(pulse.Qtotal)){ //only alphas
    // 	i += n - 1;
    // 	goto reset;
    //   }
    // }


    //find prompt alpha candidates
    //////////////////////////////////
    n = 1;
    double dt = kEnd - kStart;
    while(i-n >= 0){
      tree->GetEntry(i-n);
      ++n;

      //continue if not to start window yet
      if((rnpo.delayed.t - pulse.t) < kStart)continue;

      //beyond time window?
      if((rnpo.delayed.t - pulse.t) > kEnd)
	break; 

      //pass criteria?
      if(pulse.Qtotal > kQlow //lower charge cut
	 && pulse.Qtotal < kQhigh //upper charge cut
	 && pulse.PSD > f->Eval(pulse.Qtotal)){ //only alphas
	  rnpo.prompt.push_back(pulse);
	  ++rnpo.size;
	  ++rnpo.mult_p;
      }
    }

    //find time-displaced alpha candidates for background subtraction
    //////////////////////////////////////////////////////////////////
    n = 1;
    int nh = 0;
    while(i+n < N){
      tree->GetEntry(i+n);
      ++n;
      //continue if not yet to far time window
      if(pulse.t - rnpo.delayed.t < kFarStart) continue;
      //stop if past end of far time window
      if(pulse.t - rnpo.delayed.t > kFarEnd) break;
      
	
      //overlay multiple far windows to accumulate statistics
      pulse.t -= rnpo.delayed.t + (nHalf + nHalfSep) * kPo215HalfLife + kStart;
      int nh = int(pulse.t / dt);
      //cout<<nh<<endl;
      pulse.t -= nh * dt;
      pulse.t += kStart + rnpo.delayed.t;
      //pass criteria?
      if(pulse.Qtotal > kQlow //lower charge cut
	 && pulse.Qtotal < kQhigh //upper charge cut
	 && pulse.PSD > f->Eval(pulse.Qtotal)){ //only alphas
	rnpo.far.push_back(pulse);
	++rnpo.size;
	++rnpo.mult_f;
      }
    }
    if(rnpo.size < 2) continue;

    //fill values for branch pointers
    /////////////////////////////////
    nPrompt = rnpo.mult_p;
    nFar = rnpo.mult_f;
    Delayed = rnpo.delayed;
    for(n = (int)rnpo.prompt.size();n > 0; --n)
      vPrompt.push_back(rnpo.prompt[n-1]);
    if(vPrompt.size()==0)vPrompt.push_back(nullevt);
    for(n = 0; n < (int)rnpo.far.size(); ++n)
      vFar.push_back(rnpo.far[n]);
    if(vFar.size()==0)vFar.push_back(nullevt);
    
    RnPoTree->Fill();
    
  reset:
    //reset and clear
    n = 0;
  }
  cout<<RnPoTree->GetEntries()<<" total candidate RnPo events.\n";
  RnPoTree->Write("", TObject::kOverwrite);
  newfile->Close();
  file->Close();
  return 0;
}
