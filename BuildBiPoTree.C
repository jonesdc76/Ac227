#include <iostream>
#include <list>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TString.h"
#include "TSystem.h"

using namespace std;

//UNIX timestamp in seconds of beginning of P50D/series007
const double kTStart = 1467830099;

//number of cells
const int kNCells = 2;

//maximum distance (mm) between Bi214 beta and Po-214 alpha
const int kMaxDist = 100;

//minimum energy threshold (MeV) in adjacent cell to be considered 
//multi-cell event. Below this may be light leakage
const double EThresh = 0.05;

//half life of Po214 in ns
const double kPo214HalfLife = 164000;

//earliest time (ns) to look for correlated partner 
const double kStart = 10000; 

//latest time (ns) to look for correlated partner
const double kEnd = kStart + 3 * kPo214HalfLife; 

//start time of far window for finding accidental correlated partners
const double kFarStart = kEnd + 10 * kPo214HalfLife; 

//factor by which far window is longer than near
const double kScale = 4.0; 

//latest time to look for correlated partner
const double kFarEnd = kFarStart + kScale * (kEnd - kStart); 



/// stackable! hit PID classifications that can be used in s_PhysPulse
enum HitTypeID {
  //electron/gamma/muon ionization region (might stack with neutron capture)
  IONI_HIT   = 1 << 0,   
     
  //neutron capture (on 6Li) region (might stack with gamma or recoil)
  NCAPT_HIT  = 1 << 1, 
  
  //nucleus recoil (fast neutron) clearly separated from gamma (might stack with neutron capture)
  RECOIL_HIT = 1 << 2, 
  
  //veto detector hit
  VETO_HIT   = 1 << 3,    
  
  //dead volume or non-detected hit (simulation only)    
  DEAD_HIT   = 1 << 4,        

  //unclassifiable event
  CRAZY_HIT  = 1 << 5,        

  //BiPo alpha region hit
  BIPO_HIT   = 1 << 6       
};


struct Pulse_t{
  Long64_t evt;
  Double_t t;
  Int_t seg;
  Float_t E;
  Float_t dt;
  Float_t y;
  Float_t PE[2];
  Float_t PSD;
  Int_t PID;
};

struct BiPo_t{
  int size; //if prompt filled size=1, prompt and delayed filled>=2,
  int mult_p;//multiplicity of prompt candidate pulses
  int mult_f;//multiplicity of far pulses
  vector<Pulse_t> prompt;
  Pulse_t delayed;
  vector<Pulse_t> far;
};



int BuildBiPoTree(TString fname, TString series = "007", 
		  TString data_set = "P50D",
		  TString data_release = "Phys_20161206", 
		  TString conv_data_dir = "/home/prospect-collab/data/converted_data") 
{
  Pulse_t nullevt;
  nullevt.evt = -99999;
  nullevt.seg = -99999;
  nullevt.E = -99999;
  nullevt.t = -99999;
  nullevt.dt = -99999;
  nullevt.y = -99999;
  nullevt.PE[0] = -99999;
  nullevt.PE[1] = -99999;
  nullevt.PSD = -99999;
  nullevt.PID = -99999;




  //Open the Root tree
  //////////////////////////
  TString fn = conv_data_dir + "/" + data_release + "/" + data_set + "/" + 
    "series" + series + "/" + fname;
  TFile *file = TFile::Open(fn.Data());
  if(file==0){
    cout<<"Cannot open file "<<fn.Data()<<". Exiting."<<endl;
    return -1;
  }
  TTree *tree = (TTree*)file->Get("PhysPulse");

  if(tree == 0){
    cout<<"Cannot find requested ROOT tree in "<<fn.Data()<<". Exiting."<<endl;
    return -1;
  }


  //parse time stamp
  //////////////////
  TString timestamp = fname;
  timestamp.Remove(timestamp.Last('_'),20);
  timestamp.Remove(0, timestamp.Length()-10);
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
  tree->SetBranchAddress("evt", &pulse.evt);
  tree->SetBranchAddress("seg", &pulse.seg);
  tree->SetBranchAddress("E", &pulse.E);  
  tree->SetBranchAddress("t", &pulse.t);
  tree->SetBranchAddress("dt", &pulse.dt);
  tree->SetBranchAddress("y", &pulse.y);
  tree->SetBranchAddress("PE", &pulse.PE);
  tree->SetBranchAddress("PSD", &pulse.PSD);
  tree->SetBranchAddress("PID", &pulse.PID);



  //Create new ROOT tree to hold BiPo Events
  //////////////////////////////////////////
  fn = gSystem->Getenv("P2X_ANALYZED");
  fn += "/" + data_release + "/" + data_set + "/" + 
    "series" + series; 
  //if the directory doesn't exist create it
  TString mkdircmd = "mkdir -p " + fn;
  system(mkdircmd);
  fn += "/biPoEvents_" + fname;
  TFile *newfile = TFile::Open(fn.Data(), "RECREATE");
  if(newfile == 0){
    cout<<"Cannot open file "<<fn.Data()<<". Exiting."<<endl;
    return -1;
  }

  Int_t nFar = 0, nPrompt = 0;
  Double_t abs_time;
  Pulse_t Delayed;
  vector<Pulse_t>vPrompt, vFar;
  TTree *BiPoTree = new TTree("BiPo", "Tree containing Bi-Po candidate events");
  BiPoTree->Branch("prompt", &vPrompt);
  BiPoTree->Branch("delayed", &Delayed, "evt/L:t/D:seg/I:E/F:dt/F:y/F:PE[2]/F:PSD/F:PID/I");
  BiPoTree->Branch("abs_time", &abs_time);//time in seconds since kTStart
  BiPoTree->Branch("far", &vFar);
  BiPoTree->Branch("nFar", &nFar, "nFar/I");
  BiPoTree->Branch("nPrompt", &nPrompt, "nPrompt/I");
  BiPoTree->Branch("abs_time", &abs_time, "abs_time/D");




  //loop over tree looking for correlated pairs
  /////////////////////////////////////////////
  BiPo_t bipo, prev_bipo;
  bipo.size = 0;
  bipo.mult_p = 0;
  bipo.mult_f = 0;
  prev_bipo = bipo;
  vector<BiPo_t> vBiPo;
  tree->GetEntry(0);
  long s_t = pulse.t;
  int excl_p = 0, excl_d = 0, excl_f = 0, n_p = 0, n_f = 0;
  int n = 1;
  bool pass = false;
  for(int i=1;i<tree->GetEntries();++i){
    prev_pulse = pulse;
    tree->GetEntry(i);
    int seg = pulse.seg;
    double pos = pulse.y;
    abs_time = pulse.t*1.0e-9 + ts - kTStart;
    if(i%100000==0 && 0)
      cout<<i<<": "<<pulse.evt<<" "<<pulse.seg<<" "<<pulse.t<<" "<<pulse.dt
	  <<" "<<pulse.E<<" "<<" "<<pulse.y<<" "<<pulse.PE[0]<<" "
	  <<pulse.PE[1]<<" "<<pulse.PID<<" "<<pulse.PSD<<endl;

    //look for nuclear recoil-like particles as alpha candidates
    if(!(pulse.PID & RECOIL_HIT))continue;
    bipo.delayed = pulse;
    ++bipo.size;


    //*******************************************************************
    //NOTE: no "continue" statements at this scope after this in for loop.
    //We are now searching for electron candidate since we have found a
    //recoil-like event and must reach end of loop naturally to reset
    //all variables and containers.
    //*******************************************************************

    n = 1;
    while(i-n >= 0 && i+n < tree->GetEntries()){
      //backward search for evidence of same event seen in other cell(s)
      tree->GetEntry(i-n);
      ++n;
      if(pulse.evt != bipo.delayed.evt)break;
      else if(pulse.seg != bipo.delayed.seg///exclude if seen in multiple cells
	&& pulse.E > EThresh){//with enough energy not just leakage
	++excl_d;
	goto reset;
      }
    }
    n = 1;
    while(i+n < tree->GetEntries()){
      //forward search for evidence of same event seen in other cell(s)
      tree->GetEntry(i+n);
      ++n;
      if(pulse.evt != bipo.delayed.evt)break;
      else if(pulse.seg != bipo.delayed.seg//exclude if seen in multiple cells
	      && pulse.E > EThresh){//with enough energy not just leakage
	++excl_d;
	goto reset;
      }
    }




    //find prompt beta candidates
    //////////////////////////////////
    n = 1;
    pass = false;
    while(i-n >= 0){
      prev_pulse = pulse;
      tree->GetEntry(i-n);
      ++n;

      //if passes all tests and not a multi-cell event then store
      if(pass){
	if(pulse.evt != prev_pulse.evt){//not same event as next entry
	  bipo.prompt.push_back(prev_pulse);
	  ++bipo.size;
	  ++bipo.mult_p;
	}else{
	  ++excl_p;
	}
      }

      //beyond time window?
      if((bipo.delayed.t - pulse.t) > kEnd)
	break; 

      //pass criteria?
      if(pulse.PID & IONI_HIT //only electrons
	 && pulse.seg == bipo.delayed.seg //in the same segment/cell
	 && abs(pulse.y - bipo.delayed.y) < kMaxDist //close to recoil position
	 && (bipo.delayed.t - pulse.t) > kStart){ //inside time window
	if(pulse.evt != prev_pulse.evt){//not same event as previous entry
	  pass = true;
	}else{
	  pass = false;
	  ++excl_p;
	}
      }else pass = false;
	 
    }
    

    //find betas in time-displaced window for background subtraction
    /////////////////////////////////////////////////////////////////
    if(bipo.size > 1){
      n = 1;
      pass = false;
      while(i+n < tree->GetEntries()){
	prev_pulse = pulse;
	tree->GetEntry(i+n);
	++n;
	//if passes all tests and not a multi-cell event then store
	if(pass){
	  if(pulse.evt != prev_pulse.evt){//not same event as next entry
	    bipo.far.push_back(pulse);
	    ++bipo.mult_f;
	  }else{
	    ++excl_f;
	  }
	}

	//beyond time window?
	if((pulse.t - bipo.delayed.t) > kFarEnd)break; 
	
	//pass criteria?
	if(pulse.PID & IONI_HIT //only electrons
	   && pulse.seg == bipo.delayed.seg //in the same segment/cell
	   && abs(pulse.y - bipo.delayed.y) < kMaxDist //close to recoil event
	   && (pulse.t - bipo.delayed.t) > kFarStart){ //inside far time window
	  if(pulse.evt != prev_pulse.evt){//not same event as previous entry
	    pass = true;
	  }else{
	    pass = false;
	    ++excl_f;
	  }
	}else pass = false;	 
      }
    

      //fill values for branch pointers
      /////////////////////////////////
      nPrompt = bipo.mult_p;
      nFar = bipo.mult_f;
      Delayed = bipo.delayed;
      for(n = (int)bipo.prompt.size();n > 0; --n)
	vPrompt.push_back(bipo.prompt[n-1]);
      for(n = (int)bipo.far.size();n > 0; --n)
	vFar.push_back(bipo.far[n-1]);
      if(vFar.size()==0)vFar.push_back(nullevt);
      
      if(0)
	cout<<i<<": "<<bipo.prompt[0].evt<<" "<<bipo.delayed.evt<<" "
	    <<(bipo.prompt[0].seg-bipo.delayed.seg)<<" "
	    <<(-bipo.prompt[0].t+bipo.delayed.t)*1e-3<<" "
	    <<abs(-bipo.prompt[0].y+bipo.delayed.y)<<" "
	    <<bipo.prompt[0].E<<" "<<bipo.delayed.E<<" "
	    <<bipo.prompt[0].PID<<" "<<bipo.delayed.PID<<endl;
      
      BiPoTree->Fill();
      n_f += nFar;
      n_p += nPrompt;
    }

    reset:
    //reset and clear
    vFar.clear();
    vPrompt.clear();
    bipo.prompt.clear();
    bipo.far.clear();
    bipo.size = 0;
    bipo.mult_p = 0;
    bipo.mult_f = 0;
  }

  cout<<BiPoTree->GetEntries()<<" events stored.\n";
  cout<<excl_p<<" events excluded for prompt events seen in multiple cells.\n";
  cout<<excl_d<<" events excluded for delayed events seen in multiple cells.\n";
  cout<<excl_f<<" events excluded for far events seen in multiple cells.\n";
  cout<<"Average ratio of far to prompt: "<<double(n_f)/double(n_p)<<endl;
  BiPoTree->Write("", TObject::kOverwrite);
  newfile->Close();
  file->Close();
  return 0;

}
