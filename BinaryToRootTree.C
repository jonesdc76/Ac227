#include <bitset>
#include <string>
#include <vector>
#include <string>
#include <utility>
#include <stdio.h>
#include <tuple>
#include <iostream>
#include <fstream>

#include "TSystem.h"
#include "TTree.h"
#include "TMath.h"
#include "TFile.h"
#include "TROOT.h"
#include "TLeaf.h"
#include "TBranch.h"
#include "TString.h"
#include "TGraph.h"
#include "TF1.h"

////////////////////////////////////////////////////////////////
//BinaryToRootTree.C                                          //
//Purpose:  Convert binary file from Caen ADCs to ROOT trees. //
//Author:   Donald Jones                                      //
//Contact:  donald.jones@temple.edu                           //
//Date:     August 2016                                       //
//See README in this folder for information about contents    //
//of the ROOT tree.                                           //
//                                                            //
//************************************************************//
//                                                            //
//Program requires environment variables to be set.           //
//Add these lines to your .bashrc file in your home directory.//
//export ROOT_TREE_DIR=/home/prospect-collab/data/ROOTtrees   //
//export BINARY_DIR=/home/prospect-collab/data                //
//export CONFIG_DIR=dir_where_channel_map_is_located          //
//                                                            //
//************************************************************//
//                                                            //
//Usage:                                                      //
//BinaryToRootTree(int series, int file_num, int time_stamp,  //
//		     const char *set_name)                    //
//     series:  number of series data comes from              //
//   file_num:  file number in a given series                 //
// time_stamp:  UNIX time stamp of file creation date         //
//   set_name:  name given to data set eg. P50C etc           //
////////////////////////////////////////////////////////////////

using namespace std;
const double SECONDS_PER_CLOCK_CYCLE = 8.0e-9;
const double ADC_SAMPLE_RATE = 5.0e8;//ADC sample rate in Hz
const uint16_t N_PEDESTAL_SAMPLES = 10;//# of samples in pre-pulse pedestal calc

/////////////////////////////////////////////////////////////////////////
//Structure for holding header information included before each waveform.
//Headers consists of 4 4-byte words 
/////////////////////////////////////////////////////////////////////////
typedef struct header_t{
  uint16_t n_samples = 0;
  uint32_t channels=0;
  uint32_t counter=0;
  uint32_t timestamp=0;
}header_t;

typedef struct stat_t{
  double val = 0;
  double rms = 0;
  uint16_t n = N_PEDESTAL_SAMPLES;
}stat_t;

///////////////////////////////////////////////
//Return name assigned to given channel. 
///////////////////////////////////////////////
bool GetChannelNames(vector<TString> &chNames){
  chNames.resize(16);
  ifstream f(Form("%s/channel_map.config", gSystem->Getenv("CONFIG_DIR")));
  char name[255];
  string line;
  uint16_t ch, n = 0;
  getline(f, line);
  if(!(f.is_open() && f.good())){
    cout<<"Channel map configuration file not found. Exiting.\n";
    return false;
  }
  cout<<"Channel map configuration file found.\n";
  while(!f.eof()){
    f>>ch>>name;
    if(ch>chNames.size()-1)chNames.resize(ch + 1);
    cout<<"Channel "<<ch<<" mapped to "<<name<<endl;
    chNames[ch] = name;
    f.peek();
    ++n;
  }
  if(n<chNames.size()){
    cout<<"Some intermediate channel maps missing. May be problems ahead.\n";
  }
  return true;//return true if successful
}
///////////////////////////////////////////////


///////////////////////////////////////////////
//Read header which consists of 4 4-byte words 
///////////////////////////////////////////////
bool ReadHeader(FILE *f, header_t *h){
  uint32_t header[4];
  uint16_t size[2];

  //Read first 4-bytes as two 2-byte words
  fread(size, 2, 2, f);
  if(feof(f))return false;//return false if end of file reached

  //Each 14-bit digitized sample is given as a 2-byte word in binary file. 
  //Size in file header is given as number of 4-byte words including 16-byte 
  //header so subtract header size and multiply by 2 to get # of 2-byte samples 
  h->n_samples = (size[0] - 4) * 2;

  //This should always be 0xa000 so check it each time to make sure data
  //is being read properly
  if(size[1] != 0xa000){
    cout<<"Failed header read check. Exiting. "<<size[1]<<"\n";
    return false;
  }

  //Read next 3 bytes as 4-byte words. 
  fread(header, 4, 3, f);
  if(feof(f))return false;//return false if end of file reached

  h->channels = header[0];
  h->counter= header[1];
  h->timestamp = header[2];

  return true;//return true if successful
} 
//////////////////////////////////////////////////



//////////////////////////////////////////////////////
//Convert Binary file to useable data in a ROOT tree
//////////////////////////////////////////////////////
int BinaryToRootTree(int series, int file_num, int time_stamp, 
		     const char *set_name){
  int verbose = 0;//set to different values for debugging


  //Assign channel names
  vector<TString>vChName;
  if(!GetChannelNames(vChName))
    return -1;

  //Open binary file
  //-----------------------------------------------------------------------
  TString filename = Form("s%03i_f%05i_ts%i", series, file_num, time_stamp);
  FILE *file = fopen(Form("%s/%s/series%03i/%s.dat",
			  gSystem->Getenv("BINARY_DIR"), set_name, series,
			  filename.Data()), "rb");
  if(file == NULL){
    printf("%s/%s/series%o3i/%s.dat",
	   gSystem->Getenv("BINARY_DIR"), set_name, series,
	   filename.Data());
    return -1;
  }else{
    cout<<filename.Data()<<".dat found"<<endl;
  }
  //-----------------------------------------------------------------------



  //Create Root tree and file for storage. 
  //-----------------------------------------------------------------------
  //File must be created first
  TFile *treeFile = new TFile(Form("%s/%s/series%03i/%s.root",
  				   gSystem->Getenv("ROOT_TREE_DIR"), set_name,
  				   series, filename.Data()), "recreate");
  if(treeFile==0){
    cout<<"Warning. New tree file failed to open."<<endl;
    return -1;
  }
  uint16_t nSamp = 0, nChan = 0;
  vector<uint16_t>vChan;

  //Read the first header to get number of channels and samples.
  header_t h;
  if(!ReadHeader(file, &h)) return -1;
  for(uint16_t i=0;i<32;++i){//check which channels are used
    if(h.channels & (1<<i)){
      cout<<"Adding Channel "<<i<<endl;
      vChan.push_back(i);
      ++nChan;
    }
  }
  nSamp = h.n_samples / nChan;
  cout<<nSamp<<" samples per event per channel"<<endl;

  if(verbose == 1 || verbose == 2){
    cout<<"0: "<<h.n_samples<<" "<<std::bitset<8>(h.channels);
    cout<<" "<<h.counter<<" "<<h.timestamp<<endl;
  }

  //Information for each event
  uint16_t *minIdx = new uint16_t[nChan], *min = new uint16_t[nChan];
  uint16_t *maxIdx = new uint16_t[nChan], *max = new uint16_t[nChan];
  uint16_t *halfMinIdx = new uint16_t[nChan], samples_per_wform = nSamp;
  double abs_time;
  //  bool good = true;

  //pedestal calculated from pre-pulse samples in waveform
  //ped1 from first n samples, ped2 from n+1 to 2n etc 
  stat_t *ped1 = new stat_t[nChan], *ped2 = new stat_t[nChan];
  stat_t *ped3 = new stat_t[nChan], *ped4 = new stat_t[nChan];
  double *average = new double[nChan], dt = 0;
  uint16_t **wform = new uint16_t *[nChan];
  for(int i=0;i<nChan;++i){
    wform[i] = new uint16_t[nSamp];
  }
  

  TTree *tree = new TTree("tree","tree");
  //  TBranch *bGood = tree->Branch("good", &good, "good/O");
  TBranch *bAbsTime = tree->Branch("abs_time", &abs_time, "abs_time/D");
  TBranch *bDt = tree->Branch("dt", &dt, "dt/D");
  TBranch *bNsamp = tree->Branch("nSamples", &samples_per_wform, "nSamples/s"); 
  TBranch *bMin[nChan], *bMinIdx[nChan],  *bMax[nChan], *bMaxIdx[nChan]; 
  TBranch *bHalfMinIdx[nChan], *bAverage[nChan];
  TBranch *bPed1[nChan], *bPed2[nChan], *bPed3[nChan], *bPed4[nChan];
  TBranch *bWform[nChan];
  for(int ich=0;ich<nChan;++ich){
    bAverage[ich] = tree->Branch(Form("average%s", vChName.at(vChan[ich]).Data()),
				 &average[ich], 
			     Form("average%s/D", vChName.at(vChan[ich]).Data()));
    bMin[ich] = tree->Branch(Form("min%s", vChName.at(vChan[ich]).Data()),&min[ich], 
			     Form("min%s/s", vChName.at(vChan[ich]).Data()));
    bMinIdx[ich] = tree->Branch(Form("minIdx%s", vChName.at(vChan[ich]).Data()), 
				&minIdx[ich], 
				Form("minIdx%s/s", vChName.at(vChan[ich]).Data()));
    bHalfMinIdx[ich] = tree->Branch(Form("halfMinIdx%s", vChName.at(vChan[ich]).Data()), 
				    &halfMinIdx[ich], 
				    Form("halfMinIdx%s/s", vChName.at(vChan[ich]).Data()));
    bMax[ich] = tree->Branch(Form("max%s", vChName.at(vChan[ich]).Data()), &max[ich], 
			     Form("max%s/s", vChName.at(vChan[ich]).Data()));
    bMaxIdx[ich] = tree->Branch(Form("maxIdx%s", vChName.at(vChan[ich]).Data()), 
				&maxIdx[ich], 
				Form("maxIdx%s/s", vChName.at(vChan[ich]).Data()));
    bPed1[ich] = tree->Branch(Form("ped1%s", vChName.at(vChan[ich]).Data()), &ped1[ich],
			      "val/D:rms:n/s");
    bPed2[ich] = tree->Branch(Form("ped2%s", vChName.at(vChan[ich]).Data()), &ped2[ich],
			      "val/D:rms:n/s");
    bPed3[ich] = tree->Branch(Form("ped3%s", vChName.at(vChan[ich]).Data()), &ped2[ich],
			      "val/D:rms:n/s");
    bPed4[ich] = tree->Branch(Form("ped4%s", vChName.at(vChan[ich]).Data()), &ped2[ich],
			      "val/D:rms:n/s");
    bWform[ich] = tree->Branch(Form("%s_wform", vChName.at(vChan[ich]).Data()), 
			       wform[ich],
			       Form("%s_wform[nSamples]/s", vChName.at(vChan[ich]).Data()));
  }

  //-----------------------------------------------------------------------



  //Read all event waveforms and store in Root tree
  //-----------------------------------------------------------------------
  int n=1, nrollover = 0, prev_timestamp = h.timestamp;
  double prev_abs_time = 0, N = (double)N_PEDESTAL_SAMPLES;
  while(!feof(file)){
    if(n%50000==0)cout<<"Reading event "<<n<<endl;
    //Read the waveform(s)
    for(int ich=0; ich<nChan; ++ich){
      //initialize variables
      min[ich] = 1<<14;//larger than the 14-bit ADC range
      max[ich] = 0; 
      minIdx[ich] = 0;
      maxIdx[ich] = 0; 
      average[ich] = 0;
    
      double mean = 0, meanSq = 0, prev_mean = 0;
      for(int isamp=0; isamp<nSamp; ++isamp){
	uint16_t temp[1] = {0};
	fread(temp, 2 ,1, file);
	wform[ich][isamp] = temp[0];
	if(wform[ich][isamp] < min[ich]){
	  min[ich] = wform[ich][isamp]; 
	  minIdx[ich] = isamp;
	}
	if(wform[ich][isamp] > max[ich]){
	  max[ich] = wform[ich][isamp]; 
	  maxIdx[ich] = isamp;
	}
	average[ich] += wform[ich][isamp];

	if(isamp < 4 * N_PEDESTAL_SAMPLES){
	  if(isamp>=N_PEDESTAL_SAMPLES)
	    prev_mean = wform[ich][isamp-N_PEDESTAL_SAMPLES];
	  mean += wform[ich][isamp] - prev_mean;//rolling mean over latest N samples
	  meanSq +=  pow(wform[ich][isamp],2) - pow(prev_mean, 2);//rolling mean square
	  
	  switch(isamp){ 
	  case N_PEDESTAL_SAMPLES-1: 
	    ped1[ich].val = mean / N;
	    ped1[ich].rms = sqrt(meanSq / N - pow(mean / N, 2));
	  case 2*N_PEDESTAL_SAMPLES-1: 
	    ped2[ich].val = mean / N;			      
	    ped2[ich].rms = sqrt(meanSq / N - pow(mean / N, 2));
	  case 3*N_PEDESTAL_SAMPLES-1: 
	    ped3[ich].val = mean / N;			      
	    ped3[ich].rms = sqrt(meanSq / N - pow(mean / N, 2));
	  case 4*N_PEDESTAL_SAMPLES-1: 
	    ped4[ich].val = mean / N;			      
	    ped4[ich].rms = sqrt(meanSq / N - pow(mean / N, 2));
	  default: break;
	  }
	}
	if(wform[ich][isamp]>16383)return -1;
	if(verbose == 2){
	  if(isamp%8==0)cout<<endl<<isamp<<") Ch. "<<vChan[ich]<<": ";
	  printf("%6i ", wform[ich][isamp]);
	}
      }
      average[ich] /= (double)nSamp;
    }
    for(int ich=0; ich<nChan; ++ich){
      //use pedestal with smallest rms
      stat_t ped[4] = {ped1[ich], ped2[ich], ped3[ich], ped4[ich]};
      uint16_t iped = 0;
      for(int i=1;i<nChan;++i)
	if(ped[iped].rms < ped[i].rms)iped = i;
      int idx = (minIdx[ich]>0 ? minIdx[ich] - 1 : 0);
      while(wform[ich][idx] < 0.5*(min[ich] + ped[iped].val) && idx>0)--idx;
      halfMinIdx[ich] = idx;
    }

    if(verbose == 2)cout<<endl;

    tree->Fill();

    //Read next header
    header_t htemp;
    if(!ReadHeader(file, &htemp))break;
    if(verbose == 1 || verbose == 2){
      cout<<n<<": "<<htemp.n_samples<<" "<<std::bitset<8>(htemp.channels);
      cout<<" "<<htemp.counter<<" "<<htemp.timestamp<<endl;
    } 
    if(htemp.timestamp < (double)prev_timestamp)++nrollover;
    abs_time = htemp.timestamp + nrollover * (pow(2, 31)-1.0);

    abs_time *= SECONDS_PER_CLOCK_CYCLE; 
    dt = abs_time - prev_abs_time;
    
    prev_timestamp = htemp.timestamp;    
    prev_abs_time = abs_time;

    //Channel map and sample size should be the same for every event. 
    if(!(htemp.n_samples == h.n_samples && htemp.channels == h.channels)){
      cout<<"Configuration change in middle of run not allowed. Exiting.\n";
      return -1;
    }

    ++n;
  }
  //-----------------------------------------------------------------------

  cout<<"Total Number of Events: "<<n<<endl;
  tree->Write("",TObject::kOverwrite);
  fclose(file);
  treeFile->Close();
  delete treeFile;
  return 0;
}
