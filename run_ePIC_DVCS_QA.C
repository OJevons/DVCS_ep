// Running macro for ePIC DVCS analysis
#include "ePIC_DVCS_QA.cxx"
const float fMass_proton{0.938272};
const Float_t fMass_electron{0.000511};

void run_ePIC_DVCS_QA(TString camp="Camp", TString energy="10x100", TString sett="test", TString comment="X"){

  std::cout<<"----------------------------------------------------"<<std::endl;
  std::cout<<"                 ePIC DVCS Analysis                 "<<std::endl;
  std::cout<<"----------------------------------------------------"<<std::endl;
  std::cout<<std::endl;

  // Initialize DVCS analysis object
  std::cout<<"Settings:"<<std::endl;
  std::cout<<"\tCampaign - "<<camp<<std::endl;
  std::cout<<"\tBeam energy - "<<energy<<" GeV"<<std::endl;
  std::cout<<"\tBeam setting - "<<sett<<std::endl;

  ePIC_DVCS_TASK *objDVCS = new ePIC_DVCS_TASK(camp,energy,sett);
  
  TString sInFileList;
  // If testing, load test file list
  if(sett != "hiAcc" && sett != "hiDiv") sInFileList="./filelists/inputFileList_test.list";
  // Else, use file list as directed
  else sInFileList = "./filelists/inputFileList_ePIC_"+camp+"_"+energy+"_"+sett+".list";
  objDVCS->setInFileList(sInFileList);
  TString sOutFileName;
  if(comment=="X") sOutFileName = "$EIC_WORK_DIR/DVCS_Analysis/RootFiles/ePIC_DVCS_"+camp+"_"+energy+"_"+sett+".root";
  else sOutFileName = "$EIC_WORK_DIR/DVCS_Analysis/RootFiles/ePIC_DVCS_"+camp+"_"+energy+"_"+comment+".root";
 
  objDVCS->setOutFile(sOutFileName);

  // Calculate maximum expected t for Roman Pots from beam energies
  Int_t beamP{0};
  if(energy == "5x41") beamP = 41;
  else if(energy == "10x100") beamP = 100;
  else if(energy == "18x275") beamP = 275;
  else beamP = 100;
  Float_t momCross2 = TMath::Power(beamP,4) + TMath::Power(fMass_proton,4) - 2*TMath::Power(beamP,2)*TMath::Power(fMass_proton,2);
  Float_t momDot = TMath::Power(beamP,2) - TMath::Sqrt(momCross2)*TMath::Cos(5e-3);
  Float_t tRPmax = 2*(TMath::Power(fMass_proton,2) - momDot);

  // Set DVCS cut values
  objDVCS->setMin_Q2(1);                         // GeV^2
  objDVCS->setMax_tRP(TMath::Abs(tRPmax)*1.25);  // GeV^2
  //objDVCS->setMax_tRP(2);  // GeV^2
  objDVCS->setMax_M2miss(1);                     // GeV^2

  // Set other behaviours
  objDVCS->setUsePID(kFALSE);
  objDVCS->setUseExplicitMatch(kTRUE);
  objDVCS->setUseEventBeams(kFALSE);

  objDVCS->doAnalysis();
}

