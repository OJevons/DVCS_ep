// Running macro for ePIC DVCS analysis
#include "ePIC_DVCS_TASK.cxx"

void run_ePIC_DVCS(TString camp="Camp", TString energy="10x100", TString sett="test"){

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
  
  TString sInFileList = "./filelists/inputFileList_ePIC_"+camp+"_"+energy+"_"+sett+".list";
  objDVCS->setInFileList(sInFileList);
  TString sOutFileName = "$EIC_WORK_DIR/DVCS_Analysis/RootFiles/ePIC_DVCS_"+camp+"_"+energy+"_"+sett+".root";
  objDVCS->setOutFile(sOutFileName);

  // Set DVCS cut values
  objDVCS->setMin_Q2(1);     // GeV^2
  objDVCS->setMax_tRP(0.3);  // GeV^2
  objDVCS->setMax_M2miss(1); // GeV^2

  // Set other behaviours
  objDVCS->setUsePID(kFALSE);
  objDVCS->setUseExplicitMatch(kTRUE);

  objDVCS->doAnalysis();
}

