
// header files
// c++ headers
#include <iostream>
using namespace std;

// root headers
#include "TGrid.h"
#include "TGridResult.h"
#include "TFileMerger.h"

// main program
void GetFiles(Int_t opt = 0)
// get a specifc set of files from grid
{

  //connect to the GRID
  TGrid::Connect("alien://");

  TGridResult* result = NULL;
  if (opt == 0) result = gGrid->Query("/alice/cern.ch/user/j/jgcn/LHC16r_kIncohJpsiToMu_fwd/myOutputDir/","AnalysisResults.root");
  if (opt == 1) result = gGrid->Query("/alice/cern.ch/user/j/jgcn/LHC16s_kIncohJpsiToMu_fwd/myOutputDir/","AnalysisResults.root");
  TFileMerger m;
  //  if (opt == 0) m.OutputFile("AnalysisResults_LHC16r_fwd_kIncohJpsiToMu_nopdca.root");
  //  if (opt == 1) m.OutputFile("AnalysisResults_LHC16s_fwd_kIncohJpsiToMu_nopdca.root");
  if (opt == 0) m.OutputFile("AnalysisResults_LHC16r_fwd_kIncohJpsiToMu.root");
  if (opt == 1) m.OutputFile("AnalysisResults_LHC16s_fwd_kIncohJpsiToMu.root");
   
  Int_t i = 0;
  //Loop over the TGridResult entries and add them to the TFileMerger
  while(result->GetKey(i,"turl")) {
    cout << endl << endl;
    cout << " adding " << result->GetKey(i,"turl") << endl;
    m.AddFile(result->GetKey(i,"turl"));
    i++;
  }

  cout << endl << endl;
  cout << endl << endl;
  cout << endl << endl;

  cout << " ********** MERGING ************ " << endl << endl << endl;
  //Merge
  m.Merge();

}

	
	
	
