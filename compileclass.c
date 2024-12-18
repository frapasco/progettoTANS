void Compileclass (TString myopt="fast"){
  TString opt;
  if(myopt.Contains("force")){opt="kfg";} else {opt ="kg"; }
  
  gSystem->SetBuildDir("./build",true);
  gSystem->AddIncludePath("-I./include");
  gSystem->CompileMacro("./src/Vertex.cxx", opt.Data());
  gSystem->CompileMacro("./src/Track.cxx", opt.Data()); 
  gSystem->CompileMacro("./src/Hit.cxx", opt.Data()); 
  gSystem->CompileMacro("./macros/simulation.cxx", opt.Data()); 
  gSystem->CompileMacro ("./macros/reconstruction.cxx", opt.Data()); 
  gSystem->CompileMacro ("./macros/canvasses.cxx", opt.Data());
}

void Clean(TString myopt="build"){
  gSystem->Exec("rm -rf ./build");
  if(myopt.Contains("all")){
    gSystem->Exec("rm -rf ./data.txt");
    gSystem->Exec("rm -rf ./simulation.root");
    gSystem->Exec("rm -rf ./reconstruction.root");
    gSystem->Exec("rm -rf ./histo.root");
  }
}
