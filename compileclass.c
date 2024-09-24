void compilevertex (TString myopt="fast"){
  TString opt;
  if(myopt.Contains("force")){opt="kfg";} else {opt ="kg"; }
  
  gSystem->SetBuildDir("./build",true);
  gSystem->AddIncludePath("-I./include");
  gSystem->CompileMacro("./src/Vertex.cxx", opt.Data());
  gSystem->CompileMacro("./src/Track.cxx", opt.Data()); 
  gSystem->CompileMacro("./src/Intpoint.cxx", opt.Data()); 
  gSystem->CompileMacro("./macros/simulation.cxx", opt.Data()); 
  //gSystem->CompileMacro ("./macros/reconstraction.C", opt.Data()); 
  //gSystem->CompileMacro ("./src/canvasses.C", opt.Data()); 
}
