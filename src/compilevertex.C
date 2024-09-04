//funzione per compilare classi e macro, da eseguire in modalità interpretata


void compilevertex (TString myopt="fast"){
  TString opt;
  if(myopt.Contains("force")){opt="kfg";} else {opt ="kg"; }
  gSystem->CompileMacro ("Vertex.cxx", opt.Data());
  gSystem->CompileMacro ("Track.cxx", opt.Data()); 
  gSystem->CompileMacro ("Intpoint.cxx", opt.Data()); 
  gSystem->CompileMacro ("simulazione.C", opt.Data()); 
  gSystem->CompileMacro ("ricostruzione.C", opt.Data()); 
  gSystem->CompileMacro ("canvasses.C", opt.Data()); 
}
