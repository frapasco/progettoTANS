#include "TApplication.h"
#include <iostream>

#include "./include/Track.h"
#include "./include/Hit.h"
#include "./include/Vertex.h"
#include "./include/simulation.h"

int main(int argc,char *argv[]){
  TApplication *myApp = new TApplication("myApp",&argc,argv); 
  for(int i=0;i<10;i++)
    std::cout<<"test "<<i<<std::endl;
  myApp->Run();
  simulazione();
  delete(myApp);
  return 0;
}
