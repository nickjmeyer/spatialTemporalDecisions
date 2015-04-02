#include "test.hpp"
#include <omp.h>

int main(int argc, char ** argv){
  njm::sett.set(argc,argv);  

  njm::sett.clean();
  return 0;
}
