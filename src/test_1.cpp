#include <iostream>
#include <ctime>
#include "micro.h"

using namespace std;

int main (int argc, char *argv[])
{

  try {

    double start, end, t_assembly;
    Problem micropp (argc, argv);

    // assembly
    start = clock();
    micropp.assembly_A();
    end = clock();
    t_assembly = double(end - start) / CLOCKS_PER_SEC;

  } catch (int &e) {
    cerr << "Error : " << e << endl;
    return 1;
  }

  return 0;
}
