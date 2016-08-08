#include <iostream>
#include "GRN.h"

int main(int argc, char** argv)
{
  GRN grn = GRN(argv[1]);
  std::cout << grn.toString() << std::endl;
}
