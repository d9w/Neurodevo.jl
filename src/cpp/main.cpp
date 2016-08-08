#include <iostream>
#include <cstdio>
#include <time.h>
#include <mpi.h>
#include "easylogging++.h"
#include "color.h"
#include "GRN.h"
#include "GA.h"

INITIALIZE_EASYLOGGINGPP

int main(int argc, char** argv)
{
  MPI_Init(NULL, NULL);
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  srand(time(NULL));

  el::Configurations conf("ga_log.conf");

  //conf.set(el::Level::Info, el::ConfigurationType::Filename,
  //         "logs/ga_mnist.log");
  el::Loggers::reconfigureAllLoggers(conf);

  GA ga("grns", world_size, world_rank);
  ga.run();

  MPI_Finalize();
  return 0;
}
