#include <petsc.h>

MPI_Comm get_petsc_comm_self() { return PETSC_COMM_SELF; }
