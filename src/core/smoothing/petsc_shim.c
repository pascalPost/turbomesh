// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

#include <petsc.h>

MPI_Comm get_petsc_comm_self() { return PETSC_COMM_SELF; }
