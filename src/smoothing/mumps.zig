// build with
// // mumps
// tests.linkSystemLibrary("dmumps");
// tests.linkSystemLibrary2("mumps_common", .{ .needed = true });
// tests.linkSystemLibrary2("cmumps", .{ .needed = true });
// tests.linkSystemLibrary2("esmumps", .{ .needed = true });
// tests.linkSystemLibrary2("ptesmumps", .{ .needed = true });
// tests.linkSystemLibrary2("gfortran", .{ .needed = true });
// tests.linkSystemLibrary2("lapack", .{ .needed = true });
// tests.linkSystemLibrary2("scalapack", .{ .needed = true });
// tests.linkSystemLibrary2("blas", .{ .needed = true });
// tests.linkSystemLibrary2("pord", .{ .needed = true });
// tests.linkSystemLibrary2("scotch", .{ .needed = true });
// tests.linkSystemLibrary2("scotcherrexit", .{ .needed = true });
// tests.linkSystemLibrary2("metis", .{ .needed = true });
// tests.linkSystemLibrary2("mpi_usempif08", .{ .needed = true });
// tests.linkSystemLibrary2("mpi_mpifh", .{ .needed = true });
// tests.linkSystemLibrary2("mpi", .{ .needed = true });

const std = @import("std");
const c = @cImport({
    @cInclude("dmumps_c.h");
    @cInclude("mpi.h");
});

// extern const MPI_COMM_WORLD: c.ompi_mpi_comm_world;
// extern var MPI_COMM_WORLD_test: opaque {};
// const p_external_variable = &MPI_COMM_WORLD_test;
// const MPI_COMM_WORLD: c.MPI_Comm = MPI_COMM_WORLD_ptr.*;

// const MPI_COMM_WORLD = &c.ompi_mpi_comm_world;
// const MPI_COMM_WORLD :  = &c.ompi_mpi_comm_world;

test "mumps" {
    // const MPI_COMM_WORLD: c.MPI_Comm = @as(*c.MPI_Comm, @ptrCast(@alignCast(p_external_variable))).*;

    // const tes = c.MPI_COMM_WORLD();
    // _ = tes;

    if (c.MPI_Init(null, null) != c.MPI_SUCCESS) {
        return error.MpiInitFailed;
    }

    var rank: c_int = undefined;
    if (c.MPI_Comm_rank(@ptrCast(&c.ompi_mpi_comm_world), &rank) != c.MPI_SUCCESS) {
        return error.MpiCommRankFailed;
    }

    var rhs = [_]f64{ 1, 4 };
    var a = [_]f64{ 1, 2 };
    var irn = [2]c.MUMPS_INT{ 1, 2 }; // nonzero row indices
    var jcn = [_]c.MUMPS_INT{ 1, 2 }; // nonzero col indices

    var mumps_par = c.DMUMPS_STRUC_C{};

    mumps_par.sym = 0; // unsymmetric matrix
    mumps_par.par = 1; // host involved  in parallel steps of factorization and solve
    mumps_par.comm_fortran = -987654; // magic constant representing MPI_COMM_WORLD

    // init a mumps instance
    mumps_par.job = -1; // mumps job init
    c.dmumps_c(&mumps_par);

    // solve

    if (rank == 0) {
        mumps_par.n = 2; // order of the matrix
        mumps_par.nnz = 2; // nonzeros in the matrix
        mumps_par.irn = irn[0..];
        mumps_par.jcn = jcn[0..];
        mumps_par.a = a[0..];
        mumps_par.rhs = rhs[0..];
    }

    // mumps_par.icntl[1] = -1;
    mumps_par.icntl[2] = -1;
    // mumps_par.icntl[3] = -1;
    // mumps_par.icntl[4] = 0;

    mumps_par.job = 6; // mumps job analysis + factorization + solve
    c.dmumps_c(&mumps_par);

    // terminate
    mumps_par.job = -2; // mumps job end
    c.dmumps_c(&mumps_par);

    if (rank == 0) {
        std.debug.print("solution is {} {}\n", .{ rhs[0], rhs[1] });
    }

    if (c.MPI_Finalize() != c.MPI_SUCCESS) {
        return error.MpiFinalizeFailed;
    }
}
