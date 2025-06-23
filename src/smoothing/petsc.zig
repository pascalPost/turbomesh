const std = @import("std");
const petsc = @import("petsc_import.zig");

test "petsc simple" {
    if (petsc.PetscInitialize(null, null, null, null) != petsc.PETSC_SUCCESS) {
        return error.PetscInitializeFailed;
    }

    // Create Matrix

    var A: petsc.Mat = undefined;

    if (petsc.MatCreate(petsc.PETSC_COMM_WORLD, &A) != petsc.PETSC_SUCCESS) {
        return error.MatCreateFailed;
    }

    const m: petsc.PetscInt = 2;
    const n: petsc.PetscInt = 2;
    _ = petsc.MatSetSizes(A, petsc.PETSC_DECIDE, petsc.PETSC_DECIDE, m, n);

    // matrix format (?)

    _ = petsc.MatSetFromOptions(A);

    // TODO: prealloc!

    _ = petsc.MatSetUp(A);

    // Set Matrix A = diag(1, 2)
    const rows = [_]petsc.PetscInt{ 0, 1 };
    const cols = [_]petsc.PetscInt{ 0, 1 };
    const vals = [_]petsc.PetscScalar{ 1.0, 2.0 };
    _ = petsc.MatSetValues(A, 1, &rows[0], 1, &cols[0], &vals[0], petsc.INSERT_VALUES);
    _ = petsc.MatSetValues(A, 1, &rows[1], 1, &cols[1], &vals[1], petsc.INSERT_VALUES);

    _ = petsc.MatAssemblyBegin(A, petsc.MAT_FINAL_ASSEMBLY);
    _ = petsc.MatAssemblyEnd(A, petsc.MAT_FINAL_ASSEMBLY);

    // Create RHS

    var b: petsc.Vec = undefined;
    _ = petsc.VecCreate(petsc.PETSC_COMM_WORLD, &b);
    _ = petsc.VecSetSizes(b, petsc.PETSC_DECIDE, m);
    _ = petsc.VecSetFromOptions(b);

    // Set RHS b = [1, 4]^T
    const rhs = [_]petsc.PetscScalar{ 1, 4 };
    _ = petsc.VecSetValues(b, 2, rows[0..], rhs[0..], petsc.INSERT_VALUES); // TODO: use better version

    _ = petsc.VecAssemblyBegin(b);
    _ = petsc.VecAssemblyEnd(b);

    // Solver

    var ksp: petsc.KSP = undefined;
    _ = petsc.KSPCreate(petsc.PETSC_COMM_WORLD, &ksp);
    _ = petsc.KSPSetOperators(ksp, A, A);

    _ = petsc.KSPSolve(ksp, b, b);

    var res: [*c]const petsc.PetscScalar = undefined;
    _ = petsc.VecGetArrayRead(b, &res);
    std.debug.print("{} {}\n", .{ res[0], res[1] });

    _ = petsc.KSPDestroy(&ksp);
    _ = petsc.VecDestroy(&b);
    _ = petsc.MatDestroy(&A);

    _ = petsc.PetscFinalize();
}
