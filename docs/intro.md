# Motivation

Many projects about CFD.
Here we want to focus on creating meshes for CFD, which is not so well covered.

We also want to concentrate on the implementation details.
Programming language: Zig

Previous implementation in Rust.

We can go into details why Rust experiance is not optimal (personal view):

- feels to high level
- I want to have full control over memory, b.c. I know what I am doing.
- We are not talking about security critical systems here - this is why this perhaps feels like this.
- We want to focus on how to efficiently utilise the hardware. I am not saying that the implementation
is fully optimized, but I tried to have good memory access patterns as a default.
- Motivate to look more into how to exploit computer hardware (as opposed to this trend to work more on higher level w/o
any understandind of the lower levels)
- C++ is just bad. Long compile times, hard syntax,...
- Other option C (or very trimmed down version of C++), but I like to try new things and Zig feels like a very nice language
- Very important (also compared to RUST): very nice C interopt!!!


# Structure

- Mathematical basics
- Implementation details
- Single Block smoothing
- Multi block smoothing (global)


# APIs we choose

For now we will rely on a compiler, the OS and a sparse matrix solving library.

We will not be looking into designing a GUI block structured mesher.

We can still reflect on how to do this.
