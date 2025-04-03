# Stoat_cxx

Stoat_cxx is a C++ implementation (64 bits architecture) of the Stoat tool, designed specifically for advanced Genome-Wide Association Studies (GWAS) that focus on snarl structures within pangenome graphs. This C++ version brings significant performance benefits—such as faster execution and lower memory usage—compared to its Python counterpart, making it an efficient tool for handling large-scale pangenome data.

## Dependency

Manual installation : 

- jansson 
- libbdsg
- htslib
- eigen3
- boost
- Catch2 v3

Or you can use the Dockerfile provided [here](Dockerfile)

## Building

```bash
git clone --branch stoat_cxx https://github.com/Plogeur/STOAT.git
cd stoat_cxx

mkdir build && cd build
cmake .. && make -j 4

# ./unit_tests
# ./stoat
```

## Benchmarking

Performance comparison between the **C++ Stoat (stoat_cxx)** and the **Python3 Stoat (stoat)**:
- **Matrix Construction**: The C++ implementation is approximately **135% faster** in building the matrix.
- **Snarl p-value Calculation**: Achieves a **550% speed improvement** over the Python version, significantly reducing computation time.
- **Memory Efficiency**: The C++ version consumes **150% less memory**, making it a more efficient choice for large datasets.
