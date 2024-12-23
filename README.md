# snarl_project_cxx

The snarl_project_cxx is a C++ implementation of the Snarl_project tool, designed specifically for advanced Genome-Wide Association Studies (GWAS) that focus on snarl structures within pangenome graphs. Unlike traditional GWAS tools, which typically analyze linear genome variants, Snarl_project is built to handle the complexity of pangenome graphs by extracting and analyzing "snarl" regions from VCF (Variant Call Format) files. These snarls represent complex structural variations, capturing nested and overlapping variant patterns across genomes.

By focusing on these snarl structures, snarl_project_cxx allows researchers to investigate genetic variations with a higher level of detail, particularly in the context of diverse populations and complex traits. This C++ version brings significant performance benefits—such as faster execution and lower memory usage—compared to its Python counterpart, making it an efficient tool for handling large-scale pangenome data.

## Benchmarking

Performance comparison between the **C++ snarl_project** and the **Python3 snarl_project**:
- **Matrix Construction**: The C++ implementation is approximately **135% faster** in building the matrix.
- **Snarl p-value Calculation**: Achieves a **550% speed improvement** over the Python version, significantly reducing computation time.
- **Memory Efficiency**: The C++ version consumes **150% less memory**, making it a more efficient choice for large datasets.

