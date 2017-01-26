# GraGLeS_3D
The brand new parallel level-set grain growth solver - faster, more effitient and bigger networks than ever before

Grain Growth Level Set 3dim (GraGLeS_3D) is a software kit to simulate the grain growth and related phenoma such as recrystallization in polycrytalline materials. The software is designed to account for anisotropic grain boundaries, finite, in particular, low triple junction mobilities (ongoing work) and any additional driving force such as bulk energy densities resulting from anisotropic magnetic
susceptibilities or orientation dependent stored elastic energy densities.

The algorithm considers pysical units as input argumenst but utilizes normilized quantities internally. To mimic very specific material properties, such as grain size, GB mobilities, GB energies and bulk energies or even microstructure with certain orientation gradients or stored elastic energies related to orientations of grains or subgrains, you are encouraged to use the IMM microstructure generator:
https://github.com/GraGLeS/IMM_MicrostructureGenerator

The GraGLeS algorithm utilizes an OpenMP parallelization strategy and is optimized for ccNUMA archictecture. If you want to switch on this features set the <GrainScheduler> to 1. 0 is default.

An installation guide is provided here:
http://gragles.readthedocs.org/en/latest/

Operating the 3D version is essentially the same as the 2D. You will need to install jemalloc on your system to overload linux malloc to assure thread-local memory placement.

For detailed information about the algorithm we refer to our publications:

Articles:

A highly efficient 3D level-set grain growth algorithm tailored for ccNUMA architecture
C. Mießen, N. Velinov, G. Gottstein, L.A. Barrales-Mora
https://arxiv.org/abs/1701.06658
Cite as:	arXiv:1701.06658
This article is submitted to IOP MSMSE.

An Advanced Level Set Approach to Grain Growth - Accounting for Grain Boundary Anisotropy and Finite Triple Junction Mobility, Acta Materialia, Volume 99, 15 October 2015, Pages 39–48
DOI information: 10.1016/j.actamat.2015.07.040.

bibtex:
@article{
Miessen2015,
Author = {Mie{\ss}en, C. and Liesenjohann, M. and Barrales-Mora, L. A. and Shvindlerman, L. S. and Gottstein, G.},
Title = {An advanced level set approach to grain growth – Accounting for grain boundary anisotropy and finite triple junction mobility},
Journal = {Acta Materialia},
Volume = {99},
Pages = {39-48},
Year = {2015} 
}

(private version attached to the 2D git)



Acknowledgements:

The authors gratefully acknowledge the support from the FZJuelich and the RWTH Aachen University for granting computing time within the frame of the JARAHPC project no 6687.
