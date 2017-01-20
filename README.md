# GraGLeS_3D
The brand new parallel level-set grain growth solver

Grain Growth Level Set 3dim (GraGLeS_3D) is a software kit to simulate the grain growth phenoma as it occurs in polycrystalline materials. The software is designed to account for anisotropic grain boundaries and finite, in particular, low triple junction mobilities. The algorithm utilizes an OpenMP parallelization strategy and is optimized (and ongoing work) for ccNUMA archictecture.

An installation guide is provided here:

http://gragles.readthedocs.org/en/latest/

Operating the 3D version is esentially the same as the 2D. You will need to install jemalloc on your system to overload the linux allocator.

For detailed information about the physics we refer to our Acta Materialia paper (private version attached to the git):

Articles:

A highly efficient 3D level-set grain growth algorithm tailored for ccNUMA architecture
C. Mießen, N. Velinov, G. Gottstein, L.A. Barrales-Mora
https://arxiv.org/submit/1781715/view
This article is submitted to IOP MSMSE.

An Advanced Level Set Approach to Grain Growth - Accounting for Grain Boundary Anisotropy and Finite Triple Junction Mobility, Acta Materialia, Volume 99, 15 October 2015, Pages 39–48
DOI information: 10.1016/j.actamat.2015.07.040.
bibtex:

@article{Miessen2015,
Author = {Mie{\ss}en, C. and Liesenjohann, M. and Barrales-Mora, L. A. and Shvindlerman, L. S. and Gottstein, G.},
Title = {An advanced level set approach to grain growth – Accounting for grain boundary anisotropy and finite triple junction mobility},
Journal = {Acta Materialia},
Volume = {99},
Pages = {39-48},
Year = {2015} 
}

Acknowledgements:

The authors gratefully acknowledge the support from the FZJuelich and the RWTH Aachen University for granting computing time within the frame of the JARAHPC project no 6687.
