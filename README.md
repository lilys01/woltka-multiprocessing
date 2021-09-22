# woltka-multiprocessing
Simple multiprocessing implemented in the classification step in [Woltka](https://github.com/qiyunzhu/woltka). When classifying on a directory of files, there is over a 2x speedup using 4 cores and over a 5x speedup using 12 cores. Multiprocessing was done using the [Pathos](https://github.com/uqfoundation/pathos) multiprocessing library which utilizies dill for synchronization. Only the classify function in workflow.py is uploaded here. 

## References: 

M.M. McKerns, L. Strand, T. Sullivan, A. Fang, M.A.G. Aivazis,
"Building a framework for predictive science", Proceedings of
the 10th Python in Science Conference, 2011;
http://arxiv.org/pdf/1202.1056

Michael McKerns and Michael Aivazis,
"pathos: a framework for heterogeneous computing", 2010- ;
https://uqfoundation.github.io/project/pathos

Zhu Q, Huang S, Gonzalez A, McGrath I, McDonald D, Haiminen N, et al. OGUs enable effective, phylogeny-aware analysis of even shallow metagenome community structures. bioRxiv. 2021. doi: https://doi.org/10.1101/2021.04.04.438427.
