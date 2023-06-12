## BioMiCo

This is a repository for a forked version of our code that implements **a novel hierarchical model for Bayesian inference of microbial communities (BioMiCo)**. The model takes abundance data derived from environmental DNA, and models the composition of each sample by a two-level hierarchy of mixture distributions constrained by Dirichlet priors. BioMiCo is supervised, using known features for samples and appropriate prior constraints to overcome the challenges posed by many variables, sparse data, and large numbers of rare species. The model is trained on a portion of the data, where it learns how assemblages of species are mixed to form communities and how assemblages are related to the known features of each sample. Training yields a model that can predict the features of new samples.

The paper describing the method is:

>Shafiei, M., Dunn, K. A., Boon, E., MacDonald, S. M., Walsh, D. A., Gu, H., & Bielawski, J. P. (2015). BioMiCo: a supervised Bayesian model for inference of microbial community structure. Microbiome, 3(1), 1-15.

Please cite this paper if you use this repository or the orginal (below).

The orginal code was written by Mahdi Shafiei (Bielawski Group postdoc). The repository for the original version is here:  https://sourceforge.net/projects/biomico/
