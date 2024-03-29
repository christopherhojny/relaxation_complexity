This repository contains code used in the articles

   Computational Aspects of Relaxation Complexity: Possibilities and Limitations
   by Gennadiy Averkov, Christopher Hojny, and Matthias Schymura,

   Efficient MIP Techniques for Computing the Relaxation Complexity
   by Gennadiy Averkov, Christopher Hojny, and Matthias Schymura,

   and

   The role of rationality in integer-programming relaxations
   by Manuel Aprile, Gennadiy Averkov, Marco Di Summa, and Christopher Hojny.


The repository consists of three parts, each part has its own README file
with more detailed information on how to use the corresponding code.

Part 1 (computer-aided-proofs)
------------------------------

Part 1 contains scripts to execute the computer-based parts of the proof of
Lemma 8 in the article "Computational Aspects of Relaxation Complexity:
Possibilities and Limitations" as well as of the proof of Proposition 1 in
the article "The role of rationality in integer-programming relaxations".

Part 2 (sage-code)
------------------

Part 2 contains a Sage/SCIP implementation for computing observers of a
finite lattice-convex set and, based on the computed observers, the relaxation
complexity of finite lattice-convex sets (provided the set of observers is
finite).

Part 3 (scip-code)
------------------

Part 3 contains a C/C++ implementation to compute rc(X, Y) for a finite
full-dimensional lattice-convex set X and a set Y of integer points not
contained in X.
