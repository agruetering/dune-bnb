# Instances Overview
------------------
The instances of the numerical experiments in [BGM] can be found in this folder. 
The subfolders correspond to the following experiments: 
  - alpha_beta: Influence of the Tikhonov parameter and the penalty term of the box
constraints on the branch-and-bound algorithm [Table 1, BGM],
  - obj_red: Impact of the balance between branching and cutting plane iterations on the
branch-and-bound algorithm [Table 2, BGM],
  - tolerance: Influence of the relative allowed deviation (TOL) from the optimum on the
branch-and-bound algorithm  [Table 3, BGM], and
  - s\*max*: Performance of the branch-and-bound algorithm for instances generated with θ
switching points, allowing at most σ switching [Table 4, BGM].

The branch-and-bound tree in [Figure 2, BGM] is based on the instance s3max1/run4. 

### Bibliography
[BGM] C. Buchheim, A. Grütering, and C. Meyer, Parabolic optimal control problems with combinatorial switching constraints – Part III: Branch-and-bound 
    algorithm, arXiV preprint arXiv:2401.10018 (2024)
