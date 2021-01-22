# RefiNA: fast, accurate, and robust <em>post hoc</em> refinement of any network alignment solution
This repository contains a Python implementation of RefiNA, an unsupervised method for refining initial alignments between the nodes of two graphs.  RefiNA's goal is to increase the matched neighborhood consistency of the alignment solution produced by any of the many existing methods for network alignment, thereby improving its quality. 

**Paper**: Mark Heimann, Xiyuan Chen, Fatemeh Vahedian, and Danai Koutra. <a href="https://gemslab.github.io/papers/heimann-2021-RefiNA.pdf">"Refining Network Alignment to Improve Matched Neighborhood Consistency."</a> SIAM International Conference on Data Mining (SDM), 2021.

Included in this paper version is the full supplementary material accompanying the original submission, although SDM only officially publishes the main body of the paper.  

**Versioning and Dependencies**
Our implementation has been tested with Python 3.5.2.  The main technical stack includes NumPy, SciPy, scikit-learn, and NetworkX.  See requirements.txt for exact versions of all installed packages in our tested environment. 

**Example**
We include initial solutions obtained from base network alignment methods REGAL and NetAlign on two datasets used in our paper.  Both datasets are given as edgelists of a combined network containing the adjacency matrices of both networks to align as disconnected components.  The first dataset is the Arenas email dataset, where a simulated alignment scenario has been created by aligning one graph to a permuted, noisy copy of itself: both graphs have the same number of nodes, but the second graph has 5 percent noise introduced in the form of edge removal. The second dataset consists of the PPI networks of two different species: the first 1124 nodes in the combined network belong to the first network, and the remainder belong to the second network.  The first comes with a dictionary of true alignments between the node IDs of the two graphs, whereas the second graph has no underlying node correspondence.  See our paper for more details about the datasets and the base network alignment methods.

**Usage**
python main.py   (Performs dense refinement of REGAL's solution of the Arenas dataset with 5 percent noise.  See main.py for a list of flags that can be used to perform sparse refinement, change the dataset, or change the initial alignment solution.)

Please consider citing this paper if you find the code helpful. 
```bibtex
@inproceedings{refina,
 title={Refining Network Alignment to Improve Matched Neighborhood Consistency},
 author={Heimann, Mark and Chen, Xiyuan and Vahedian, Fatemeh and Koutra, Danai},
 booktitle={SIAM International Conference on Data Mining},
 organization={SIAM},
 year={2021},
}
```





