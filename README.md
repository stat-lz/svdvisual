# SVD visualization matlab, python and R functions

(updated readme on Jan 30, 2024)

This is the site that holds the SVD visualization functions. See two major references related to this site. 

- Zhang, L, Marron, J. S., Shen, H., and Zhu, Z. (2007) Singular Value Decomposition and its Visualization, *Journal of Computational and Graphical Statistics*, 16(4), pp. 833-854.

- Zhang, L. (2007) Functional Singular Value Decomposition and Multi-Resolution Anomaly Detection, PhD dissertation, UNC-CH. 

A major function is

svd3dplot

related functions include

columnmean, rowmean, ad doublemean functions for calculating different types of means. 
svdls is a modified function to calculate a unique set of singular value decomposition
imagels is an inner function to show top-view image view of a matrix, that is associate with SVD

Please make sure MATLAB has svds function available. 

We have an early version of R package. You can download it and just install it into your R system manually.

Update on Jan 30. We have uploaded a python function svd3dplot.py in the system, this is an alpha version that can automatically handle 3 components. We will design a more flexible version later. 
