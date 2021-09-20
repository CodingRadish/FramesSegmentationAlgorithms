# FramesSegmentationAlgorithms
Repository contains files with two segmentation algorithms developed at the Warsaw University of Technology as part of an engineering thesis. The files are use as a plugin in the Frames software for point cloud processing.

Frames is a software which uses custom data structures and functions as well as the Eigen library.

Two small custom functions were added (these operations were more frequently used).

The cpp file contains the two segmentation algorithms and a parameter calculating algorithm. They are executed by a data structure called a "method" which is used to export the algorithm as a plugin in the Frames UI.
The three methods contain:
1. Parameter calculation: local plane coefficient
2. Segmentation by region growing
3. Segmentation using a normal vector histogram.

More information can be found in the engineering thesis entitled "Development of algorithm for cloud of point's segmentation based on orientation similarity".
The thesis is available at the WUT.