GLSVD: Global and Local Singular Value Decomposition
________________________________________________________________________
________________________________________________________________________

This is a toolkit of MPI software for top-N recommendation which implements 
the methods: rGLSVD, sGLSVD, rLSVD and sLSVD presented in the paper: 
Local Latent Space Models for Top-N Recommendation.

Building & Installing
----------------------
The code has been tested on Linux (x86_64) and it is MPI-dependent.
In order to compile the code, gcc (https://gcc.gnu.org/), 
CMake (https://cmake.org/), 
and MPI (https://www.open-mpi.org/) are required.

The code can run on MacOS as well, if the compiler used is gcc.
Also, the Command Line Tools will be of use: xcode-select --install

In order to compile the code, please run the command:

       $ ./build.sh 

Then, the executables rlvsvd, slsvd, rglsvd, sglsvd and gpredict are created 
in glsvd/build/src.

How to Run GLSVD Methods
-------------------------

rLSVD
------

	$ mpirun -np 4 ./rlsvd -train_file=train_1 -test_file=test_1 
	  -dimensions_file=dimensions -participation_file=assignment 
	  -num_clusters=10 [-topn=20]

* This is an example command on how to run rLSVD with 4 cores for the 
  train_file 'train_1' and the test_file 'test_1', 
  with the local rank of each cluster specified in dimensions_file 'dimensions',
  and with the initial assignment of users to clusters in 
  participation_file 'assignment' for 10 clusters. 
  The size of the recommendation list is set to be 20 
  (instead of the default which is 10).
*
  This command performs both the training and the testing for rLSVD.
 
sLSVD
------

        $ mpirun -np 4 ./slsvd -train_file=train_1
          -dimensions_file=dimensions -participation_file=assignment
          -num_clusters=10 

* This is an example command on how to train sLSVD with 4 cores for the
  train_file 'train_1',
  with the local rank of each cluster specified in dimensions_file 'dimensions',
  and with the initial assignment of users to clusters in
  participation_file 'assignment' for 10 clusters.

*  Running the above command will output a file named new_${participation_file},
    which in our example would be 'new_assignment', that 
   contains the new assignments of users to clusters.

* In order to see how this method did in the test set, we run the following:


       $ mpirun -np 4 ./rlsvd -train_file=train_1 -test_file=test_1
          -dimensions_file=dimensions -participation_file=new_assignment
          -num_clusters=10 [-topn=20]

rGLSVD
------

        $ mpirun -np 4 ./rglsvd -train_file=train_1
          -dimensions_file=dimensions -num_dimensions=10
	  -participation_file=assignment -gu_file=gu_file 
	  -num_clusters=10 

* This is an example command on how to train rGLSVD with 4 cores 
  for the  train_file 'train_1',
  with the local rank of each cluster specified in dimensions_file 'dimensions',
  with the global rank being set to 10, 
  with the initial assignment of users to clusters in
  participation_file 'assignment' for 10 clusters, 
  and with the initial personalized weight of each user in file 'gu_file'.

*  Running the above command will output a file named new_${gu_file},
    which in our example would be 'new_gu_file', that
   contains the new personalized weights of users.

* In order to see how this method did in the test set, we run the following:

        $ mpirun -np 4 ./gpredict -train_file=train_1 -test_file=test_1
          -dimensions_file=dimensions -num_dimensions=10
          -participation_file=assignment -gu_file=new_gu_file
          -num_clusters=10 [-topn=20]


sGLSVD
------

        $ mpirun -np 4 ./sglsvd -train_file=train_1
          -dimensions_file=dimensions -num_dimensions=10
          -participation_file=assignment -gu_file=gu_file
          -num_clusters=10

* This is an example command on how to train sGLSVD with 4 cores
  for the  train_file 'train_1',
  with the local rank of each cluster specified in dimensions_file 'dimensions',
  with the global rank being set to 10,
  with the initial assignment of users to clusters in
  participation_file 'assignment' for 10 clusters,
  and with the initial personalized weight of each user in file 'gu_file'.

*  Running the above command will output a file named new_${gu_file},
    which in our example would be 'new_gu_file', that
   contains the new personalized weights of users.
   It will also output a file named new_${participation_file}, which 
   in our example would be 'new_assignment', that contains the new 
   assignment of users to clusters.

* In order to see how this method did in the test set, we run the following:

        $ mpirun -np 4 ./gpredict -train_file=train_1 -test_file=test_1
          -dimensions_file=dimensions -num_dimensions=10
          -participation_file=new_assignment -gu_file=new_gu_file
          -num_clusters=10 [-topn=20]


To see a complete list of the available parameters, along with a small 
description, please open the manpage for rLSVD - sLSVD - rGLSVD - sGLSVD.

         $ ./rlsvd -help
	 $ ./slsvd -help
	 $ ./rglsvd -help
	 $ ./sglsvd -help

File Formats
-------------
The train, and test files are stored in CSR
(https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_.28CSR.2C_CRS_or_Yale_format.29) format.

We present the file formats for a toy example of 5 users and 3 items, using 2 
clusters, with the algorithm converging in 3 iterations.

Also please look at the example files provided in the folder 'example_files'. 
The example given is using 10 clusters.

* Train File
  
  The train file has a separate row per user. Every row has #item_id #rating 
  #item_id #rating e.t.c. Since we are using implicit data, all our ratings are 1.
  The example train file below shows that the first user purchased items 2 and 3,
  the second user purchased item 3, the third user purchased item 1, the fourth 
  user puchased items 1 and 2, and finally the fifth user purchased items 1 and 3.
  
  2 1 3 1
  3 1
  1 1
  1 1 2 1
  1 1 3 1

* Test File
  
  The test file has a separate row per user. Every row has #item_id #rating, 
  specifying which items have been hidden per user.
  The ratings are '1' here, too. 
  The example test file shown below shows that for the first user, the first 
  item was put in the test set, for the second, the third and the fifth users,
  the second item was put in the test set, while for the fourth user the third 
  item was hidden.

  1 1
  2 1
  2 1
  3 1
  2 1

* Participation File
  
  The participation file has a separate row per user.
  Every row contains the id of the cluster the user is assigned to: 
  {0,1,...,#num_clusters-1}. The number of clusters presented in this file should 
  be the same as the number of clusters specified by the parameter num_clusters.
  The example participation file for 5 users is shown below:

  0
  0
  1
  1
  0

* Gu File
  
  The gu file has a separate row per user, showing the personalized weight of the 
  user, which takes values between 0 and 1. The example gu file for 5 users 
  is shown below:
  
  0.3
  0.7
  0
  1
  0

* Dimensions file

  The dimensions file has a separate row per cluster, containing the local rank 
  that will be used for that cluster for the truncated SVD. An example file for 
  3 clusters would be:
  
  10
  5
  20

Acknowledgements
-----------------
We would like to thank Dr Douglas Rohde (http://tedlab.mit.edu/~dr/SVDLIBC/) for 
providing the SVDLIBC library.

We would also like to thank Professor George Karypis (http://glaros.dtc.umn.edu/)
for providing the GKlib library.

License
-------
Please look at the file LICENSE.

Citation
---------
If you use our code in your research, please cite us:

@inproceedings{christakopoulou2018local,
  title={Local latent space models for top-n recommendation},
  author={Christakopoulou, Evangelia and Karypis, George},
  booktitle={Proceedings of the 24th ACM SIGKDD International Conference on Knowledge Discovery \& Data Mining},
  pages={1235--1243},
  year={2018},
  organization={ACM}
}

Contact Information
-------------------
GLSVD was written by Evangelia Christakopoulou and George Karypis.
If you have any issues, please send us an email to: evangel@cs.umn.edu, 
karypis@cs.umn.edu. 
