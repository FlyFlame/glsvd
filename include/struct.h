/**************************************************************/
/*! \file
  
    \brief This file contains all the necessary data structures
*/
/**************************************************************/


#ifndef __STRUCT_H__
#define __STRUCT_H__



/**************************************************************/
/*!   
    \brief A data structure for ctrl parameters
 */
/**************************************************************/
typedef struct {

  /*! a file name that contains the training data in csr format */
  char *train_file;
  /*! a file name that contains the testing data in csr format */
  char *test_file;
  /*a file name containing for each cluster the local rank*/
  char *dimensions_file;
  /*! a file name containing the clustering assignment for each user*/
  char *participation_file;
  /* a file name containing the weight for each user */
  char *gu_file;
  /*! the number of recommendations to be recommended */
  int topn;
  /*! the number of clusters */
  int num_clusters;
  /*! the number of dimensions for global train */
  int num_dimensions;
  /*! the number of processors */
  int num_procs;
  /*! the id of the processor */
  int id;

} ctrl_t;


#endif
