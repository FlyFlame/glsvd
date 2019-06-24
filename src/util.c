/**************************************************************/
/*! \file 
  
    \brief This file contains all the utility routines. 
*/
/**************************************************************/

#include<slim.h>

/**************************************************************/
/*! 
  \brief Create a ctrl structure wich contains all the default
         parameters
 
  \return ctrl_t* A pointer to a created ctrl structure
*/
/**************************************************************/
ctrl_t *create_ctrl()
{

  ctrl_t *ctrl = gk_malloc(sizeof(ctrl_t), "malloc ctrl");

  ctrl->train_file = NULL;
  ctrl->test_file = NULL;
  ctrl->participation_file = NULL;
  ctrl->gu_file = NULL;
  ctrl->dimensions_file = NULL;
  ctrl->num_clusters = 5;
  ctrl->num_dimensions = 10;
  ctrl->topn = 10;
  ctrl->num_procs = 1;
  ctrl->id = 0;

  return ctrl;

}


/**************************************************************/
/*! 
  \brief Free a ctrl structure
  
  \param[in] ctrl A pointer to a ctrl structure to be freed
*/
/**************************************************************/
void free_ctrl(ctrl_t * ctrl)
{

  gk_free((void **) &ctrl->train_file, LTERM);
  gk_free((void **) &ctrl->dimensions_file, LTERM);
  gk_free((void **) &ctrl->test_file, LTERM);
  gk_free((void **) &ctrl->gu_file, LTERM);
  gk_free((void **) &ctrl->participation_file, LTERM);

  gk_free((void **) &ctrl, LTERM);

}


/*****************************************************************************/
