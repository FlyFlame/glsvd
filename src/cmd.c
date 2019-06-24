/**************************************************************/
/*! \file
  
    \brief This file contains all the routines for parameter
           setup from the user
*/
/**************************************************************/

#include<slim.h>

/**************************************************************/
/*!
  \brief A structure for command-line options
 */
/**************************************************************/
static struct gk_option glsvd_options[] = {
  {"train_file", 1, 0, CMD_TRAIN_FILE},
  {"dimensions_file", 1, 0, CMD_DIMENSIONS_FILE},
  {"test_file", 1, 0, CMD_TEST_FILE},
  {"participation_file", 1, 0, CMD_PARTICIPATION_FILE},
  {"gu_file", 1, 0, CMD_GU_FILE},
  {"topn", 1, 0, CMD_TOPN},
  {"num_clusters", 1, 0, CMD_NUM_CLUSTERS},
  {"num_dimensions", 1, 0, CMD_NUM_DIMENSIONS},
  {"num_procs", 1, 0, CMD_NUM_PROCS},
  {"id", 1, 0, CMD_ID},
  {"help", 0, 0, CMD_HELP},
  {0, 0, 0, 0}
};



/**************************************************************/
/*! \brief Mini help
 */
/**************************************************************/
static char helpstr[][512] = {
  " 	 -train_file=string",
  " 		Specifies the input file which contains the training data. ",
  " 		This file should be in .csr format. ",
  " 		",
  " 	 -test_file=string",
  " 		Specifies the input file which contains the testing data. ",
  " 		This file should be in .csr format.",
  " 	",
  "      -participation_file=string",
  "              Specifies the file which contrains the clustering assignment.",
  "              The file has for each user the cluster they belong to (from 0 to #clusters-1).",
  " ",
  "      -gu_file=string",
  "              Specifies the file which contrains the user weights.",
  "              The file is an array of doubles from 0 to 1. ",
  " ",
  "        -topn=int",
  "               Specifies the number of recommendations to be produced for ",
  "               each user. The default value is 10.",
  " ",
  "        -num_clusters=int",
  "               Specifies the number of clusters. ",
  "               The default value is 5.",
  " ",
  "       -num_dimensions=int",
  "              Specifies the number of global dimensions for global svd. ",
  "              The default value is 10.",
  " ",
  "        -help",
  "               Print this message.",
  " 			",
  ""
};


/**************************************************************/
/*! \brief Entry point of the command-line argument parsing
  
    \param[out] ctrl  A ctrl structure to be filled out
    \param[in]  argc  Number of arguments
    \param[in]  argv  A list of arguments
*/
/**************************************************************/
void parse_cmdline(ctrl_t * ctrl, int argc, char *argv[])
{

  int c = -1, option_index = -1;

  if (ctrl == NULL)
    ctrl = create_ctrl();

  while ((c =
	  gk_getopt_long_only(argc, argv, "", glsvd_options,
			      &option_index)) != -1) {
    switch (c) {

    case CMD_TRAIN_FILE:
      ctrl->train_file = gk_strdup(gk_optarg);
      break;

    case CMD_DIMENSIONS_FILE:
      ctrl->dimensions_file = gk_strdup(gk_optarg);
      break;
      
    case CMD_TEST_FILE:
      ctrl->test_file = gk_strdup(gk_optarg);
      break;

    case CMD_PARTICIPATION_FILE:
      ctrl->participation_file = gk_strdup(gk_optarg);
      break;

    case CMD_GU_FILE:
      ctrl->gu_file = gk_strdup(gk_optarg);
      break;

    case CMD_TOPN:
      ctrl->topn = atoi(gk_optarg);
      break;

    case CMD_NUM_CLUSTERS:
      ctrl->num_clusters = atoi(gk_optarg);
      break;

    case CMD_NUM_DIMENSIONS:
      ctrl->num_dimensions = atoi(gk_optarg);
      break;
      
    case CMD_NUM_PROCS:
      ctrl->num_procs = atoi(gk_optarg);
      break;

    case CMD_ID:
      ctrl->id = atoi(gk_optarg);
      break;

    case CMD_HELP:
      printf("%s\n %s %s %s \n %s \n", " ", "   Usage: " ,argv[0], "[options] ","");
      for (int i = 0; strlen(helpstr[i]) > 0; i++)
	printf("%s\n", helpstr[i]);
      exit(0);

    case '?':
    default:
      printf("Illegal command-line option(s) %s\n", gk_optarg);
      exit(0);

    }
  }

  if (argc - gk_optind != 0 || argc == 1) {
    printf("%s\n %s %s %s \n %s %s %s \n %s \n", " ", "   Usage: " ,argv[0], "[options] ", "          use '" ,argv[0], " -help' for a summary of the options.","");
    exit(0);
  }


}
