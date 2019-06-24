/**************************************************************/
/*! \file
  
    \brief This file contains all prototypes.
*/
/**************************************************************/


#ifndef __PROTO_H__
#define __PROTO_H__

/* cmd.c */
void parse_cmdline(ctrl_t * ctrl, int argc, char *argv[]);

/* util.c */
ctrl_t *create_ctrl();
void free_ctrl(ctrl_t * ctrl);

#endif
