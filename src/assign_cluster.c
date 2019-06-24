/**************************************************************/
/*! \file assign_cluster.c
    \brief This file computes the algorithm sLSVD.
    (Learning the model, gu and participation)

    \author    Evangelia Christakopoulou
    \version   1.0
    \date      2017
    \copyright University of Minnesota, 2017
    \license   GPL v2
*/
/**************************************************************/

#include<slim.h>

int main(int argc, char *argv[])
{

  ctrl_t *ctrl;
  gk_csr_t *train;
  int my_id, my_num_procs;
  int i;
  int j;
  int d;
  int u;
  size_t lnlen;
  char *line=NULL;
  FILE *fpin=NULL;
  FILE *fpout=NULL;
  double value = 0;
  int cluster_id = 0;
  gk_csr_t* local_train;
  int user;
  int k;
  int min_assignment;
  double min_error;
  SVDRec result;
  int iter;
  int users_diff;
  int pos, pos2;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &my_num_procs);
  
  // parse command line 
  ctrl = create_ctrl();
  parse_cmdline(ctrl, argc, argv);

  ctrl->id = my_id;
  ctrl->num_procs = my_num_procs;

  // Reading input files
  train = gk_csr_Read(ctrl->train_file, GK_CSR_FMT_CSR, 1, 1);
  gk_csr_CreateIndex(train, GK_CSR_COL);

  fpin = gk_fopen(ctrl->dimensions_file, "r", "gk_readfile");

  int dimensions[ctrl->num_clusters];
  int total_dimensions = 0;
  cluster_id = 0;
  while(gk_getline(&line, &lnlen, fpin) != -1) {
    sscanf(line, "%le", &value);
    dimensions[cluster_id] = value;
    total_dimensions += dimensions[cluster_id];
    cluster_id++;
  }
  
  gk_fclose(fpin);

  int assignment[train->nrows];
  int prev_assignment[train->nrows];
  fpin = gk_fopen(ctrl->participation_file, "r", "gk_readfile");
  u = 0;
  while(gk_getline(&line, &lnlen, fpin) != -1) {
    sscanf(line, "%le", &value);
    assignment[u] = value;
    prev_assignment[u] = value; 
    u++;
  }
  gk_fclose(fpin);
  

  //Setup for MPI 
  double *proj_usv = gk_malloc(train->ncols*sizeof(double),"usv");
  double *error= gk_malloc(train->ncols*sizeof(double),"error");
  double *local_error = gk_malloc((ctrl->num_clusters)*sizeof(double),"error");
  int *sizes = NULL;
  int *displays = NULL;
  int *cluster_sizes = NULL;
  int *cluster_displays = NULL;
  int *vcluster_sizes = NULL;
  int *vcluster_displays = NULL;


  iter = 0;
  users_diff = train->nrows;

  int datasize = train->nrows;
  int step = (datasize / ctrl->num_procs) +
    (ctrl->id < (datasize % ctrl->num_procs) ? 1 : 0);
  int startu = ((datasize / ctrl->num_procs) * ctrl->id) +
    gk_min(ctrl->id, datasize % ctrl->num_procs);
  int endu = startu + step;
  if ((endu < datasize) && (ctrl->id == ctrl->num_procs - 1)) {
    endu = datasize;
    step = datasize - startu;
  }

  int clustersize = ctrl->num_clusters;
  int clusterstep = (clustersize / ctrl->num_procs) +
    (ctrl->id < (clustersize % ctrl->num_procs) ? 1 : 0);
  int startcluster = ((clustersize / ctrl->num_procs) * ctrl->id) +
    gk_min(ctrl->id, clustersize % ctrl->num_procs);
  int endcluster = startcluster + clusterstep;
  if ((endcluster < clustersize) && (ctrl->id == ctrl->num_procs - 1)) {
    endcluster = clustersize;
    clusterstep = clustersize - startcluster;
  }

  int *partial_participation = gk_imalloc(step, "malloc new g");
  if (ctrl->id == 0) {
    sizes = gk_ismalloc(ctrl->num_procs, 0, "malloc sizes");
    displays = gk_ismalloc(ctrl->num_procs, 0, "malloc displays");
    cluster_sizes = gk_ismalloc(ctrl->num_procs, 0, "malloc sizes");
    cluster_displays = gk_ismalloc(ctrl->num_procs, 0, "malloc displays");
    vcluster_sizes = gk_ismalloc(ctrl->num_procs, 0, "malloc sizes");
    vcluster_displays = gk_ismalloc(ctrl->num_procs, 0, "malloc displays");
  }
  MPI_Gather(&step, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (ctrl->id == 0) {
    displays[0] = 0;
    for (i = 0; i < ctrl->num_procs; i++) {
      if (i != ctrl->num_procs - 1)
        displays[i + 1] = displays[i] + sizes[i];
    }
  }
  
  int sclustersize = 0;
  for(i=startcluster; i<endcluster; i++){
    sclustersize += dimensions[i];
  }
  int vclustersize = sclustersize*train->ncols;
  double *partial_S = gk_malloc(sclustersize*sizeof(double),"gk_S");
  double *partial_V = gk_malloc(vclustersize*sizeof(double),"gk_V");
  double *partial_SV = gk_malloc(vclustersize*sizeof(double),"gk_SV");


  MPI_Gather(&sclustersize, 1, MPI_INT, cluster_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&vclustersize, 1, MPI_INT, vcluster_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (ctrl->id == 0) {
    cluster_displays[0] = 0;
    for (i = 0; i < ctrl->num_procs; i++) {
      if (i != ctrl->num_procs - 1)
        cluster_displays[i + 1] = cluster_displays[i] + cluster_sizes[i];
    }

    vcluster_displays[0] = 0;
    for (i = 0; i < ctrl->num_procs; i++) {
      if (i != ctrl->num_procs - 1)
        vcluster_displays[i + 1] = vcluster_displays[i] + vcluster_sizes[i];
    }

  }

  
  double *long_V = gk_malloc(total_dimensions*train->ncols*sizeof(double),"malloc long V");
  double *long_SV = gk_malloc(total_dimensions*train->ncols*sizeof(double),"malloc long V");
  double *long_S = gk_malloc(total_dimensions*sizeof(double),"malloc long S");
  double *long_projection = gk_malloc(total_dimensions*sizeof(double),"malloc long S");

  
  while(users_diff >= 0.01*(train->nrows)){
    if(ctrl->id==0)
      printf("iter %d\n",iter);
    
    int spos = 0;
    for(cluster_id=startcluster; cluster_id<endcluster; cluster_id++){
      int nrows=0;
      int nnz=0;
      
      for(i=0;i<train->nrows;i++){
	if(assignment[i]==cluster_id){
	  nrows++;
	  nnz+=train->rowptr[i+1]-train->rowptr[i];
	}
      }
      
      //Create the local training matrix for cluster_id
      local_train = gk_csr_Create();
      local_train->nrows = nrows;
      local_train->ncols = train->ncols;
      local_train->rowptr = gk_zmalloc(local_train->nrows+1, "gk_csr_ExtractPartition: rowptr");
      local_train->rowind = gk_imalloc(nnz, "gk_csr_ExtractPartition: rowind");
      local_train->rowval = gk_fmalloc(nnz, "gk_csr_ExtractPartition: rowval");
      
      nnz = 0;
      k = 0; 
      local_train->rowptr[0] = 0;
      for (i=0;i<train->nrows;i++){
	if(assignment[i]==cluster_id){
	  for(j=train->rowptr[i]; j<train->rowptr[i+1]; j++){
	    local_train->rowind[nnz] = train->rowind[j];
	    local_train->rowval[nnz] = train->rowval[j];
	    nnz++;
	  }
	  k++;
	  local_train->rowptr[k] = nnz;
	}
      }
      
      
      gk_csr_CreateIndex(local_train, GK_CSR_COL);          

      SMat svd_local_train = (SMat) calloc(1, sizeof(struct smat));
      svd_local_train->rows = local_train->nrows;
      svd_local_train->cols = train->ncols;
      svd_local_train->vals = nnz;
      
      svd_local_train->pointr = gk_malloc((svd_local_train->cols + 1)*sizeof(long), "svdNewSMat: pointr");
      svd_local_train->rowind = gk_malloc(nnz*sizeof(long), "svdNewSMat: rowind");
      svd_local_train->value  = gk_malloc(nnz*sizeof(double), "svdNewSMat: value");
      
      nnz = 0;
      svd_local_train->pointr[0] = 0;
      for (i=0;i<train->ncols;i++){
	for(j=local_train->colptr[i]; j<local_train->colptr[i+1]; j++){
	  svd_local_train->rowind[nnz] = (long)local_train->colind[j];
	  svd_local_train->value[nnz] = (double)local_train->colval[j];
	  nnz++;
	}
	svd_local_train->pointr[i+1] = nnz;
      }
      gk_csr_Free(&local_train);

      //Compute the SVD for the local training matrix for the cluster cluster_id
      result = svdLAS2A(svd_local_train, dimensions[cluster_id]);

      svdFreeSMat(svd_local_train);
      k= spos*train->ncols;
      for(d = 0; d < dimensions[cluster_id]; d++){
	for(j=0; j<train->ncols; j++){
	  partial_SV[k] =  result->Vt->value[d][j]*result->S[d];
	  k++;
	}
      }
      
      
      memcpy(partial_S + spos,result->S,dimensions[cluster_id]*sizeof(double));
      
      i = spos*train->ncols;
      for(d=0; d<dimensions[cluster_id]; d++){
	memcpy(partial_V + i,result->Vt->value[d],train->ncols*sizeof(double));
	i+= train->ncols;
      }
      
      
      svdFreeSVDRec(result);
      spos += dimensions[cluster_id];
      
    }

    MPI_Gatherv(partial_S, sclustersize, MPI_DOUBLE, long_S,
		cluster_sizes, cluster_displays, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Bcast(long_S, total_dimensions, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Gatherv(partial_V, vclustersize, MPI_DOUBLE, long_V,
		vcluster_sizes, vcluster_displays, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Bcast(long_V, total_dimensions*(train->ncols), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Gatherv(partial_SV, vclustersize, MPI_DOUBLE, long_SV,
		vcluster_sizes, vcluster_displays, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Bcast(long_SV, total_dimensions*(train->ncols), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //Assign every user to the cluster with the smallest error
    for(user= startu; user<endu; user++){
      pos = 0;
      pos2 = 0;
      for(cluster_id=0; cluster_id<ctrl->num_clusters; cluster_id++){
	local_error[cluster_id] = 0;
	if(cluster_id > 0){
	  pos += dimensions[cluster_id-1]*(train->ncols);
	  pos2 += dimensions[cluster_id-1];
	}
	for(d = 0; d < dimensions[cluster_id]; d++){
	  long_projection[pos2 + d] = 0;
	  for(j=train->rowptr[user]; j<train->rowptr[user+1]; j++){
	    k = train->rowind[j];
	    long_projection[pos2 + d] +=  long_V[pos + d*train->ncols + k]*(1/long_S[pos2 + d]);  
	  }
	}
	
	for(i=0;i<train->ncols;i++){
	  proj_usv[i] = 0;
	  for(d=0; d<dimensions[cluster_id]; d++){
	    proj_usv[i] += long_projection[pos2 + d]*long_SV[pos + d*train->ncols + i];
	  }
	}
	
	
	for(i=0; i<train->ncols; i++){
	  error[i]=proj_usv[i];
	}
	for(i=train->rowptr[user];i<train->rowptr[user+1];i++){
	  j=train->rowind[i];
	  error[j] = proj_usv[j] - 1;
	}
	for(i=0; i<train->ncols; i++)
	  local_error[cluster_id] += error[i]*error[i];
      }//end loop of clusters

      min_assignment = 0; 
      min_error = local_error[0];
      
      for(cluster_id=1; cluster_id <ctrl->num_clusters; cluster_id++){
	if(local_error[cluster_id]<min_error){
	  min_error = local_error[cluster_id];
	  min_assignment = cluster_id;
	}
      }
      partial_participation[user-startu] = min_assignment;
    }

    MPI_Gatherv(partial_participation, step, MPI_INT, assignment,
		sizes, displays, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(assignment, train->nrows, MPI_INT, 0, MPI_COMM_WORLD);

    users_diff = 0;
    for(user=0; user<train->nrows; user++){
      if(assignment[user]!=prev_assignment[user])
	users_diff++;
    }
    if(ctrl->id==0)
      printf("Users jumping clusters are %d\n",users_diff);
    
  
    for(user=0;user<train->nrows; user++){
      prev_assignment[user] = assignment[user];
    }
    
    if(users_diff < 0.01*(train->nrows)){
      fpout = gk_fopen("new_assignment", "w+","gk_writefile");
      for(user=0; user<train->nrows; user++){
	fprintf(fpout, "%d\n",assignment[user]);
      }
      gk_fclose(fpout);
    }
    iter++;
  }

     
  gk_free((void **)&long_S, &long_V, &long_projection, &long_SV, LTERM);
  gk_free((void **)&sizes, &displays, &cluster_sizes, &cluster_displays, &vcluster_sizes, &vcluster_displays, LTERM);
  gk_free((void **)&partial_S, &partial_V, &partial_SV, LTERM);
  gk_free((void **)&partial_participation, LTERM);
  gk_free((void **)&local_error, &error, &proj_usv, LTERM);
  gk_csr_Free(&train);
  free_ctrl(ctrl);
  MPI_Finalize();
}
