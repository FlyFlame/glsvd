/**************************************************************/
/*! \file calculate_gu.c
    \brief This file implements the algorithm rGLSVD.

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
  gk_csr_t *train, *orig_train;
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
  char new_gu_name[200]="";
  int k;
  SVDRec result;
  SVDRec global_result;
  int iter;
  int global_nnz;
  int users_diff;
  double nom, denom;
  int pos;
  
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

  orig_train = gk_csr_Dup(train);
  gk_csr_CreateIndex(orig_train, GK_CSR_COL);

  fpin = gk_fopen(ctrl->dimensions_file, "r", "gk_readfile");

  int dimensions[ctrl->num_clusters];
  int assignment[train->nrows];
  double gu[train->nrows];
  double new_gu[train->nrows];

  if(ctrl->id==0){
    cluster_id = 0;
    while(gk_getline(&line, &lnlen, fpin) != -1) {
      sscanf(line, "%le", &value);
      dimensions[cluster_id] = value;
      cluster_id++;
    }
    
    gk_fclose(fpin);
    
    fpin = gk_fopen(ctrl->participation_file, "r", "gk_readfile");
    u = 0;
    while(gk_getline(&line, &lnlen, fpin) != -1) {
      sscanf(line, "%le", &value);
      assignment[u] = value;
      u++;
    }
    gk_fclose(fpin);
    
    fpin = gk_fopen(ctrl->gu_file, "r", "gk_readfile gu");
    u = 0;
    while(gk_getline(&line, &lnlen, fpin) != -1) {
      sscanf(line, "%le", &value);
      gu[u] = value;
      new_gu[u] = value;
      u++;
    }
    gk_fclose(fpin);
  }

  //Setup for MPI  
  MPI_Bcast(dimensions, ctrl->num_clusters, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(assignment, train->nrows, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(gu, train->nrows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(new_gu, train->nrows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  int total_dimensions = 0;
  for(i=0;i<ctrl->num_clusters;i++){
    total_dimensions += dimensions[i];
  }
  
  iter = 0;
  users_diff = train->nrows;

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

  int *cluster_sizes = NULL;
  int *cluster_displays = NULL;
  int *vcluster_sizes = NULL;
  int *vcluster_displays = NULL;

  if (ctrl->id == 0) {
    cluster_sizes = gk_ismalloc(ctrl->num_procs, 0, "malloc sizes");
    cluster_displays = gk_ismalloc(ctrl->num_procs, 0, "malloc displays");
    vcluster_sizes = gk_ismalloc(ctrl->num_procs, 0, "malloc sizes");
    vcluster_displays = gk_ismalloc(ctrl->num_procs, 0, "malloc displays");
  }

  int *rows = gk_malloc(ctrl->num_clusters*sizeof(int),"malloc rows");

  for(i=0; i<ctrl->num_clusters; i++){
    rows[i] = 0;
    for(u=0; u<train->nrows; u++){
      if(assignment[u]==i)
        rows[i] ++;
    }
  }
  
  int sclustersize = 0;
  int usclustersize = 0;
  for(i=startcluster; i<endcluster; i++){
    sclustersize += dimensions[i];
    usclustersize += rows[i]*dimensions[i];
  }
  
  int vclustersize = sclustersize*train->ncols;
  double *partial_US = gk_malloc(usclustersize*sizeof(double),"gk_S");
  double *partial_V = gk_malloc(vclustersize*sizeof(double),"gk_V");


  MPI_Gather(&usclustersize, 1, MPI_INT, cluster_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&vclustersize, 1, MPI_INT, vcluster_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int us_total_size;
  if (ctrl->id == 0) {
    us_total_size = 0;
    cluster_displays[0] = 0;
    for (i = 0; i < ctrl->num_procs; i++) {
      us_total_size += cluster_sizes[i];
      if (i != ctrl->num_procs - 1)
        cluster_displays[i + 1] = cluster_displays[i] + cluster_sizes[i];
    }

    vcluster_displays[0] = 0;
    for (i = 0; i < ctrl->num_procs; i++) {
      if (i != ctrl->num_procs - 1)
        vcluster_displays[i + 1] = vcluster_displays[i] + vcluster_sizes[i];
    }
    
  }

  MPI_Bcast(&us_total_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  double *long_V = gk_malloc(total_dimensions*train->ncols*sizeof(double),"malloc long V");
  double *long_US = gk_malloc(us_total_size*sizeof(double),"malloc long US");

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

  int *local_positions = gk_malloc(ctrl->num_clusters*sizeof(int),"malloc local positions");
  int *local_displays = gk_malloc(ctrl->num_clusters*sizeof(int),"malloc local displays");
  int *local_corr = gk_malloc(train->nrows*sizeof(int),"malloc local correspondence");
  
  local_displays[0] = 0;

  for(cluster_id=0; cluster_id < ctrl->num_clusters; cluster_id++){
    local_positions[cluster_id] = 0;
    if(cluster_id>0){
      local_displays[cluster_id] = local_displays[cluster_id-1]+
        rows[cluster_id-1]*dimensions[cluster_id-1];
    }
  }
  for(u=0; u < train->nrows; u++){
    for(cluster_id=0; cluster_id<ctrl->num_clusters; cluster_id++){
      if(assignment[u] == cluster_id){
	local_corr[u] = local_displays[cluster_id] + local_positions[cluster_id];
	local_positions[cluster_id]+= dimensions[cluster_id];
      }
    }
  }
  
  gk_free((void **)&rows,LTERM);
    
  double *partial_gu = gk_malloc(step*sizeof(double),"malloc gu");
  double error;

  double *total_error = NULL;
  int *sizes = NULL;
  int *displays = NULL;
  int *local_indices = gk_malloc(orig_train->ncols*sizeof(int),"malloc index");
  double *global = gk_malloc(orig_train->ncols*sizeof(double),"malloc global");
  double *local = gk_malloc(orig_train->ncols*sizeof(double),"malloc local");
  
  if(ctrl->id==0){
    sizes = gk_ismalloc(ctrl->num_procs, 0, "malloc sizes");
    displays = gk_ismalloc(ctrl->num_procs, 0, "malloc displays");
    total_error = gk_dsmalloc(ctrl->num_procs,0,"malloc error");
  }

  MPI_Gather(&step, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (ctrl->id == 0) {
    displays[0] = 0;
    for (i = 0; i < ctrl->num_procs; i++) {
      if (i != ctrl->num_procs - 1)
	displays[i + 1] = displays[i] + sizes[i];
    }
  }

  
  double *global_US = gk_dmalloc(train->nrows*ctrl->num_dimensions,"malloc global US");
  double *global_V = gk_dmalloc(train->ncols*ctrl->num_dimensions, "malloc global V");

  
  while(users_diff >= 0.01*(orig_train->nrows)){

    if(ctrl->id==0)
      printf("iter %d\n",iter);
    
    if(ctrl->id==0){
      gk_csr_Free(&train);
      
      int nrows = orig_train->nrows;
      int nnz = orig_train->rowptr[orig_train->nrows] - orig_train->rowptr[0];
      
      train = gk_csr_Create();
      train->nrows = nrows;
      train->ncols = orig_train->ncols;
      train->rowptr = gk_zmalloc(train->nrows+1, "gk_csr_ExtractPartition: rowptr");
      train->rowind = gk_imalloc(nnz, "gk_csr_ExtractPartition: rowind");
      train->rowval = gk_fmalloc(nnz, "gk_csr_ExtractPartition: rowval");
      
      nnz = 0;
      k = 0; 
      train->rowptr[0] = 0;
      for (i=0;i<orig_train->nrows;i++){
	for(j=orig_train->rowptr[i]; j<orig_train->rowptr[i+1]; j++){
	  train->rowind[nnz] = orig_train->rowind[j];
	  train->rowval[nnz] = orig_train->rowval[j]*gu[i];
	  nnz++;
	}
	k++;
	train->rowptr[k] = nnz;
      }
      
      gk_csr_CreateIndex(train, GK_CSR_COL);
      
      global_nnz = train->rowptr[train->nrows] - train->rowptr[0];
      SMat svd_train = (SMat) calloc(1, sizeof(struct smat));
      svd_train->rows = train->nrows;
      svd_train->cols = orig_train->ncols;
      svd_train->vals = global_nnz;
      
      svd_train->pointr = gk_malloc((svd_train->cols + 1)*sizeof(long), "svdNewSMat: pointr");
      svd_train->rowind = gk_malloc(global_nnz*sizeof(long), "svdNewSMat: rowind");
      svd_train->value  = gk_malloc(global_nnz*sizeof(double), "svdNewSMat: value");
      
      global_nnz = 0;
      svd_train->pointr[0] = 0;
      for (i=0;i<orig_train->ncols;i++){
	for(j=train->colptr[i]; j<train->colptr[i+1]; j++){
	  svd_train->rowind[global_nnz] = (long)train->colind[j];
	  svd_train->value[global_nnz] = (double)train->colval[j];
	  global_nnz++;
	}
	svd_train->pointr[i+1] = global_nnz;
      }

      //Compute SVD of rank=ctrl->num_dimensions on the global training matrix
      global_result = svdLAS2A(svd_train, ctrl->num_dimensions);
    
      int index = 0;
      int index2 = 0;
    
      for(j=0; j<train->nrows; j++){
	for(d = 0; d < ctrl->num_dimensions; d++){
	  global_US[index] =  global_result->Ut->value[d][j]*global_result->S[d];
	  index++;
	}
      }
    
      for(d=0; d<ctrl->num_dimensions; d++){
	memcpy(global_V + index2,global_result->Vt->value[d],train->ncols*sizeof(double));
	index2+= train->ncols;
      }

      svdFreeSVDRec(global_result);
      svdFreeSMat(svd_train);
    }
    
    MPI_Bcast(global_US, ctrl->num_dimensions*orig_train->nrows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(global_V, ctrl->num_dimensions*orig_train->ncols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    
    int index = 0;
    int index2 = 0;

    for(cluster_id=startcluster; cluster_id<endcluster; cluster_id++){
      int nrows=0;
      int nnz=0;
      
      for(i=0;i<orig_train->nrows;i++){
	if((assignment[i]==cluster_id)){
	  nrows++;
	  nnz+=orig_train->rowptr[i+1]-orig_train->rowptr[i];
	}
      }
      //Compute the training matrix for cluster_id.
      local_train = gk_csr_Create();
      local_train->nrows = nrows;
      local_train->ncols = orig_train->ncols;
      local_train->rowptr = gk_zmalloc(local_train->nrows+1, "gk_csr_ExtractPartition: rowptr");
      local_train->rowind = gk_imalloc(nnz, "gk_csr_ExtractPartition: rowind");
      local_train->rowval = gk_fmalloc(nnz, "gk_csr_ExtractPartition: rowval");
      
      nnz = 0;
      k = 0; 
      local_train->rowptr[0] = 0;
      for (i=0;i<orig_train->nrows;i++){
	if(assignment[i]==cluster_id){
	  for(j=orig_train->rowptr[i]; j<orig_train->rowptr[i+1]; j++){
	    local_train->rowind[nnz] = orig_train->rowind[j];
	    local_train->rowval[nnz] = orig_train->rowval[j]*(1-gu[i]);
	    nnz++;
	  }
	  k++;
	  local_train->rowptr[k] = nnz;
	}
      }
      
      
      gk_csr_CreateIndex(local_train, GK_CSR_COL);          

      SMat svd_local_train = (SMat) calloc(1, sizeof(struct smat));
      svd_local_train->rows = local_train->nrows;
      svd_local_train->cols = orig_train->ncols;
      svd_local_train->vals = nnz;
      
      svd_local_train->pointr = gk_malloc((svd_local_train->cols + 1)*sizeof(long), "svdNewSMat: pointr");
      svd_local_train->rowind = gk_malloc(nnz*sizeof(long), "svdNewSMat: rowind");
      svd_local_train->value  = gk_malloc(nnz*sizeof(double), "svdNewSMat: value");
      
      nnz = 0;
      svd_local_train->pointr[0] = 0;
      for (i=0;i<orig_train->ncols;i++){
	for(j=local_train->colptr[i]; j<local_train->colptr[i+1]; j++){
	  svd_local_train->rowind[nnz] = (long)local_train->colind[j];
	  svd_local_train->value[nnz] = (double)local_train->colval[j];
	  nnz++;
	}
	svd_local_train->pointr[i+1] = nnz;
      }

      //Compute SVD for the training matrix of cluster cluster_id
      result = svdLAS2A(svd_local_train, dimensions[cluster_id]);
      
      for(j=0; j<nrows; j++){
	for(d = 0; d < dimensions[cluster_id]; d++){
	  partial_US[index] =  result->Ut->value[d][j]*result->S[d];
	  index++;
	}
      }
    
      for(d=0; d<dimensions[cluster_id]; d++){
	memcpy(partial_V + index2,result->Vt->value[d],train->ncols*sizeof(double));
	index2+= train->ncols;
      }
      svdFreeSVDRec(result);
      
      svdFreeSMat(svd_local_train);
      gk_csr_Free(&local_train);
    }
    
    MPI_Gatherv(partial_V, vclustersize, MPI_DOUBLE, long_V,
		vcluster_sizes, vcluster_displays, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Bcast(long_V, total_dimensions*(train->ncols), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Gatherv(partial_US, usclustersize, MPI_DOUBLE, long_US,
		cluster_sizes, cluster_displays, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Bcast(long_US, us_total_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    error = 0;
    for (u=startu;u<endu;u++){
      for(i=0;i<orig_train->ncols;i++)
	local_indices[i] = 0;
      
      for(i=orig_train->rowptr[u]; i<orig_train->rowptr[u+1]; i++){
	local_indices[orig_train->rowind[i]] = 1;
      }
      
      cluster_id = assignment[u];
      pos = 0;
      for(k=0; k<cluster_id; k++){
	pos += train->ncols*dimensions[k];
      }

      
      //learning gu
      nom = 0;
      denom = 0;
      for(i=0;i<orig_train->ncols;i++){
	local[i] = 0;
	global[i] = 0;
	for(d=0; d<dimensions[cluster_id]; d++){ 
	  if(gu[u]!=1) 
	    local[i] += (1/(1-gu[u]))*long_US[local_corr[u] + d]*long_V[pos +d*train->ncols + i]; 
	} 
	for(d=0; d<ctrl->num_dimensions; d++){ 
	  if(gu[u]!=0) 
	    global[i] += (1/(gu[u]))*global_US[u*ctrl->num_dimensions + d]*global_V[d*train->ncols + i];
	}
	nom += (local[i] - global[i]) * (local[i] - local_indices[i]); 
	denom += (local[i] - global[i])*(local[i] - global[i]); 
	error += (gu[u]*global[i] + (1-gu[u])*local[i] - local_indices[i])*
	  (gu[u]*global[i] + (1-gu[u])*local[i] - local_indices[i]);
      } 
      if(denom!=0)
	partial_gu[u-startu] = nom/denom;
      if(partial_gu[u-startu] < 0)
	partial_gu[u-startu] = 0;
      if(partial_gu[u-startu] > 1)
	partial_gu[u-startu] = 1;
      
    }//end of user loop
    

    MPI_Gatherv(partial_gu, step, MPI_DOUBLE, new_gu,
                sizes, displays, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  

    if(ctrl->id==0){
      users_diff = 0;
      for(u=0; u<orig_train->nrows; u++){
	if((gu[u]-new_gu[u]>0.05)||(gu[u]-new_gu[u]<-0.05))
	  users_diff++;
      }
      
      printf("Users changing their gu are %d\n",users_diff);
      
      for(u=0; u<orig_train->nrows; u++){
	gu[u] = new_gu[u];
      }
    
      
      if(users_diff < 0.01*(orig_train->nrows)){
	sprintf(new_gu_name,"new_%s",ctrl->gu_file);
	fpout = gk_fopen(new_gu_name, "w+","gk_writefile");
	for(u=0; u<orig_train->nrows; u++){
	  fprintf(fpout, "%f\n",gu[u]);
	}
	gk_fclose(fpout);
      }
    
    }

    MPI_Gather(&error, 1, MPI_DOUBLE, total_error, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if(ctrl->id==0){
      double sum_error = 0;
      for(i=0; i<ctrl->num_procs; i++)
	sum_error += total_error[i];
    
      printf("The training error is %f\n",sum_error);
    }
    
    MPI_Bcast(gu, orig_train->nrows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&users_diff, 1, MPI_INT, 0, MPI_COMM_WORLD);
    iter++;
  }

  gk_free((void **)&cluster_sizes, &cluster_displays, &vcluster_sizes, &vcluster_displays, LTERM);
  gk_free((void **)&partial_US, &partial_V, &long_V, &long_US, &global_US, &global_V, LTERM);
  gk_free((void **)&local_positions, &local_displays, &local_corr, LTERM);
  gk_free((void **)&partial_gu,&sizes, &displays, &total_error, 
	  &local_indices, &global, &local, LTERM);
  gk_csr_Free(&train);
  gk_csr_Free(&orig_train);
  free_ctrl(ctrl);
  MPI_Finalize();
}
