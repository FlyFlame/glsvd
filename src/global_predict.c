/**************************************************************/
/*! \file global_predict.c
    \brief This file is used to test the performance of rGLSVD and 
    sGLSVD on the test set.

    \author    Evangelia Christakopoulou
    \version   1.0
    \date      2017
    \copyright University of Minnesota, 2017
    \license GPL v2
*/
/**************************************************************/

#include<slim.h>

int main(int argc, char *argv[])
{

  ctrl_t *ctrl;
  gk_csr_t *train;
  gk_csr_t *test;
  int my_id, my_num_procs;
  int i;
  int u;
  int d;
  int j;
  int k;
  int jj;
  int kk;
  size_t lnlen;
  char *line=NULL;
  FILE *fpin=NULL;
  double value = 0;
  int cluster_id = 0;
  gk_csr_t* local_train;
  SVDRec result;
  int global_nnz;
  SVDRec global_result;
  int pos, pos2;
  double global_usv, local_usv;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &my_num_procs);

  // parse command line 
  ctrl = create_ctrl();
  parse_cmdline(ctrl, argc, argv);

  ctrl->id = my_id;
  ctrl->num_procs = my_num_procs;

  // Reading train and test 
  train = gk_csr_Read(ctrl->train_file, GK_CSR_FMT_CSR, 1, 1);
  gk_csr_CreateIndex(train, GK_CSR_COL);

  test = gk_csr_Read(ctrl->test_file, GK_CSR_FMT_CSR, 1, 1);
  gk_csr_CreateIndex(test, GK_CSR_COL);
  
  fpin = gk_fopen(ctrl->dimensions_file, "r", "gk_readfile");

  int dimensions[ctrl->num_clusters];
  int assignment[train->nrows];
  double gu[train->nrows];

  if(ctrl->id == 0){
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
      u++;
    }
    gk_fclose(fpin);
  }
  
  MPI_Bcast(dimensions, ctrl->num_clusters, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(assignment, train->nrows, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(gu, train->nrows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  int total_dimensions = 0;
  for(i=0;i<ctrl->num_clusters;i++){
    total_dimensions += dimensions[i];
  }


  double *global_US = gk_dmalloc(train->nrows*ctrl->num_dimensions,"malloc global US");
  double *global_V = gk_dmalloc(train->ncols*ctrl->num_dimensions, "malloc global V");
  
  if(ctrl->id == 0){
    global_nnz = train->rowptr[train->nrows] - train->rowptr[0];
    SMat svd_train = (SMat) calloc(1, sizeof(struct smat));
    svd_train->rows = train->nrows;
    svd_train->cols = train->ncols;
    svd_train->vals = global_nnz;
    
    svd_train->pointr = gk_malloc((svd_train->cols + 1)*sizeof(long), "svdNewSMat: pointr");
    svd_train->rowind = gk_malloc(global_nnz*sizeof(long), "svdNewSMat: rowind");
    svd_train->value  = gk_malloc(global_nnz*sizeof(double), "svdNewSMat: value");
    
    global_nnz = 0;
    svd_train->pointr[0] = 0;
    for (i=0;i<train->ncols;i++){
      for(j=train->colptr[i]; j<train->colptr[i+1]; j++){
	svd_train->rowind[global_nnz] = (long)train->colind[j];
	svd_train->value[global_nnz] = (double)train->colval[j]*gu[train->colind[j]];
	global_nnz++;
      }
      svd_train->pointr[i+1] = global_nnz;
    }
    
    global_result = svdLAS2A(svd_train, ctrl->num_dimensions);
    svdFreeSMat(svd_train);

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

  }
  
  
  MPI_Bcast(global_US, ctrl->num_dimensions*train->nrows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(global_V, ctrl->num_dimensions*train->ncols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  

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

  int index = 0;
  int index2 = 0;


  for(cluster_id=startcluster; cluster_id<endcluster; cluster_id++){
    int nrows=0;
    int nnz=0;
      
    for(i=0;i<train->nrows;i++){
      if(assignment[i]==cluster_id){
	nrows++;
	nnz+=train->rowptr[i+1]-train->rowptr[i];
      }
    }
    

    local_train = gk_csr_Create();
    local_train->nrows = nrows;
    local_train->ncols = train->ncols;
    local_train->rowptr = gk_zmalloc(nrows+1, "gk_csr_ExtractPartition: rowptr");
    local_train->rowind = gk_imalloc(nnz, "gk_csr_ExtractPartition: rowind");
    local_train->rowval = gk_fmalloc(nnz, "gk_csr_ExtractPartition: rowval");
    
  
    
    nnz = 0;
    k = 0; 
    local_train->rowptr[0] = 0;
    for (i=0;i<train->nrows;i++){
      if(assignment[i]==cluster_id){
	  for(j=train->rowptr[i]; j<train->rowptr[i+1]; j++){
	    local_train->rowind[nnz] = train->rowind[j];
	    local_train->rowval[nnz] = train->rowval[j]*(1-gu[i]);
	    nnz++;
	  }
	k++;
	local_train->rowptr[k] = nnz;
      }
    }  
    
    gk_csr_CreateIndex(local_train,GK_CSR_COL);
    
      
    nnz = local_train->rowptr[local_train->nrows] - local_train->rowptr[0];
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

    gk_csr_Free(&local_train);
    svdFreeSMat(svd_local_train);
  }//end of cluster loop
  

  MPI_Gatherv(partial_V, vclustersize, MPI_DOUBLE, long_V,
	      vcluster_sizes, vcluster_displays, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
  MPI_Bcast(long_V, total_dimensions*(train->ncols), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
  MPI_Gatherv(partial_US, usclustersize, MPI_DOUBLE, long_US,
	      cluster_sizes, cluster_displays, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
  MPI_Bcast(long_US, us_total_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  gk_free((void **)&cluster_sizes, &cluster_displays, &vcluster_sizes, &vcluster_displays, LTERM);
  gk_free((void **)&partial_US, &partial_V, LTERM);
  
  
  double *USV = gk_malloc(train->ncols*sizeof(double), "gk_readfile: US");
  int *local_indices = gk_malloc(train->ncols*sizeof(int),"malloc index"); 
  gk_dkv_t *rcmd = gk_dkvmalloc(train->ncols, "malloc rcmd");
  int *hr_hits = gk_malloc(ctrl->num_clusters*sizeof(int),"malloc hr_hits");
  double *arhr_hits = gk_malloc(ctrl->num_clusters*sizeof(double),"malloc arhr hits");
  double *zero_error = gk_malloc(ctrl->num_clusters*sizeof(double),"malloc zero error");
  double *one_error = gk_malloc(ctrl->num_clusters*sizeof(double),"malloc one error");
  double *rating_norm = gk_malloc(ctrl->num_clusters*sizeof(double),"malloc rating norm");

  for(cluster_id=0; cluster_id < ctrl->num_clusters; cluster_id++){
    hr_hits[cluster_id] = 0;
    arhr_hits[cluster_id] = 0;
    zero_error[cluster_id] = 0;
    one_error[cluster_id] = 0;
    rating_norm[cluster_id] = 0;
  }
  
  
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
  int *local = gk_malloc(train->nrows*sizeof(int),"malloc local correspondence");
  
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
	local[u] = local_displays[cluster_id] + local_positions[cluster_id];
	local_positions[cluster_id]+= dimensions[cluster_id];
      }
    }
  }
  

  for(u=startu; u<endu;u++){
    cluster_id = assignment[u];
    for(i=0;i<train->ncols;i++)
      local_indices[i] = 0;
    
    for(i=train->rowptr[u]; i<train->rowptr[u+1]; i++){
      local_indices[train->rowind[i]] = 1;
    }
    int nrcmd = 0;
    pos = 0;
    pos2 = 0;
    for(k=0; k<cluster_id; k++){
      pos += train->ncols*dimensions[k];
      pos2 += rows[k]*dimensions[k];
    }
    
    for(i=0; i<train->ncols; i++){
      USV[i] = 0;
      global_usv = 0; local_usv = 0; 
      for(d=0; d<ctrl->num_dimensions; d++){
	global_usv += global_US[u*ctrl->num_dimensions + d]*global_V[d*train->ncols + i];
      }
      for(d=0; d<dimensions[cluster_id]; d++){
	local_usv += long_US[local[u] + d]*long_V[pos + d*train->ncols + i];
      }
      USV[i] = global_usv + local_usv;
      if(local_indices[i] == 0){
	rcmd[nrcmd].key = USV[i];
	rcmd[nrcmd].val = i;
	nrcmd++;
	zero_error[cluster_id] += (local_usv+global_usv)*(local_usv+global_usv);
	rating_norm[cluster_id] += global_usv*global_usv;
      } 
      else{
	one_error[cluster_id] += (local_usv - 1 + global_usv)*(local_usv - 1 + global_usv);
	rating_norm[cluster_id] += (1-global_usv)*(1-global_usv);
      }
    }
    gk_dkvsortd(nrcmd, rcmd);
    int nrcmd2 = gk_min(nrcmd, ctrl->topn);
    for (jj = 0; jj < nrcmd2; jj++) {
      for (kk = test->rowptr[u]; kk < test->rowptr[u + 1]; kk++) {
	if (rcmd[jj].val == test->rowind[kk]) {
	  hr_hits[cluster_id] += 1;
	  arhr_hits[cluster_id] += 1.0/(double) (jj + 1);
	}
      }
    }
  }
  

  int *total_hits = NULL;
  double *total_arhr_hits = NULL;
  double *total_zero_error = NULL;
  double *total_one_error = NULL;
  double *total_rating_norm = NULL;
  
  if(ctrl->id==0){
    total_hits = gk_malloc((ctrl->num_clusters*ctrl->num_procs)*sizeof(int),"malloc hits");
    total_arhr_hits = gk_malloc((ctrl->num_clusters*ctrl->num_procs)*sizeof(double),"malloc arhr hits");
    total_zero_error = gk_malloc(ctrl->num_clusters*ctrl->num_procs*sizeof(double),"malloc zero error");
    total_one_error = gk_malloc(ctrl->num_clusters*ctrl->num_procs*sizeof(double),"malloc one error");
    total_rating_norm = gk_malloc(ctrl->num_clusters*ctrl->num_procs*sizeof(double),"malloc rating norm");
  }
  
  
  MPI_Gather(zero_error, ctrl->num_clusters, MPI_DOUBLE, total_zero_error, ctrl->num_clusters, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(one_error, ctrl->num_clusters, MPI_DOUBLE, total_one_error, ctrl->num_clusters, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(hr_hits, ctrl->num_clusters, MPI_INT, total_hits, ctrl->num_clusters, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(arhr_hits, ctrl->num_clusters, MPI_DOUBLE, total_arhr_hits, ctrl->num_clusters, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(rating_norm, ctrl->num_clusters, MPI_DOUBLE, total_rating_norm, ctrl->num_clusters, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  int *hr = gk_malloc(ctrl->num_clusters*sizeof(int),"malloc hr_hits");
  double *arhr = gk_malloc(ctrl->num_clusters*sizeof(double),"malloc arhr hits");

  if(ctrl->id==0){
    for(cluster_id = 0; cluster_id < ctrl->num_clusters; cluster_id++){
      hr[cluster_id] = 0;
      arhr[cluster_id] = 0;
      for(i=0; i<ctrl->num_procs; i++){
	hr[cluster_id] += total_hits[i*ctrl->num_clusters+cluster_id];
	arhr[cluster_id] += total_arhr_hits[i*ctrl->num_clusters+cluster_id];
      }
      printf("Cluster %d HR is %d ARHR is %f\n",cluster_id,hr[cluster_id],arhr[cluster_id]);
    }
  }
  
  
  gk_free((void **)&long_V,&long_US,LTERM);
  gk_free((void **)&local_positions,&local_displays,&local,LTERM);
  gk_free((void **)&local_indices, &rcmd, LTERM);
  
  gk_csr_Free(&train);
  gk_csr_Free(&test);
  free_ctrl(ctrl);
  MPI_Finalize();
}
  
