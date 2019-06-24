/**************************************************************/
/*! \file predict.c
    \brief This file implements the algorithm rLSVD.
    It is also used to evaluate the performance of either rLSVD or 
    sLSVD on the test set.

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
  gk_csr_t* local_test;
  SVDRec result;

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

  test = gk_csr_Read(ctrl->test_file, GK_CSR_FMT_CSR, 1, 1);
  gk_csr_CreateIndex(test, GK_CSR_COL);
  
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
  fpin = gk_fopen(ctrl->participation_file, "r", "gk_readfile");
  u = 0;
  while(gk_getline(&line, &lnlen, fpin) != -1) {
    sscanf(line, "%le", &value);
    assignment[u] = value;
    u++;
  }
  gk_fclose(fpin);

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
    
    local_test = gk_csr_Create();
    local_test->nrows = nrows;
    local_test->ncols = train->ncols;
    local_test->rowptr = gk_zmalloc(local_test->nrows+1, "gk_csr_ExtractPartition: rowptr");
    local_test->rowind = gk_imalloc(nrows, "gk_csr_ExtractPartition: rowind");
    local_test->rowval = gk_fmalloc(nrows, "gk_csr_ExtractPartition: rowval");

    
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
    
    gk_csr_CreateIndex(local_train,GK_CSR_COL);
    
    nnz = 0;
    k = 0; 

    local_test->rowptr[0] = 0;
    for (i=0;i<test->nrows;i++){
      if(assignment[i]==cluster_id){
        for(j=test->rowptr[i]; j<test->rowptr[i+1]; j++){
          local_test->rowind[nnz] = test->rowind[j];
          local_test->rowval[nnz] = test->rowval[j];
          nnz++;
        }
        k++;
        local_test->rowptr[k] = nnz;
      }
    }
    
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
    svdFreeSMat(svd_local_train);

    double *USV = gk_malloc(train->ncols*sizeof(double), "gk_readfile: US");
    int *local_indices = gk_malloc(train->ncols*sizeof(int),"malloc index"); 
    gk_dkv_t *rcmd = gk_dkvmalloc(train->ncols, "malloc rcmd");
    int hr_hits = 0;
    double arhr_hits = 0;
    
    for(u=0; u<nrows;u++){
      for(i=0;i<train->ncols;i++)
	local_indices[i] = 0;
      
      for(i=local_train->rowptr[u]; i<local_train->rowptr[u+1]; i++){
	local_indices[local_train->rowind[i]] = 1;
      }
      int nrcmd = 0;
      for(i=0; i<train->ncols; i++){
	USV[i] = 0;
	if(local_indices[i] == 0){
	  for(d=0; d<dimensions[cluster_id]; d++){
	    USV[i] += result->Ut->value[d][u]*result->S[d]*result->Vt->value[d][i];
	  }
	  rcmd[nrcmd].key = USV[i];
	  rcmd[nrcmd].val = i;
	  nrcmd++;
	}   
      }
      gk_dkvsortd(nrcmd, rcmd);

      int nrcmd2 = gk_min(nrcmd, ctrl->topn);
      for (jj = 0; jj < nrcmd2; jj++) {
	for (kk = local_test->rowptr[u]; kk < local_test->rowptr[u + 1]; kk++) {
	  if (rcmd[jj].val == local_test->rowind[kk]) {
	    hr_hits += 1;
	    arhr_hits += 1.0/(double) (jj + 1);
	  }
	}
      }
    }
    
    printf("Cluster_id: %d HR_hits: %d ARHR_hits: %f\n",cluster_id+1, hr_hits, arhr_hits);
    
    svdFreeSVDRec(result);
    gk_csr_Free(&local_train);
    gk_csr_Free(&local_test);
    gk_free((void **)&USV, LTERM);
    gk_free((void **)&local_indices, &rcmd, LTERM);
  }
  gk_csr_Free(&train);
  gk_csr_Free(&test);
  free_ctrl(ctrl);
  MPI_Finalize();
}
  
