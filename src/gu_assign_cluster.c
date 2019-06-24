/**************************************************************/
/*! \file gu_assign_cluster.c
    \brief This file implements the algorithm sGLSVD.

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
  char new_participation_name[200]="";
  int user;
  int k, nnz;
  int min_assignment;
  double min_error;
  SVDRec result;
  SVDRec global_result;
  int iter;
  int users_diff;
  int pos, pos2;
  double nom, denom, weight, min_gu;
  double obj;
  
    
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
  int prev_assignment[train->nrows];
  double gu[train->nrows];

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
      prev_assignment[u] = value; 
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

  //Setup MPI  
  MPI_Bcast(dimensions, ctrl->num_clusters, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(assignment, train->nrows, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(prev_assignment, train->nrows, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(gu, train->nrows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  
  double *gprojection = gk_dmalloc(ctrl->num_dimensions,"projection for global");
  int *local_indices = gk_malloc(orig_train->ncols*sizeof(int),"malloc index");
  double *global = gk_malloc(orig_train->ncols*sizeof(double),"global");
  double *local = gk_malloc(orig_train->ncols*sizeof(double),"local");
  double *error= gk_malloc(orig_train->nrows*sizeof(double),"error");
  double *local_error = gk_malloc((ctrl->num_clusters)*sizeof(double),"error");
  double *cluster_gu = gk_malloc((ctrl->num_clusters)*sizeof(double),"cluster gu");
  int *sizes = NULL;
  int *displays = NULL;
  int *cluster_sizes = NULL;
  int *cluster_displays = NULL;
  int *vcluster_sizes = NULL;
  int *vcluster_displays = NULL;


  iter = 0;
  users_diff = orig_train->nrows;

  int datasize = orig_train->nrows;
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
  double *partial_gu = gk_dmalloc(step, "malloc partial gu");
  double *partial_error = gk_dmalloc(step, "malloc partial error");
  if (ctrl->id == 0) {
    displays = gk_ismalloc(ctrl->num_procs, 0, "malloc displays");
    sizes = gk_ismalloc(ctrl->num_procs, 0, "malloc sizes");
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

  double *global_S = gk_dmalloc(ctrl->num_dimensions,"malloc global US");
  double *global_V = gk_dmalloc(orig_train->ncols*ctrl->num_dimensions, "malloc global V");

  
  while(users_diff >= 0.01*(orig_train->nrows)){
    int total_dimensions = 0;
    for(i=0;i<ctrl->num_clusters;i++){
      total_dimensions += dimensions[i];
    }
    
    int sclustersize = 0;
    for(i=startcluster; i<endcluster; i++){
      sclustersize += dimensions[i];
    }
    int vclustersize = sclustersize*orig_train->ncols;
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
    gk_csr_Free(&train);
  
    double *long_V = gk_malloc(total_dimensions*orig_train->ncols*sizeof(double),"malloc long V");
    double *long_SV = gk_malloc(total_dimensions*orig_train->ncols*sizeof(double),"malloc long V");
    double *long_S = gk_malloc(total_dimensions*sizeof(double),"malloc long S");
    double *long_projection = gk_malloc(total_dimensions*sizeof(double),"malloc long S");
    
    
    if(ctrl->id==0)
      printf("iter %d\n",iter);
    
    if(ctrl->id == 0){
      
      train = gk_csr_Create();
      train->nrows = orig_train->nrows;
      train->ncols = orig_train->ncols;
      train->rowptr = gk_zmalloc(train->nrows+1, "gk_csr_ExtractPartition: rowptr");
      train->rowind = 
	gk_imalloc(orig_train->rowptr[orig_train->nrows] - orig_train->rowptr[0], 
		   "gk_csr_ExtractPartition: rowind");
      train->rowval = gk_fmalloc(orig_train->rowptr[orig_train->nrows] - 
				 orig_train->rowptr[0], "gk_csr_ExtractPartition: rowval");
      
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
      
      int global_nnz = train->rowptr[train->nrows] - train->rowptr[0];
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
      
      //Compute SVD of rank ctrl->num_dimensions on the global training matrix.
      global_result = svdLAS2A(svd_train, ctrl->num_dimensions);
  

      int index2 = 0;

      for(d=0; d<ctrl->num_dimensions; d++){
	global_S[d] = global_result->S[d];
        memcpy(global_V + index2,global_result->Vt->value[d],train->ncols*sizeof(double));
        index2+= train->ncols;
      }

      gk_csr_Free(&train);
      svdFreeSMat(svd_train);
      svdFreeSVDRec(global_result);
      
      
    }
    
    MPI_Bcast(global_S, ctrl->num_dimensions, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(global_V, ctrl->num_dimensions*orig_train->ncols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    int spos = 0;
    for(cluster_id=startcluster; cluster_id<endcluster; cluster_id++){
      int nrows=0;
      int nnz=0;
      
      for(i=0;i<orig_train->nrows;i++){
	if(assignment[i]==cluster_id){
	  nrows++;
	  nnz+=orig_train->rowptr[i+1]-orig_train->rowptr[i];
	}
      }
      
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
      gk_csr_Free(&local_train);

      //Compute SVD on the training matrix for cluster cluster_id
      result = svdLAS2A(svd_local_train, dimensions[cluster_id]);

      svdFreeSMat(svd_local_train);
      k= spos*orig_train->ncols;
      for(d = 0; d < dimensions[cluster_id]; d++){
	for(j=0; j<orig_train->ncols; j++){
	  partial_SV[k] =  result->Vt->value[d][j]*result->S[d];
	  k++;
	}
      }
      
      
      memcpy(partial_S + spos,result->S,dimensions[cluster_id]*sizeof(double));
      
      i = spos*orig_train->ncols;
      for(d=0; d<dimensions[cluster_id]; d++){
	memcpy(partial_V + i,result->Vt->value[d],orig_train->ncols*sizeof(double));
	i+= orig_train->ncols;
      }
      
      
      svdFreeSVDRec(result);
      spos += dimensions[cluster_id];
      
    }//end cluster loop


    MPI_Gatherv(partial_S, sclustersize, MPI_DOUBLE, long_S,
		cluster_sizes, cluster_displays, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Bcast(long_S, total_dimensions, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Gatherv(partial_V, vclustersize, MPI_DOUBLE, long_V,
		vcluster_sizes, vcluster_displays, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Bcast(long_V, total_dimensions*(orig_train->ncols), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Gatherv(partial_SV, vclustersize, MPI_DOUBLE, long_SV,
		vcluster_sizes, vcluster_displays, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Bcast(long_SV, total_dimensions*(orig_train->ncols), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //Assign user to the cluster with the smallest training error and update
    //his corresponding weight
    for(user= startu; user<endu; user++){
      pos = 0;
      pos2 = 0;
      
      for(d = 0; d < ctrl->num_dimensions; d++){
	gprojection[d] = 0;
	if(gu[user]!=0){
	  for(j=orig_train->rowptr[user]; j<orig_train->rowptr[user+1]; j++){
	    k = orig_train->rowind[j];
	    gprojection[d] +=  global_V[d*orig_train->ncols + k]*(1/global_S[d]);  
	  }
	}
      }
      
      for(i=0; i<orig_train->ncols; i++){
	global[i] = 0;
	for(d=0; d<ctrl->num_dimensions; d++){
	  global[i] += gprojection[d]*global_S[d]*global_V[d*orig_train->ncols + i];
	}
      }
      
      for(i=0;i<orig_train->ncols;i++)
	local_indices[i] = 0;
      
      for(i=orig_train->rowptr[user]; i<orig_train->rowptr[user+1]; i++)
	local_indices[orig_train->rowind[i]] = 1;
    
      //Trying the different clusters
      for(cluster_id=0; cluster_id<ctrl->num_clusters; cluster_id++){
	if(cluster_id > 0){
	  pos += dimensions[cluster_id-1]*(orig_train->ncols);
	  pos2 += dimensions[cluster_id-1];
	}
	for(d = 0; d < dimensions[cluster_id]; d++){
	  long_projection[pos2 + d] = 0;
	  if(gu[user]!=1){
	    for(j=orig_train->rowptr[user]; j<orig_train->rowptr[user+1]; j++){
	      k = orig_train->rowind[j];
	      long_projection[pos2 + d] +=  long_V[pos + d*orig_train->ncols + k]*(1/long_S[pos2 + d]);  
	    }
	  }
	}
	
	for(i=0;i<orig_train->ncols;i++){
	  local[i] = 0;
	  for(d=0; d<dimensions[cluster_id]; d++){
	    local[i] += long_projection[pos2 + d]*long_SV[pos + d*orig_train->ncols + i];
	  }
	}
	
	nom = 0;
	denom = 0;
	for(i=0; i<orig_train->ncols; i++){
	  nom += (local[i] - global[i]) * (local[i] - local_indices[i]);
	  denom += (local[i] - global[i])*(local[i] - global[i]); 
	}
	
	if(denom!=0)
	  weight = nom/denom;
	else
	  weight = 0;
	if(weight < 0)
	  weight = 0;
	if(weight > 1)
	  weight = 1;

	cluster_gu[cluster_id] = weight;
	
	local_error[cluster_id] = 0;
	for(i=0; i<orig_train->ncols; i++){
	  local_error[cluster_id] += 
	    (weight*global[i] + (1-weight)*local[i] - local_indices[i])*
	    (weight*global[i] + (1-weight)*local[i] - local_indices[i]);
	}
      }//end loop of clusters

      min_assignment = assignment[user]; 
      min_error = local_error[min_assignment];
      min_gu = cluster_gu[min_assignment];
      
      for(cluster_id=0; cluster_id <ctrl->num_clusters; cluster_id++){
	if(local_error[cluster_id]<min_error){
	  min_error = local_error[cluster_id];
	  min_assignment = cluster_id;
	  min_gu = cluster_gu[cluster_id];
	}
      }
      
      partial_error[user-startu] = min_error; 
      partial_participation[user-startu] = min_assignment;
      partial_gu[user-startu] = min_gu;

    }//end loop for users
    
    MPI_Gatherv(partial_error, step, MPI_DOUBLE, error,
                sizes, displays, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Gatherv(partial_participation, step, MPI_INT, assignment,
		sizes, displays, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(assignment, orig_train->nrows, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Gatherv(partial_gu, step, MPI_DOUBLE, gu, 
		sizes, displays, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(gu, orig_train->nrows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    users_diff = 0;
    for(user=0; user<orig_train->nrows; user++){
      if(assignment[user]!=prev_assignment[user])
	users_diff++;
    }
    
    if(ctrl->id==0){
      obj = 0;
      for(user=0; user<orig_train->nrows; user++){
	obj += error[user];
      }
      printf("Objective is %f\n",obj);
      printf("Users jumping clusters are %d\n",users_diff);
      
    }    
    
    for(user=0;user<orig_train->nrows; user++){
      prev_assignment[user] = assignment[user];
    }
    
    
    if(ctrl->id==0){
      sprintf(new_participation_name,"new_%s",ctrl->participation_file); 
      //    sprintf(new_participation_name,"new_%s_%d",ctrl->participation_file,iter);
      fpout = gk_fopen(new_participation_name, "w+","gk_writefile");
      for(user=0; user<orig_train->nrows; user++){
	fprintf(fpout, "%d\n",assignment[user]);
      }
      gk_fclose(fpout);
      
      sprintf(new_gu_name,"new_%s",ctrl->gu_file); 
      //      sprintf(new_gu_name,"new_%s_%d",ctrl->gu_file,iter);
      fpout = gk_fopen(new_gu_name, "w+","gk_writefile");
      for(user=0; user<orig_train->nrows; user++){
        fprintf(fpout, "%f\n",gu[user]);
      }
      gk_fclose(fpout);

    } 

 
    iter++;
    gk_free((void **)&long_S, &long_V, &long_projection, &long_SV, LTERM);
    gk_free((void **)&partial_S, &partial_V, &partial_SV, LTERM);
    
  }

  gk_free((void **)&global_S, &global_V, LTERM);
  gk_free((void **)&local_indices, &global, &local, &gprojection, LTERM);
  gk_free((void **)&sizes, &displays, &cluster_sizes, &cluster_displays, &vcluster_sizes, &vcluster_displays, LTERM);
  gk_free((void **)&partial_participation, &partial_gu, &partial_error, LTERM);
  gk_free((void **)&local_error, &error, &cluster_gu, LTERM);
  gk_csr_Free(&orig_train);
  free_ctrl(ctrl);
  MPI_Finalize();
}
