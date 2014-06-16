#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
//#include <string.h>
#include <R_ext/Utils.h> /* for R_CheckUserInterrupt() */
#include <R_ext/Boolean.h> /* for R_CheckUserInterrupt() */

#ifdef _OPENMP
#include <omp.h>
#endif

#define ULONGLONG uint64_t

Rboolean R_ToplevelExec(void (*fun)(void *), void *data);
int checkInterrupt(void);
static void chkIntFn(void *);
int FLAG_interrupt=0;
ULONGLONG interrupt_counter=0;


extern "C"
void rec(int **x, int *dim, int *pos, int *margins, ULONGLONG *local_countTables, double *local_probTables, ULONGLONG *local_nO000, double *preCalcFact, int *rowsums, int *colsums, double diff_lfmarginsTotal_lfN, int *O000, double *p0, int this_thread){

	if(FLAG_interrupt){return;} //2**18-1
	if( this_thread == 0 && (((++interrupt_counter) & 262143) == 0 ) && (checkInterrupt()==1)){ //check if thread is master to check for user interrupt
		FLAG_interrupt=1; interrupt_counter=0; return;
	}
	if(interrupt_counter > 1000000000){interrupt_counter=0;}


	int some_columns_left=0;
	int current_row=pos[0];
	int current_col=pos[1];
	if( dim[1] - 1 - current_col == 0 ){ // reached last column
		int diff_row=margins[current_row] - rowsums[current_row];
		if(diff_row<0){return;}

		if( dim[0] - 1 - current_row == 0 ){ // reached last row

			int diff_col=margins[dim[0]+current_col] - colsums[current_col];
			if( diff_row == diff_col && diff_row >= 0){
				x[current_row][current_col]=diff_col;
				//rowsums[current_row] += diff_col;
				//colsums[current_col] += diff_col;

				double lmm=0.0;
				for(int i=0; i<dim[0]; ++i){for(int j=0; j<dim[1]; ++j){lmm += preCalcFact[ x[i][j] ];}}

				double prob=exp( (diff_lfmarginsTotal_lfN - lmm) );
				local_probTables[0] += prob; ++local_countTables[0];
				if( x[0][0] == *O000 ){ ++local_nO000[0]; }
				if( x[0][0] <= *O000 ){ local_probTables[1] += prob; ++local_countTables[1]; } //less
				if( x[0][0] >= *O000 ){ local_probTables[2] += prob; ++local_countTables[2]; } //greater
				if( prob <= *p0 ){ local_probTables[3] += prob; ++local_countTables[3]; } //two sided, minimum likelihood
			}

			return;
		}else{ // some rows still left

			x[current_row][current_col]=diff_row;
			rowsums[current_row] += diff_row;
			colsums[current_col] += diff_row;

			++current_row; // go to next row
			current_col=0; // go to first column in next row

			some_columns_left=1;
		}
	}else{ // some columns still left
		some_columns_left=1;
	}


	if(some_columns_left ){ // explore next column recursively
		int *rowsums0=(int*)malloc(dim[0]*sizeof(int)); //new int[dim[0]];
		int *colsums0=(int*)malloc(dim[1]*sizeof(int)); //new int[dim[1]];
		int *pos0=(int*)malloc(2*sizeof(int)); //new int[2];
		pos0[0]=current_row;
		pos0[1]=current_col+1;
		int is_last_row=0; if(dim[0]-1-current_row == 0){is_last_row=1;}

		int diff_row=margins[ current_row ] - rowsums[current_row];
		int diff_col=margins[ dim[0]+current_col ] - colsums[current_col];
		int minAtCurrentPosition=( (diff_row < diff_col) ? diff_row : diff_col );
		for(int i = 0; i<=minAtCurrentPosition; i++ ){
			x[current_row][current_col]=i;
			memcpy((void*)colsums0, (void*)colsums, dim[1]*sizeof(int)); colsums0[current_col] += i;

			// check for valid column margin as early as possible to avoid unnecassary recursive function calls
			int tmp=colsums0[current_col] - margins[dim[0]+current_col];
			if( (is_last_row ) && ( tmp != 0 ) ){
				if( tmp < 0 ){continue;}
				break;
			}

			memcpy((void*)rowsums0, (void*)rowsums, dim[0]*sizeof(int)); rowsums0[current_row] += i;

			rec( x, dim, pos0, margins,local_countTables, local_probTables,local_nO000, preCalcFact, rowsums0, colsums0 ,diff_lfmarginsTotal_lfN, O000, p0, this_thread);
		}
		free(rowsums0); free(colsums0); free(pos0);
	}

return;
}


extern "C"
void hypergeom_IxJ(int *O000, int *N, int *margins, double *p0, double *n0, double *Prob, double *Freq, int *Hist_n000, int *dim, int *nthreads){
	int NN=*N;
	int h=0;
	double *preCalcFact= new double[NN+1];
	double lfN=0.0;
	double lfmargins[2];
	double lfmarginsTotal=0.0;

	//results
	ULONGLONG countTables[4];
	double probTables[4];
	ULONGLONG nO000=0;

	for( h=0; h<4; h++ ){ countTables[h]=0; probTables[h]=0.0; }

	preCalcFact[0]=0;
	for( h=1; h<= NN; ++h ){ preCalcFact[h]=preCalcFact[h-1]+log(h); }
	lfN=preCalcFact[NN];
	lfmargins[0]=0.0;
	for( h=0; h<dim[0]; ++h ){ lfmargins[0] += preCalcFact[ margins[h] ]; }
	lfmargins[1]=0.0;
	for( h=0; h<dim[1]; ++h ){ lfmargins[1] += preCalcFact[ margins[h+dim[0]] ]; }
	lfmarginsTotal=lfmargins[0]+lfmargins[1];

	double diff_lfmarginsTotal_lfN=lfmarginsTotal-lfN;

	#ifdef _OPENMP
	if(*nthreads <= 1){ *nthreads=1; }else{ *nthreads=(*nthreads < omp_get_max_threads()) ? (*nthreads) : (omp_get_max_threads()); }
	#endif


	// state is stored after last run; re-initialize
	FLAG_interrupt=interrupt_counter=0;

	int minAtFirstPosition=( (margins[ 0 ] < margins[ dim[1] ]) ? margins[ 0 ] : margins[ dim[1] ] );

	int *aux_rowsums=new int[dim[0]];
	int *aux_colsums=new int[dim[1]];
	for( int hh=0; hh<dim[0]; ++hh ){aux_rowsums[hh]=0;}
	for( int gg=0; gg<dim[1]; ++gg ){ aux_colsums[gg]=0;}

	#pragma omp parallel shared(countTables, probTables, nO000) num_threads(*nthreads)
	{
		int this_thread=0;
		#ifdef _OPENMP
		this_thread=omp_get_thread_num();
		#endif

		ULONGLONG *local_nO000=new ULONGLONG[1];
		local_nO000[0]=0;
		ULONGLONG *local_countTables=new ULONGLONG[4];
		double *local_probTables=new double[4];
		for( h=0; h<4; ++h ){ local_countTables[h]=0; local_probTables[h]=0.0; }

		int i1=0;

		#pragma omp for schedule(dynamic)
		for( i1=0; i1<=minAtFirstPosition; ++i1 ){

			if(FLAG_interrupt){
				#ifdef _OPENMP
				continue;
				#else
				break;
				#endif

			}
			if( (this_thread == 0) && (((++interrupt_counter) & 262143) == 0 ) && (checkInterrupt()==1) ){ //check if thread is master to check for user interrupt
				FLAG_interrupt=1; continue;
			}

			int *pos=new int[2];
			int *rowsums=new int[dim[0]];
			int *colsums=new int[dim[1]];
			int **x=new int*[ dim[0] ];
			memcpy((void*)rowsums, (void*)aux_rowsums, (dim[0])*sizeof(int));
			memcpy((void*)colsums, (void*)aux_colsums, (dim[1])*sizeof(int));
			for( int hh=0; hh<dim[0]; hh++ ){
				x[hh]= new int[ dim[1] ];
				for( int gg=0; gg<dim[1]; gg++ ){ x[hh][gg]= 0; }
			}

			pos[0]=0; pos[1]=1; //re-initialize in each iteration

			x[0][0]=i1;
			rowsums[0]=colsums[0]=x[0][0];

			rec( x, dim, pos, margins,local_countTables, local_probTables,local_nO000, preCalcFact, rowsums, colsums ,diff_lfmarginsTotal_lfN, O000, p0, this_thread);

			for( int hh=0; hh<dim[0]; hh++ ){ delete[] x[hh]; }; delete[] x;
			delete[] pos;
			delete[] rowsums;
			delete[] colsums;
		} //end loop i1


		#pragma omp critical
		{
			for( h=0; h<4; ++h ){
				countTables[h] += local_countTables[h];
				probTables[h] += local_probTables[h];
			}
			nO000 += local_nO000[0];
		}

		delete[] local_countTables;
		delete[] local_probTables;
		delete[] local_nO000;
	}

	for( h=0; h<4; ++h ){
		Freq[h] += countTables[h];
		Prob[h] += probTables[h];
	}
	double min=Prob[1];
	ULONGLONG min_count=Freq[1];
	if( min>Prob[2] ){ min=Prob[2];min_count=Freq[2]; }
	Prob[4]=2.0*min;
	Freq[4]=2.0*min_count;

	*n0=nO000;

	delete[] preCalcFact;
	delete[] aux_colsums;
	delete[] aux_rowsums;
}

