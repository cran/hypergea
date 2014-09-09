#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>


#ifdef _OPENMP
#include <omp.h>
#endif

#define ULONGLONG uint64_t

// #include <stdio.h>
// void printTable(int **x, int *dim, int *rowsums, int *colsums, int *margins){
// 		printf("-------------------------------------------------------\n");
// 		for( int i=0; i<dim[0]; i++ ){
// 			for( int j=0; j<dim[1]; j++ ){
// 				printf("%i\t", x[i][j]);
// 			}
// 			printf("| %i\t| %i\n", rowsums[i], margins[i]);
// 		}
// 		printf("-------------------------------------------------------\n");
// 		for( int i=0; i<dim[1]; i++ ){ printf( "%i\t", colsums[i] ); }
// 		printf("\n-------------------------------------------------------\n");
// 		for( int i=0; i<dim[1]; i++ ){ printf( "%i\t", margins[dim[0]+i] ); }
// 		printf("\n=======================================================\n");
// 
// 
// }


struct StructConstants{
	int *dim;
	int *margins;
	double *preCalcFact;
	double diff_lfmarginsTotal_lfN;
	double *p0;
	int *O000;
};


void update( int **x, int *dim, ULONGLONG *countTables, double *probTables, ULONGLONG *nO000, double *preCalcFact, int *O000, double *p0, double diff_lfmarginsTotal_lfN ){
	double lmm=0.0;
	for(int i=0; i<dim[0]; ++i){for(int j=0; j<dim[1]; ++j){lmm += preCalcFact[ x[i][j] ];}}
	double prob=(diff_lfmarginsTotal_lfN - lmm) < -708 ? 0.0 : exp(diff_lfmarginsTotal_lfN - lmm); // Taylor series approximation beyond x=-708 always zero (at double precision), but extremely slow
	probTables[0] += prob; ++countTables[0];
	if( x[0][0] == *O000 ){ ++nO000[0]; }
	if( x[0][0] <= *O000 ){ probTables[1] += prob; ++countTables[1]; } //less
	if( x[0][0] >= *O000 ){ probTables[2] += prob; ++countTables[2]; } //greater
	if( prob <= *p0 ){ probTables[3] += prob; ++countTables[3]; } //two sided, minimum likelihood
}




extern "C"
void run2x2(struct StructConstants *Constants, ULONGLONG *countTables, double *probTables, ULONGLONG *nO000, int *nthreads){

	int *margins=Constants->margins;
	#pragma omp parallel shared(countTables, probTables, nO000) num_threads(*nthreads) if(*nthreads > 1)
	{
		ULONGLONG *local_nO000=new ULONGLONG[1];
		local_nO000[0]=0;
		ULONGLONG *local_countTables=new ULONGLONG[4];
		double *local_probTables=new double[4];
		for( int h=0; h<4; ++h ){ local_countTables[h]=0; local_probTables[h]=0.0; }
	
		int min00=( (margins[0] < margins[2]) ? margins[0] : margins[2] );
		#pragma omp for schedule(dynamic)
		for( int i=0; i<=min00; ++i ){
			int *rowsums=new int[2];
			int *colsums=new int[2];
			int **x=new int*[ 2 ];
			for( int hh=0; hh<2; hh++ ){
				x[hh]= new int[2 ]; rowsums[hh]=colsums[hh]=0;
			}
			x[0][0]=x[1][0]=x[0][1]=x[1][1]=-1;
			
			x[0][0]=i; 
			rowsums[0] += x[0][0]; colsums[0] += x[0][0];
			x[0][1]=( (margins[0] - i < margins[3]) ? margins[0] - i : margins[3] );
			if( x[0][0] + x[0][1] != margins[0] ){ continue; }
			rowsums[0] += x[0][1]; colsums[1] += x[0][1];
			x[1][0]=margins[2]-i;
			rowsums[1] += x[1][0]; colsums[0] += x[1][0];
			x[1][1]=margins[3]-x[0][1];
			rowsums[1] += x[1][1]; colsums[1] += x[1][1];
			
			update( x, Constants->dim, countTables, probTables, nO000
					, Constants->preCalcFact, Constants->O000, Constants->p0, Constants->diff_lfmarginsTotal_lfN );
			
			delete[] rowsums; delete[] colsums;
			delete[] x[0]; delete[] x[1]; delete[] x;
		}
		
		#pragma omp critical
		{
			for( int h=0; h<4; ++h ){
				countTables[h] += local_countTables[h];
				probTables[h] += local_probTables[h];
			}
			nO000 += local_nO000[0];
		}

		delete[] local_countTables;
		delete[] local_probTables;
		delete[] local_nO000;
	}
}



/*
 * 
 * 	CASE I x 2
 * 
 */



extern "C"
void recIx2(struct StructConstants *Constants, int **x, int i, int *rowsums, int *colsums, ULONGLONG *countTables, double *probTables, ULONGLONG *nO000){

	int I=Constants->dim[0];
	int *dim=Constants->dim;
	int *margins=Constants->margins;
	int remain_col=margins[ I ] - colsums[0]; // margins[I]==margin of first column
	if( i == I-1 ){ // reached last row
		x[ i ][ 0 ]=remain_col;
		if( margins[ i ] >= x[ i ][ 0 ] ){ 
			x[ i ][ 1 ]=margins[ i ]-x[ i ][ 0 ];
			colsums[ 0 ] += x[ i ][ 0 ];
			colsums[ 1 ] += x[ i ][ 1 ];
			rowsums[ i ] = x[i][0]+x[i][1];
			if( margins[ I+1 ] == colsums[ 1 ] ){
				update( x, dim, countTables, probTables, nO000
					, Constants->preCalcFact, Constants->O000, Constants->p0, Constants->diff_lfmarginsTotal_lfN );
			}
		}
	}else{
		int *colsums0=(int*)malloc(2*sizeof(int));
		int *rowsums0=(int*)malloc(I*sizeof(int));
		int **xx=new int*[ I ];
		for( int ii=0; ii<I; ++ii ){ 
			xx[ii]=new int[2];
			xx[ii][0]=x[ii][0];
			xx[ii][1]=x[ii][1];
		}
		for( int h = 0; h<=remain_col; ++h ){
			memcpy((void*)colsums0, (void*)colsums, 2*sizeof(int));
			memcpy((void*)rowsums0, (void*)rowsums, I*sizeof(int)); 
			if( margins[ i ] - h < 0 ){ continue; } // skip if invalid margin
			xx[i][0]=h; xx[i][1]=margins[ i ] - xx[i][0]; 
			colsums0[0] += xx[i][0]; colsums0[1] += xx[i][1];
			rowsums0[i] = xx[i][0]+xx[i][1];
			
			if( margins[ I+1 ] < colsums0[1] ){ continue; } // skip if invalid margin
			
			recIx2( Constants, xx, i+1, rowsums0, colsums0, countTables, probTables, nO000 );
		}
		delete[] rowsums0;
		delete[] colsums0;
		for( int hh=0; hh<dim[0]; ++hh ){ delete[] xx[hh]; }; delete[] xx;

	}
}


extern "C"
void runIx2(struct StructConstants *Constants, ULONGLONG *countTables, double *probTables, ULONGLONG *nO000, int *nthreads){

	int *margins=Constants->margins;
	int *dim=Constants->dim;
	int min00=( (margins[0] < margins[ dim[0] ]) ? margins[0] : margins[ dim[0] ] );
	#pragma omp parallel shared(countTables, probTables, nO000) num_threads(*nthreads) if(*nthreads > 1)
	{
		ULONGLONG *local_nO000=new ULONGLONG[1];
		local_nO000[0]=0;
		ULONGLONG *local_countTables=new ULONGLONG[4];
		double *local_probTables=new double[4];
		for( int h=0; h<4; ++h ){ local_countTables[h]=0; local_probTables[h]=0.0; }
	
		#pragma omp for schedule(dynamic)
		for( int i=0; i<=min00; ++i ){
			int *rowsums=new int[ dim[0] ];
			int *colsums=new int[2];
			int **x=new int*[ dim[0] ];
			for( int hh=0; hh<dim[0]; ++hh ){
				x[hh]= new int[ 2 ]; rowsums[hh]=0;
				for( int gg=0; gg<2; gg++ ){ x[hh][gg]= -1; }
			}
			colsums[0]=colsums[1]=0;

			x[0][0]=i; 
			rowsums[0] += x[0][0]; colsums[0] += x[0][0];
			x[0][1]=margins[0]-x[0][0];
			rowsums[0] += x[0][1]; colsums[1] += x[0][1];
			if( x[0][0] + x[0][1] != margins[0] ){ continue; }
			if( colsums[0] > margins[ dim[0] ] ){ continue; }
			if( colsums[1] > margins[ dim[0]+1 ] ){ continue; }
			recIx2( Constants, x, 1, rowsums, colsums, local_countTables, local_probTables, local_nO000 );
			
			delete[] rowsums; delete[] colsums;
			for( int hh=0; hh<dim[0]; hh++ ){ delete[] x[hh]; }; delete[] x;
		}
		
		#pragma omp critical
		{
			for( int h=0; h<4; ++h ){
				countTables[h] += local_countTables[h];
				probTables[h] += local_probTables[h];
			}
			nO000 += local_nO000[0];
		}

		delete[] local_countTables;
		delete[] local_probTables;
		delete[] local_nO000;
	}
}


/*
 * 
 * 	CASE I x J ----------------UNDER CONSTRUCTION
 * 
 */



extern "C"
void recIxJ(int *reccall, struct StructConstants *Constants, int **x, int i, int j, int *rowsums, int *colsums, ULONGLONG *countTables, double *probTables, ULONGLONG *nO000){

// 	int rc=*reccall;
// 	*reccall=(++rc);
	int I=Constants->dim[0];
	int J=Constants->dim[1];
	int *dim=Constants->dim;
	int *margins=Constants->margins;
	int some_columns_left=0;
		
	if( j == J-1 ){ // reached last col
		int remain_row=margins[ i ] - rowsums[i];
		if( remain_row < 0 ){ return; }

		x[ i ][ j ]=remain_row;
		rowsums[i] += x[ i ][ j ];
		colsums[j] += x[ i ][ j ];

		
		if( i == I-1 ){ // reached last row
			update( x, dim, countTables, probTables, nO000
				, Constants->preCalcFact, Constants->O000, Constants->p0, Constants->diff_lfmarginsTotal_lfN );
		}else{
			++i; j=0; // go to beginning of next row
			some_columns_left=1;
		}
	}else{
		some_columns_left=1;
	}		
	
	if( some_columns_left ){
		int *rowsums0=(int*)malloc(I*sizeof(int));
		int *colsums0=(int*)malloc(J*sizeof(int));
		int **xx=new int*[ I ];
		for( int ii=0; ii<I; ++ii ){ 
			xx[ii]=new int[J];
			memcpy((void*)xx[ii], (void*)x[ii], J*sizeof(int));
			//for( int jj=0; jj<J; ++jj ){ xx[ii][jj]=x[ii][jj]; }
		}
		
		int remain_row=margins[i]-rowsums[i];
		int remain_col=margins[I+j]-colsums[j];
		int min_at_current_position=( (remain_row < remain_col) ? remain_row : remain_col );
		int h_from=0; 
		// if you are in last row and see insufficient colsum
		if( i+1==I ){ 
			h_from=remain_col; 
		}
		for( int h = h_from; h<=min_at_current_position; h++ ){
			xx[i][j]=x[i][j]; // re-init xx at relevant position
//  
			if( margins[ i ] - h < 0 ){ continue; } // skip if invalid margin (row)
			if( margins[ I+j ] - h < 0 ){ continue; } // skip if invalid margin (col)
			
			memcpy((void*)rowsums0, (void*)rowsums, I*sizeof(int));
			memcpy((void*)colsums0, (void*)colsums, J*sizeof(int));
			
			xx[i][j]=h; 
			rowsums0[i] += xx[i][j];
			colsums0[j] += xx[i][j]; 
			
			if( margins[ i ] - rowsums0[i] < 0 ){ continue; } // skip if invalid margin (row)
			if( margins[ I+j ] - colsums0[j] < 0 ){ continue; } // skip if invalid margin (col)
		
			// if you are in second last col: look forward to check if next call is promising
			if( (j+2 == J) && (margins[ I+j+1 ] < colsums0[j+1] + margins[ i ] - rowsums0[i] ) ){
				//printf("skipped crit slc\n"); 
				continue;
			}
			
			recIxJ( reccall, Constants, xx, i, j+1, rowsums0, colsums0, countTables, probTables, nO000 );
		}
		
		delete[] rowsums0;
		delete[] colsums0;
		for( int hh=0; hh<dim[0]; ++hh ){ delete[] xx[hh]; }; delete[] xx;
	}
}


extern "C"
void runIxJ(struct StructConstants *Constants, ULONGLONG *countTables, double *probTables, ULONGLONG *nO000, int *nthreads){

	int *margins=Constants->margins;
	int *dim=Constants->dim;
	int zero=0;
	int *reccall=&zero;
	int min00=( (margins[0] < margins[ dim[0] ]) ? margins[0] : margins[ dim[0] ] ); // get minimum range for x_00
	#pragma omp parallel shared(countTables, probTables, nO000) num_threads(*nthreads) if(*nthreads > 1)
	{
		ULONGLONG *local_nO000=new ULONGLONG[1];
		local_nO000[0]=0;
		ULONGLONG *local_countTables=new ULONGLONG[4];
		double *local_probTables=new double[4];
		for( int h=0; h<4; ++h ){ local_countTables[h]=0; local_probTables[h]=0.0; }
	
		#pragma omp for schedule(dynamic)
		for( int i=0; i<=min00; ++i ){
			int *rowsums=new int[ dim[0] ];
			int *colsums=new int[ dim[1] ];
			int **x=new int*[ dim[0] ];
			for( int hh=0; hh<dim[0]; ++hh ){
				x[hh]= new int[ dim[1] ]; rowsums[hh]=0;
				for( int gg=0; gg<dim[1]; ++gg ){ x[hh][gg]= -1; colsums[gg]=0; }
			}

			x[0][0]=i; 
			rowsums[0] += x[0][0]; colsums[0] += x[0][0];
			
			recIxJ( reccall, Constants, x, 0, 1, rowsums, colsums, local_countTables, local_probTables, local_nO000 );
			
			delete[] rowsums; delete[] colsums;
			for( int hh=0; hh<dim[0]; ++hh ){ delete[] x[hh]; }; delete[] x;
		}
		
		#pragma omp critical
		{
			for( int h=0; h<4; ++h ){
				countTables[h] += local_countTables[h];
				probTables[h] += local_probTables[h];
			}
			nO000 += local_nO000[0];
		}

		delete[] local_countTables;
		delete[] local_probTables;
		delete[] local_nO000;
	}
	//printf("reccall: %i\n", *reccall);
}










extern "C"
void hypergeom_IxJ(int *O000, int *N, int *margins, double *p0, double *n0, double *Prob, double *Freq, int *dim, int *nthreads){

	int NN=*N;
	int h=0;
	double *preCalcFact= new double[NN+1];
	double lfN=0.0;
	double lfmargins[2];
	double lfmarginsTotal=0.0;

	//results
	ULONGLONG zeroInt=0;
	ULONGLONG *countTables=new ULONGLONG[4];
	double *probTables=new double[4];
	ULONGLONG *nO000=&zeroInt;

	
	for( h=0; h<4; h++ ){ countTables[h]=0; probTables[h]=0.0; }

	preCalcFact[0]=0;
	for( h=1; h<= NN; ++h ){ preCalcFact[h]=preCalcFact[h-1]+log((double)h); }
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

	
	struct StructConstants *Constants;
	Constants=(struct StructConstants*)malloc(sizeof(struct StructConstants));
	Constants->preCalcFact=preCalcFact;
	Constants->diff_lfmarginsTotal_lfN=diff_lfmarginsTotal_lfN;
	Constants->p0=p0;
	Constants->O000=O000;
	Constants->dim=dim;
	Constants->margins=margins;
	
	
	if( dim[0]==2 && dim[1]==2 ){
		run2x2(Constants, countTables, probTables, nO000, nthreads);
	}else{
		if( dim[0]>=2 && dim[1]==2 ){
			runIx2(Constants, countTables, probTables, nO000, nthreads);
		}else{
			runIxJ(Constants, countTables, probTables, nO000, nthreads);
		}
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

	double ll=*nO000;
	*n0=ll;

	delete[] preCalcFact;
	delete[] countTables;
	delete[] probTables;
	//free(nO000);
	free(Constants);
}

