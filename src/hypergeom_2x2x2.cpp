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


extern "C"
void hypergeom_2x2x2(int *O000, int *N, int *marg, int *margins, double *p0, double *n0, double *Prob, double *Freq, int *nthreads){

	int NN=*N;

	int h=0;
	double *preCalcFact=new double[NN+1];

	double lfN=0.0;
	double lfmargins[3];
	double lfmarginsTotal=0.0;


	//results
	double probTables[4];
	ULONGLONG countTables[4];
	ULONGLONG nO000=0;
	for( h=0; h<4; h++ ){ countTables[h]=0; probTables[h]=0.0; }

	preCalcFact[0]=0;
	for( h=1; h<= NN; h++){
		preCalcFact[h]=preCalcFact[h-1]+log((double)h);
	}
	lfN=2.0*preCalcFact[NN];
	lfmargins[0]=preCalcFact[ margins[0] ] + preCalcFact[ margins[1] ];
	lfmargins[1]=preCalcFact[ margins[2] ] + preCalcFact[ margins[3] ];
	lfmargins[2]=preCalcFact[ margins[4] ] + preCalcFact[ margins[5] ];
	lfmarginsTotal=lfmargins[0]+lfmargins[1]+lfmargins[2];

	double diff_lfmarginsTotal_lfN=lfmarginsTotal-lfN;
	int N_minus_marg2=NN - marg[2];

	int i1=0;

	#ifdef _OPENMP
	if(*nthreads <= 1){ *nthreads=1; }else{ *nthreads=(*nthreads < omp_get_max_threads()) ? (*nthreads) : (omp_get_max_threads()); }
	#endif

	#pragma omp parallel shared(countTables, probTables, nO000) num_threads(*nthreads)
	{

		ULONGLONG local_nO000=0;
		ULONGLONG local_countTables[4];
		double local_probTables[4];
		for( h=0; h<4; ++h ){ local_countTables[h]=0; local_probTables[h]=0.0; }

		#pragma omp for schedule(dynamic)
		for( i1=0; i1<=marg[0]; ++i1 ){

			int x[2][2][2];
			x[0][0][0]=i1;

			double lm1=0.0;
			lm1 += preCalcFact[ x[0][0][0] ];

			int marg1_plus_marg2_minus_N_x000=marg[1]+marg[2] - NN - x[0][0][0];

			int i2=0;
			int minSecondLoop=(marg[0] - x[0][0][0] );
			for( i2=0; i2<=minSecondLoop; ++i2 ){

				x[0][0][1]=i2;

				double lm2=preCalcFact[ x[0][0][1] ];

				int marg0_minus_x000_minus_x001=marg[0] - x[0][0][0] - x[0][0][1];
				int marg1_minus_x000_minus_x001=marg[1] - x[0][0][0] - x[0][0][1];

				int i3=0;
				for( i3=0; i3<=(marg0_minus_x000_minus_x001); ++i3 ){
					x[0][1][0]=i3;
					x[0][1][1]=marg0_minus_x000_minus_x001 - x[0][1][0];

					int marg2_minus_x000_minus_x010=marg[2]- x[0][0][0] - x[0][1][0];

					int aux=marg1_plus_marg2_minus_N_x000 + x[0][1][1];
					int m1=0; if( aux>0 ){m1=aux;}
					int m2=marg1_minus_x000_minus_x001;
					aux=marg2_minus_x000_minus_x010;
					if(m2>aux){m2=aux;}
					if( m1>m2 ){continue;}

					double lm3=preCalcFact[ x[0][1][0] ] + preCalcFact[ x[0][1][1] ];

					int i4=0;
					for( i4 = (m1+1); i4 <= (m2+1); ++i4 ){
						x[1][0][0] = i4-1;
						x[1][0][1] = marg1_minus_x000_minus_x001 - x[1][0][0];
						x[1][1][0] = marg2_minus_x000_minus_x010 - x[1][0][0];
						x[1][1][1] = N_minus_marg2 - x[0][0][1] - x[0][1][1] - x[1][0][1];

						double lm4=preCalcFact[ x[1][0][0] ]
							+ preCalcFact[ x[1][0][1] ]
							+ preCalcFact[ x[1][1][0] ]
							+ preCalcFact[ x[1][1][1] ];

						double lm=lm1+lm2+lm3+lm4;

						double prob= (diff_lfmarginsTotal_lfN - lm) < -708 ? 0.0 : exp(diff_lfmarginsTotal_lfN - lm);
						local_probTables[0] += prob; ++local_countTables[0];
						if( x[0][0][0] == *O000 ){ ++local_nO000; }
						if( x[0][0][0] <= *O000 ){ local_probTables[1] += prob; ++local_countTables[1]; } //less
						if( x[0][0][0] >= *O000 ){ local_probTables[2] += prob; ++local_countTables[2]; } //greater
						if( prob <= *p0 ){ local_probTables[3] += prob; ++local_countTables[3]; } //two sided, minimum likelihood

					} //end loop i4
				} //end loop i3
			} //end loop i2
		} //end loop i1

		#pragma omp critical
		{
			for( h=0; h<4; h++ ){
				countTables[h] += local_countTables[h];
				probTables[h] += local_probTables[h];
			}
			nO000 += local_nO000;
		}
	}

	for( h=0; h<4; h++ ){
		Freq[h] += countTables[h];
		Prob[h] += probTables[h];
	}
	double min=Prob[1];
	ULONGLONG min_count=Freq[1];
	if( min>Prob[2] ){ min=Prob[2];min_count=Freq[2]; }
	Prob[4]=2.0*min; if(Prob[4]>1.0){Prob[4]=1.0;}
	Freq[4]=2.0*min_count;

	*n0=nO000;

	delete [] preCalcFact;
}
