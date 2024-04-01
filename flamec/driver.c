#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include "FLAME.h"

/* Various constants that control what gets timed */

#define TRUE 1
#define FALSE 0

int symm_ll_unb_var1( FLA_Obj, FLA_Obj, FLA_Obj );
int syr2k_ln_unb_var1( FLA_Obj, FLA_Obj, FLA_Obj );
int Trsv_l_unb_var1( FLA_Obj, FLA_Obj, FLA_Obj );

void cat1(int write_out);
void cat2(int write_out);
void cat3(int write_out);

int main(int argc, char *argv[])
{
  //cat1(1);
  cat2(1);
  //cat3(1);
  exit(0);
}

void cat1(int write_out) {
  int n, nfirst, nlast, ninc, i, irep, nrepeats;
  int m = 123;

  double
    dtime, dtime_best, 
    diff;

  dtime_best = 0.0;

  FLA_Obj A, B, C; 

  FLA_Obj Cref, Cold;
  
  /* Initialize FLAME. */
  FLA_Init( );

  /* Every time trial is repeated "repeat" times and the fastest run in recorded */
  printf( "%% number of repeats:" );
  scanf( "%d", &nrepeats );
  printf( "%% %d\n", nrepeats );

  /* Timing trials for matrix sizes n=nfirst to nlast in increments 
     of ninc will be performed. */
  printf( "%% enter nfirst, nlast, ninc:" );
  scanf( "%d%d%d", &nfirst, &nlast, &ninc );
  printf( "%% %d %d %d \n", nfirst, nlast, ninc );
  fflush( stdout );

  i = 1;
  FILE *out = fopen("output.txt", "w");
  if (write_out && out == NULL) {
    printf("Failed to open output file\n");
    exit(1);
  }
  int dim; double ref_time; double our_time;
  for ( n=nfirst; n<= nlast; n+=ninc ){
    dim = n;
    /* Allocate space for the matrices and vectors */
    FLA_Obj_create( FLA_DOUBLE, n, m, 1, n, &C );
    FLA_Obj_create( FLA_DOUBLE, n, m, 1, n, &Cref );
    FLA_Obj_create( FLA_DOUBLE, n, m, 1, n, &Cold );
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &A );
    FLA_Obj_create( FLA_DOUBLE, n, m, 1, n, &B );

    /* Generate random matrix A, and vectors x, and y */
    FLA_Random_symm_matrix( FLA_LOWER_TRIANGULAR, A );
    FLA_Random_matrix( Cold );
    FLA_Random_matrix( B );

    for ( irep=0; irep<nrepeats; irep++ ) {
    /* Time reference implementation (from libflame) */
      FLA_Copy( Cold, Cref );
    
      /* start clock */
      dtime = FLA_Clock();
    
      //FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A, x, FLA_ONE, yref );
      FLA_Gemm(FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, B, FLA_ONE, Cref);
 
      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    printf( "data_ref( %d, 1:2 ) = [ %d %le ];\n", i, n, dtime_best );
    ref_time = dtime_best;
    fflush( stdout );
    FLA_Triangularize( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, A );

    /* Time your unblocked Variant 1 */

    for ( irep=0; irep<nrepeats; irep++ ){
      /* Copy vector yold to y */
      FLA_Copy( Cold, C );
    
      /* start clock */
      dtime = FLA_Clock();
 
      /* Comment out the below call and call your routine instead */
      symm_ll_unb_var1( A, B, C );
      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    diff = FLA_Max_elemwise_diff( C, Cref );
  
    printf( "data_unb_var1( %d, 1:3 ) = [ %d %le %le];\n", i, n,
	    dtime_best, diff  );

    our_time = dtime_best;
    fflush( stdout );

    FLA_Obj_free( &C );
    FLA_Obj_free( &Cold );
    FLA_Obj_free( &Cref );
    FLA_Obj_free( &A );
    FLA_Obj_free( &B );

    i++;
    if (write_out) {
      fprintf(out, "%d %f %f\n", dim, our_time, ref_time);
    }
  }

  FLA_Finalize( );
}

void cat2(int write_out) {
  int n, nfirst, nlast, ninc, i, irep, nrepeats;
  int m = 123;

  double
    dtime, dtime_best, 
    diff;

  dtime_best = 0.0;

  FLA_Obj A, B, C; 

  FLA_Obj Cref, Cold;
  
  /* Initialize FLAME. */
  FLA_Init( );

  /* Every time trial is repeated "repeat" times and the fastest run in recorded */
  printf( "%% number of repeats:" );
  scanf( "%d", &nrepeats );
  printf( "%% %d\n", nrepeats );

  /* Timing trials for matrix sizes n=nfirst to nlast in increments 
     of ninc will be performed. */
  printf( "%% enter nfirst, nlast, ninc:" );
  scanf( "%d%d%d", &nfirst, &nlast, &ninc );
  printf( "%% %d %d %d \n", nfirst, nlast, ninc );
  fflush( stdout );

  i = 1;
  FILE *out = fopen("output.txt", "w");
  if (write_out && out == NULL) {
    printf("Failed to open output file\n");
    exit(1);
  }
  int dim; double ref_time; double our_time;
  for ( n=nfirst; n<= nlast; n+=ninc ){
    dim = n;
    /* Allocate space for the matrices and vectors */
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &C );
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Cref );
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &Cold );
    FLA_Obj_create( FLA_DOUBLE, n, m, 1, n, &A );
    FLA_Obj_create( FLA_DOUBLE, n, m, 1, n, &B );

    /* Generate random matrix A, and vectors x, and y */
    FLA_Random_symm_matrix( FLA_LOWER_TRIANGULAR, Cold );
    FLA_Random_matrix( A );
    FLA_Random_matrix( B );

    for ( irep=0; irep<nrepeats; irep++ ) {
    /* Time reference implementation (from libflame) */
      FLA_Copy( Cold, Cref );
    
      /* start clock */
      dtime = FLA_Clock();
    
      //FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A, x, FLA_ONE, yref );
      FLA_Gemm(FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, A, B, FLA_ONE, Cref);
      FLA_Gemm(FLA_NO_TRANSPOSE, FLA_TRANSPOSE, FLA_ONE, B, A, FLA_ONE, Cref);
 
      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    printf( "data_ref( %d, 1:2 ) = [ %d %le ];\n", i, n, dtime_best );
    ref_time = dtime_best;
    fflush( stdout );
    FLA_Triangularize( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, Cold );
    FLA_Triangularize( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, Cref );

    /* Time your unblocked Variant 1 */

    for ( irep=0; irep<nrepeats; irep++ ){
      /* Copy vector yold to y */
      FLA_Copy( Cold, C );
    
      /* start clock */
      dtime = FLA_Clock();
 
      /* Comment out the below call and call your routine instead */
      syr2k_ln_unb_var1( A, B, C );

      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    diff = FLA_Max_elemwise_diff( C, Cref );
  
    printf( "data_unb_var1( %d, 1:3 ) = [ %d %le %le];\n", i, n,
	    dtime_best, diff  );

    our_time = dtime_best;
    fflush( stdout );

    FLA_Obj_free( &C );
    FLA_Obj_free( &Cold );
    FLA_Obj_free( &Cref );
    FLA_Obj_free( &A );
    FLA_Obj_free( &B );

    i++;
    if (write_out) {
      fprintf(out, "%d %f %f\n", dim, our_time, ref_time);
    }
  }
  FLA_Finalize( );
}

void cat3(int write_out) {
  int n, nfirst, nlast, ninc, i, irep, nrepeats;

  double
    dtime, dtime_best, 
    diff;

  dtime_best = 0.0;

  FLA_Obj L, x, y; 

  FLA_Obj xref, xold;
  
  /* Initialize FLAME. */
  FLA_Init( );

  /* Every time trial is repeated "repeat" times and the fastest run in recorded */
  printf( "%% number of repeats:" );
  scanf( "%d", &nrepeats );
  printf( "%% %d\n", nrepeats );

  /* Timing trials for matrix sizes n=nfirst to nlast in increments 
     of ninc will be performed. */
  printf( "%% enter nfirst, nlast, ninc:" );
  scanf( "%d%d%d", &nfirst, &nlast, &ninc );
  printf( "%% %d %d %d \n", nfirst, nlast, ninc );
  fflush( stdout );

  i = 1;
  FILE *out = fopen("output.txt", "w");
  if (write_out && out == NULL) {
    printf("Failed to open output file\n");
    exit(1);
  }
  int dim; double ref_time; double our_time;
  for ( n=nfirst; n<= nlast; n+=ninc ){
    dim = n;
    /* Allocate space for the matrices and vectors */
    FLA_Obj_create( FLA_DOUBLE, n, 1, 1, 1, &x );
    FLA_Obj_create( FLA_DOUBLE, n, 1, 1, 1, &y );
    FLA_Obj_create( FLA_DOUBLE, n, 1, 1, 1, &xref );
    FLA_Obj_create( FLA_DOUBLE, n, 1, 1, 1, &xold );
    FLA_Obj_create( FLA_DOUBLE, n, n, 1, n, &L );

    /* Generate random matrix A, and vectors x, and y */
    FLA_Random_symm_matrix( FLA_LOWER_TRIANGULAR, L );
    FLA_Triangularize( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, L );
    FLA_Random_matrix( xold );
    FLA_Random_matrix( y );

    for ( irep=0; irep<nrepeats; irep++ ) {
    /* Time reference implementation (from libflame) */
      FLA_Copy( y, xref );
    
      /* start clock */
      dtime = FLA_Clock();
    
      //FLA_Gemv( FLA_NO_TRANSPOSE, FLA_ONE, A, x, FLA_ONE, yref );
      FLA_Trsv(FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG, L, xref);

 
      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    printf( "data_ref( %d, 1:2 ) = [ %d %le ];\n", i, n, dtime_best );
    ref_time = dtime_best;
    fflush( stdout );

    /* Time your unblocked Variant 1 */

    for ( irep=0; irep<nrepeats; irep++ ){
      /* Copy vector yold to y */
      FLA_Copy( xold, x );
    
      /* start clock */
      dtime = FLA_Clock();
 
      /* Comment out the below call and call your routine instead */
      Trsv_l_unb_var1( L, x, y );

      /* stop clock */
      dtime = FLA_Clock() - dtime;
    
      if ( irep == 0 ) 
	dtime_best = dtime;
      else
	dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
    }

    diff = FLA_Max_elemwise_diff( xref, x );
    printf( "data_unb_var1( %d, 1:3 ) = [ %d %le %le];\n", i, n,
	    dtime_best, diff  );
    our_time = dtime_best;
    fflush( stdout );

    FLA_Obj_free( &L );
    FLA_Obj_free( &xold );
    FLA_Obj_free( &xref );
    FLA_Obj_free( &x );
    FLA_Obj_free( &y );

    i++;
    if (write_out) {
      fprintf(out, "%d %f %f\n", dim, our_time, ref_time);
    }
  }
  FLA_Finalize( );
}
