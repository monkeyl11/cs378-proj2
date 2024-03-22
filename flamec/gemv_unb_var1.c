
/* Copyright 2024 The University of Texas at Austin  
 
   For licensing information see
                  http://www.cs.utexas.edu/users/flame/license.html 

   Programmed by: Devangi Parikh
                  dnp@cs.utexas.edu
                                                                     */

#include "FLAME.h"

int gemv_unb_var1( FLA_Obj y, FLA_Obj A, FLA_Obj x )
{
  FLA_Obj AL,    AR,       A0,  a1,  A2;

  FLA_Obj xT,              x0,
          xB,              chi1,
                           x2;

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  FLA_Part_2x1( x,    &xT, 
                      &xB,            0, FLA_TOP );

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &a1, &A2,
                           1, FLA_RIGHT );

    FLA_Repart_2x1_to_3x1( xT,                &x0, 
                        /* ** */            /* **** */
                                              &chi1, 
                           xB,                &x2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

	FLA_Axpy( chi1, a1, y );

    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, a1, /**/ A2,
                              FLA_LEFT );

    FLA_Cont_with_3x1_to_2x1( &xT,                x0, 
                                                  chi1, 
                            /* ** */           /* **** */
                              &xB,                x2,     FLA_TOP );

  }

  return FLA_SUCCESS;
}




