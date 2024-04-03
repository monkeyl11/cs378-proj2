
/* Copyright 2024 The University of Texas at Austin  
 
   For licensing information see
                  http://www.cs.utexas.edu/users/flame/license.html 

   Programmed by: Name of author
                  Email of author
                                                                     */

#include "FLAME.h"

int Trsv_l_unb_var1( FLA_Obj L, FLA_Obj x, FLA_Obj y )
{
  FLA_Obj LTL,   LTR,      L00,  l01,      L02, 
          LBL,   LBR,      l10t, lambda11, l12t,
                           L20,  l21,      L22;

  FLA_Obj xT,              x0,
          xB,              chi1,
                           x2;

  FLA_Obj yT,              y0,
          yB,              psi1,
                           y2;

  FLA_Part_2x2( L,    &LTL, &LTR,
                      &LBL, &LBR,     0, 0, FLA_TL );

  FLA_Part_2x1( x,    &xT, 
                      &xB,            0, FLA_TOP );

  FLA_Part_2x1( y,    &yT, 
                      &yB,            0, FLA_TOP );

  while ( FLA_Obj_length( LTL ) < FLA_Obj_length( L ) ){

    FLA_Repart_2x2_to_3x3( LTL, /**/ LTR,       &L00,  /**/ &l01,      &L02,
                        /* ************* */   /* *************************** */
                                                &l10t, /**/ &lambda11, &l12t,
                           LBL, /**/ LBR,       &L20,  /**/ &l21,      &L22,
                           1, 1, FLA_BR );

    FLA_Repart_2x1_to_3x1( xT,                &x0, 
                        /* ** */            /* **** */
                                              &chi1, 
                           xB,                &x2,        1, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( yT,                &y0, 
                        /* ** */            /* **** */
                                              &psi1, 
                           yB,                &y2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    /*                       update line 1                        */

    FLA_Obj temp;
    FLA_Obj_create( FLA_DOUBLE, 1, 1, 1, 1, &temp );
    FLA_Copy(psi1, temp);
    FLA_Dots( FLA_MINUS_ONE, l10t, x0, FLA_ONE, temp );
    FLA_Inv_scal( lambda11, temp );
    FLA_Copy(temp, chi1);

    FLA_Obj_free( &temp );
    /*                       update line n                        */

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &LTL, /**/ &LTR,       L00,  l01,      /**/ L02,
                                                     l10t, lambda11, /**/ l12t,
                            /* ************** */  /* ************************* */
                              &LBL, /**/ &LBR,       L20,  l21,      /**/ L22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &xT,                x0, 
                                                  chi1, 
                            /* ** */           /* **** */
                              &xB,                x2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &yT,                y0, 
                                                  psi1, 
                            /* ** */           /* **** */
                              &yB,                y2,     FLA_TOP );

  }

  return FLA_SUCCESS;
}


