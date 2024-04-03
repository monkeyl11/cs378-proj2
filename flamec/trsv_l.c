/* Copyright 2024 The University of Texas at Austin  
 
   For licensing information see
                  http://www.cs.utexas.edu/users/flame/license.html 

   Programmed by: Name of author
                  Email of author
                                                                     */

#include "FLAME.h"

int trsv_l_unb_var1( FLA_Obj L, FLA_Obj y )
{
  FLA_Obj LTL,   LTR,      L00,  l01,      L02, 
          LBL,   LBR,      l10t, lambda11, l12t,
                           L20,  l21,      L22;

  FLA_Obj yT,              y0,
          yB,              psi1,
                           y2;

  FLA_Part_2x2( L,    &LTL, &LTR,
                      &LBL, &LBR,     0, 0, FLA_TL );

  FLA_Part_2x1( y,    &yT, 
                      &yB,            0, FLA_TOP );

  while ( FLA_Obj_length( LTL ) < FLA_Obj_length( L ) ){

    FLA_Repart_2x2_to_3x3( LTL, /**/ LTR,       &L00,  /**/ &l01,      &L02,
                        /* ************* */   /* *************************** */
                                                &l10t, /**/ &lambda11, &l12t,
                           LBL, /**/ LBR,       &L20,  /**/ &l21,      &L22,
                           1, 1, FLA_BR );

    FLA_Repart_2x1_to_3x1( yT,                &y0, 
                        /* ** */            /* **** */
                                              &psi1, 
                           yB,                &y2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    FLA_Dots( FLA_MINUS_ONE, l10t, y0, FLA_ONE, psi1 );
    FLA_Inv_scal( lambda11, psi1 );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &LTL, /**/ &LTR,       L00,  l01,      /**/ L02,
                                                     l10t, lambda11, /**/ l12t,
                            /* ************** */  /* ************************* */
                              &LBL, /**/ &LBR,       L20,  l21,      /**/ L22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &yT,                y0, 
                                                  psi1, 
                            /* ** */           /* **** */
                              &yB,                y2,     FLA_TOP );

  }

  return FLA_SUCCESS;
}


/* Copyright 2024 The University of Texas at Austin  
 
   For licensing information see
                  http://www.cs.utexas.edu/users/flame/license.html 

   Programmed by: Name of author
                  Email of author
                                                                     */

int trsvl_l_blk( FLA_Obj L, FLA_Obj y, int nb_alg )
{
  FLA_Obj LTL,   LTR,      L00, L01, L02, 
          LBL,   LBR,      L10, L11, L12,
                           L20, L21, L22;

  FLA_Obj yT,              y0,
          yB,              y1,
                           y2;

  int b;

  FLA_Part_2x2( L,    &LTL, &LTR,
                      &LBL, &LBR,     0, 0, FLA_TL );

  FLA_Part_2x1( y,    &yT, 
                      &yB,            0, FLA_TOP );

  while ( FLA_Obj_length( LTL ) < FLA_Obj_length( L ) ){

    b = min( FLA_Obj_length( LBR ), nb_alg );

    FLA_Repart_2x2_to_3x3( LTL, /**/ LTR,       &L00, /**/ &L01, &L02,
                        /* ************* */   /* ******************** */
                                                &L10, /**/ &L11, &L12,
                           LBL, /**/ LBR,       &L20, /**/ &L21, &L22,
                           b, b, FLA_BR );

    FLA_Repart_2x1_to_3x1( yT,                &y0, 
                        /* ** */            /* ** */
                                              &y1, 
                           yB,                &y2,        b, FLA_BOTTOM );

    /*------------------------------------------------------------*/
    
    FLA_Gemv( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, L10, y0, FLA_ONE, y1 );
    if (FLA_Obj_length( LBR ) > nb_alg) {
      trsvl_l_blk( L11, y1, nb_alg);
    }
    else {
      trsv_l_unb_var1(L11, y1);      
    }
    
    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &LTL, /**/ &LTR,       L00, L01, /**/ L02,
                                                     L10, L11, /**/ L12,
                            /* ************** */  /* ****************** */
                              &LBL, /**/ &LBR,       L20, L21, /**/ L22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &yT,                y0, 
                                                  y1, 
                            /* ** */           /* ** */
                              &yB,                y2,     FLA_TOP );

  }

  return FLA_SUCCESS;
}


