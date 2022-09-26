/*! \file GaIn.h
    \brief function declarations for linking with libGaIn
*/
#ifndef LIB_GAIN_H
#define LIB_GAIN_H

// definition of fortran interface of the library
#include "GaIn.FC.h"

// definition of fortran interface of the library
#define CC_Overlap_CC                   GAIN_FC_FUNC(cc_overlap_cc,CC_OVERLAP_CC)
#define YY_Overlap_YY                   GAIN_FC_FUNC(yy_overlap_yy,YY_OVERLAP_YY)
#define CC_Overlap_C                    GAIN_FC_FUNC(cc_overlap_c,CC_OVERLAP_C)
#define CC_Overlap_Y                    GAIN_FC_FUNC(cc_overlap_y,CC_OVERLAP_Y)
#define CY_Overlap_C                    GAIN_FC_FUNC(cy_overlap_c,CY_OVERLAP_C)
#define YC_Overlap_C                    GAIN_FC_FUNC(yc_overlap_c,YC_OVERLAP_C)
#define YY_Overlap_Y                    GAIN_FC_FUNC(yy_overlap_y,YY_OVERLAP_Y)
#define YY_Overlap_C                    GAIN_FC_FUNC(yy_overlap_c,YY_OVERLAP_C)
#define YC_Overlap_Y                    GAIN_FC_FUNC(yc_overlap_y,YC_OVERLAP_Y)
#define CY_Overlap_Y                    GAIN_FC_FUNC(cy_overlap_y,CY_OVERLAP_Y)
#define C_R2Moment_C                    GAIN_FC_FUNC(c_r2moment_c,C_R2MOMENT_C)
#define CC_Coulomb_Ion                  GAIN_FC_FUNC(cc_coulomb_ion,CC_COULOMB_ION)
#define CC_Coulomb_CC                   GAIN_FC_FUNC(cc_coulomb_cc,CC_COULOMB_CC)
#define CC_Coulomb_C                    GAIN_FC_FUNC(cc_coulomb_c,CC_COULOMB_C)
#define CC_Coulomb_Y                    GAIN_FC_FUNC(cc_coulomb_y,CC_COULOMB_Y)
#define CY_Coulomb_C                    GAIN_FC_FUNC(cy_coulomb_c,CY_COULOMB_C)
#define YC_Coulomb_C                    GAIN_FC_FUNC(yc_coulomb_c,YC_COULOMB_C)
#define YY_Coulomb_Y                    GAIN_FC_FUNC(yy_coulomb_y,YY_COULOMB_Y)
#define YY_Coulomb_C                    GAIN_FC_FUNC(yy_coulomb_c,YY_COULOMB_C)
#define YC_Coulomb_Y                    GAIN_FC_FUNC(yc_coulomb_y,YC_COULOMB_Y)
#define CY_Coulomb_Y                    GAIN_FC_FUNC(cy_coulomb_y,CY_COULOMB_Y)
#define C_Coulomb_C                     GAIN_FC_FUNC(c_coulomb_c,C_COULOMB_C)
#define C_Coulomb_Y                     GAIN_FC_FUNC(c_coulomb_y,C_COULOMB_Y)
#define Y_Coulomb_C                     GAIN_FC_FUNC(y_coulomb_c,Y_COULOMB_C)
#define Y_Coulomb_Y                     GAIN_FC_FUNC(y_coulomb_y,Y_COULOMB_Y)
#define C_Overlap_C                     GAIN_FC_FUNC(c_overlap_c,C_OVERLAP_C)
#define Y_Overlap_Y                     GAIN_FC_FUNC(y_overlap_y,Y_OVERLAP_Y)
#define Y_Overlap_C                     GAIN_FC_FUNC(y_overlap_c,Y_OVERLAP_C)
#define C_Overlap_Y                     GAIN_FC_FUNC(c_overlap_y,C_OVERLAP_Y)
#define Y_R2Moment_Y                    GAIN_FC_FUNC(y_r2moment_y,Y_R2MOMENT_Y)
#define Y_Laplacian_Y                   GAIN_FC_FUNC(y_laplacian_y,Y_LAPLACIAN_Y)
#define Y_Laplacian_C                   GAIN_FC_FUNC(y_laplacian_c,Y_LAPLACIAN_C)
#define C_Laplacian_Y                   GAIN_FC_FUNC(c_laplacian_y,C_LAPLACIAN_Y)
#define C_Laplacian_C                   GAIN_FC_FUNC(c_laplacian_c,C_LAPLACIAN_C)
#define YY_Coulomb_Ion                  GAIN_FC_FUNC(yy_coulomb_ion,YY_COULOMB_ION)
#define YY_Coulomb_YY                   GAIN_FC_FUNC(yy_coulomb_yy,YY_COULOMB_YY)
#define Y_Value                         GAIN_FC_FUNC(y_value,Y_VALUE)
#define Overlap_upper_bound             GAIN_FC_FUNC(overlap_upper_bound,OVERLAP_UPPER_BOUND)
#define C_to_D                          GAIN_FC_FUNC(c_to_d,C_TO_D)
#define Y_to_D                          GAIN_FC_FUNC(y_to_d,Y_TO_D)
#define CC_to_D                         GAIN_FC_FUNC(cc_to_d,CC_TO_D)
#define YC_to_D                         GAIN_FC_FUNC(yc_to_d,YC_TO_D)
#define CY_to_D                         GAIN_FC_FUNC(cy_to_d,CY_TO_D)
#define YY_to_D                         GAIN_FC_FUNC(yy_to_d,YY_TO_D)
#define D_Coulomb_D                     GAIN_FC_FUNC(d_coulomb_d,D_COULOMB_D)
#define D_Coulomb_Y                     GAIN_FC_FUNC(d_coulomb_y,D_COULOMB_Y)
#define D_Coulomb_C                     GAIN_FC_FUNC(d_coulomb_c,D_COULOMB_C)
#define D_Coulomb_Y_shell               GAIN_FC_FUNC(d_coulomb_y_shell,D_COULOMB_Y_SHELL)
#define D_Coulomb_C_shell               GAIN_FC_FUNC(d_coulomb_c_shell,D_COULOMB_C_SHELL)
#define D_Coulomb_Y_shells              GAIN_FC_FUNC(d_coulomb_y_shells,D_COULOMB_Y_SHELLS)
#define D_Coulomb_C_shells              GAIN_FC_FUNC(d_coulomb_c_shells,D_COULOMB_C_SHELLS)
#define Y_Rotation_Matrix               GAIN_FC_FUNC(y_rotation_matrix,Y_ROTATION_MATRIX)
#define R_Rotation_Matrix               GAIN_FC_FUNC(r_rotation_matrix,R_ROTATION_MATRIX)
#define Potential_from_C                GAIN_FC_FUNC(potential_from_c,POTENTIAL_FROM_C)
#define Potential_from_Y                GAIN_FC_FUNC(potential_from_y,POTENTIAL_FROM_Y)
#define Field_from_C                    GAIN_FC_FUNC(field_from_c,FIELD_FROM_C)
#define Field_from_Y                    GAIN_FC_FUNC(field_from_y,FIELD_FROM_Y)
#define Electrostatics_from_C           GAIN_FC_FUNC(electrostatics_from_c,ELECTROSTATICS_FROM_C)
#define Electrostatics_from_Y           GAIN_FC_FUNC(electrostatics_from_y,ELECTROSTATICS_FROM_Y)
#define Cumulative_electrostatics_on_D  GAIN_FC_FUNC(cumulative_electrostatics_on_d,CUMULATIVE_ELECTROSTATICS_ON_D)
#define Cumulative_electrostatics_on_C  GAIN_FC_FUNC(cumulative_electrostatics_on_c,CUMULATIVE_ELECTROSTATICS_ON_C)
#define Cumulative_electrostatics_on_Y  GAIN_FC_FUNC(cumulative_electrostatics_on_y,CUMULATIVE_ELECTROSTATICS_ON_Y)
#define Cumulative_electrostatics_on_YY GAIN_FC_FUNC(cumulative_electrostatics_on_yy,CUMULATIVE_ELECTROSTATICS_ON_YY)


#ifdef __cplusplus
extern "C"
  {
#endif
      double C_Overlap_C   (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2);
      double C_Overlap_Y   (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *l2 , const int *m2);
      double Y_Overlap_C   (const double *a1, const double *r1, const int *l1 , const int *m1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2);
      double Y_Overlap_Y   (const double *a1, const double *r1, const int *l1 , const int *m1, 
                            const double *a2, const double *r2, const int *l2 , const int *m2);
      double CC_Overlap_C  (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2,
                            const double *a3, const double *r3, const int *nx3, const int *ny3, const int *nz3);
      double CC_Overlap_Y  (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2,
                            const double *a3, const double *r3, const int *l3 , const int *m3);
      double CY_Overlap_C  (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *l2 , const int *m2,
                            const double *a3, const double *r3, const int *nx3, const int *ny3, const int *nz3);
      double YC_Overlap_C  (const double *a1, const double *r1, const int *l1 , const int *m1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2,
                            const double *a3, const double *r3, const int *nx3, const int *ny3, const int *nz3);
      double YY_Overlap_C  (const double *a1, const double *r1, const int *l1, const int *m1, 
                            const double *a2, const double *r2, const int *l2, const int *m2,
                            const double *a3, const double *r3, const int *nx3, const int *ny3, const int *nz3);
      double YC_Overlap_Y  (const double *a1, const double *r1, const int *l1, const int *m1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2,
                            const double *a3, const double *r3, const int *l3, const int *m3);
      double CY_Overlap_Y  (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *l2, const int *m2,
                            const double *a3, const double *r3, const int *l3, const int *m3);
      double YY_Overlap_Y  (const double *a1, const double *r1, const int *l1, const int *m1, 
                            const double *a2, const double *r2, const int *l2, const int *m2,
                            const double *a3, const double *r3, const int *l3, const int *m3);
      double CC_Overlap_CC (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2,
                            const double *a3, const double *r3, const int *nx3, const int *ny3, const int *nz3,
                            const double *a4, const double *r4, const int *nx4, const int *ny4, const int *nz4);
      double YY_Overlap_YY (const double *a1, const double *r1, const int *l1, const int *m1, 
                            const double *a2, const double *r2, const int *l2, const int *m2,
                            const double *a3, const double *r3, const int *l3, const int *m3,
                            const double *a4, const double *r4, const int *l4, const int *m4);
      double C_Laplacian_C (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2);
      double C_Laplacian_Y (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *l2 , const int *m2);
      double Y_Laplacian_C (const double *a1, const double *r1, const int *l1 , const int *m1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2);
      double Y_Laplacian_Y (const double *a1, const double *r1, const int *l1, const int *m1, 
                            const double *a2, const double *r2, const int *l2, const int *m2);
      double C_Coulomb_C   (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2);
      double C_Coulomb_Y   (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *l2 , const int *m2);
      double Y_Coulomb_C   (const double *a1, const double *r1, const int *l1 , const int *m1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2);
      double Y_Coulomb_Y   (const double *a1, const double *r1, const int *l1, const int *m1, 
                            const double *a2, const double *r2, const int *l2, const int *m2);
      double CC_Coulomb_C  (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2, 
                            const double *a3, const double *r3, const int *nx3, const int *ny3, const int *nz3);
      double CC_Coulomb_Y  (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2, 
                            const double *a3, const double *r3, const int *l3 , const int *m3);
      double CY_Coulomb_C  (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *l2 , const int *m2, 
                            const double *a3, const double *r3, const int *nx3, const int *ny3, const int *nz3);
      double YC_Coulomb_C  (const double *a1, const double *r1, const int *l1 , const int *m1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2, 
                            const double *a3, const double *r3, const int *nx3, const int *ny3, const int *nz3);
      double YY_Coulomb_C  (const double *a1, const double *r1, const int *l1, const int *m1, 
                            const double *a2, const double *r2, const int *l2, const int *m2, 
                            const double *a3, const double *r3, const int *nx3, const int *ny3, const int *nz3);
      double YC_Coulomb_Y  (const double *a1, const double *r1, const int *l1, const int *m1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2, 
                            const double *a3, const double *r3, const int *l3, const int *m3);
      double CY_Coulomb_Y  (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *l2, const int *m2, 
                            const double *a3, const double *r3, const int *l3, const int *m3);
      double YY_Coulomb_Y  (const double *a1, const double *r1, const int *l1, const int *m1, 
                            const double *a2, const double *r2, const int *l2, const int *m2, 
                            const double *a3, const double *r3, const int *l3, const int *m3);
      double CC_Coulomb_CC (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1,
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2, 
                            const double *a3, const double *r3, const int *nx3, const int *ny3, const int *nz3, 
                            const double *a4, const double *r4, const int *nx4, const int *ny4, const int *nz4);
      double YY_Coulomb_YY (const double *a1, const double *r1, const int *l1, const int *m1,
                            const double *a2, const double *r2, const int *l2, const int *m2, 
                            const double *a3, const double *r3, const int *l3, const int *m3, 
                            const double *a4, const double *r4, const int *l4, const int *m4);
      double CC_Coulomb_Ion(const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2, 
                            const double *rion);
      double YY_Coulomb_Ion(const double *a1, const double *r1, const int *l1, const int *m1, 
                            const double *a2, const double *r2, const int *l2, const int *m2, 
                            const double *rion);
      double C_R2Moment_C  (const double *a1, const double *r1, const int *nx1, const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2,
                            const int *nr2);
      double Y_R2Moment_Y  (const double *a1, const double *r1, const int *l1, const int *m1, 
                            const double *a2, const double *r2, const int *l2, const int *m2,
                            const int *nr2);
      double Y_Value       (const double *a1, const int *l1, const int *m1, const double *r);
      double Overlap_upper_bound(const double *a1, const double *r1, const int *l1,
                                 const double *a2, const double *r2, const int *l2);
      void   C_to_D        (const double *a1, const double *r1, const int *nx1  , const int *ny1, const int *nz1, 
                                  double *aD,   //< exponent (double)    for the derivatives basis
                                  double *rD,   //< center   (double[3]) for the derivatives basis
                                  double *cD,   //< derivative coefficient array (packed, see iD for corresponding derivative indices). May contain up to 455 coeffs.
                                  int    *iD,   //< derivative index array iD -> ndx,ndy,ndz. May contain up to 455 coeffs.
                                  int    *nD);  //< number of derivatives
      void   Y_to_D        (const double *a1, const double *r1, const int *l1, const int *m1, 
                                  double *aD,   //< exponent (double)    for the derivatives basis
                                  double *rD,   //< center   (double[3]) for the derivatives basis
                                  double *cD,   //< derivative coefficient array (packed, see iD for corresponding derivative indices). May contain up to 455 coeffs.
                                  int    *iD,   //< derivative index array iD -> ndx,ndy,ndz. May contain up to 455 coeffs.
                                  int    *nD);  //< number of derivatives
      void   CC_to_D       (const double *a1, const double *r1, const int *nx1  , const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2, 
                                  double *aD,   //< exponent (double)    for the derivatives basis
                                  double *rD,   //< center   (double[3]) for the derivatives basis
                                  double *cD,   //< derivative coefficient array (packed, see iD for corresponding derivative indices). May contain up to 455 coeffs.
                                  int    *iD,   //< derivative index array iD -> ndx,ndy,ndz. May contain up to 455 coeffs.
                                  int    *nD);  //< number of derivatives
      void   YC_to_D       (const double *a1, const double *r1, const int *l1, const int *m1, 
                            const double *a2, const double *r2, const int *nx2, const int *ny2, const int *nz2, 
                                  double *aD,   //< exponent (double)    for the derivatives basis
                                  double *rD,   //< center   (double[3]) for the derivatives basis
                                  double *cD,   //< derivative coefficient array (packed, see iD for corresponding derivative indices). May contain up to 455 coeffs.
                                  int    *iD,   //< derivative index array iD -> ndx,ndy,ndz. May contain up to 455 coeffs.
                                  int    *nD);  //< number of derivatives
      void   CY_to_D       (const double *a1, const double *r1, const int *nx1  , const int *ny1, const int *nz1, 
                            const double *a2, const double *r2, const int *l2, const int *m2, 
                                  double *aD,   //< exponent (double)    for the derivatives basis
                                  double *rD,   //< center   (double[3]) for the derivatives basis
                                  double *cD,   //< derivative coefficient array (packed, see iD for corresponding derivative indices). May contain up to 455 coeffs.
                                  int    *iD,   //< derivative index array iD -> ndx,ndy,ndz. May contain up to 455 coeffs.
                                  int    *nD);  //< number of derivatives
      void   YY_to_D       (const double *a1, const double *r1, const int *l1, const int *m1, 
                            const double *a2, const double *r2, const int *l2, const int *m2, 
                                  double *aD,   //< exponent (double)    for the derivatives basis
                                  double *rD,   //< center   (double[3]) for the derivatives basis
                                  double *cD,   //< derivative coefficient array (packed, see iD for corresponding derivative indices). May contain up to 455 coeffs.
                                  int    *iD,   //< derivative index array iD -> ndx,ndy,ndz. May contain up to 455 coeffs.
                                  int    *nD);  //< number of derivatives
      double D_Coulomb_D   (const double *aD1, const double *rD1, const double *cD1, const int *iD1 , const int *nD1,
                            const double *aD2, const double *rD2, const double *cD2, const int *iD2 , const int *nD2);
      bool   D_Coulomb_Y   (const double *aD1, const double *rD1, const double *cD1, const int *lD1,
                            const double *a2,  const double *r2,  const int *l2,     const int *m2, 
                                  double *value);
      bool   D_Coulomb_C   (const double *aD1, const double *rD1, const double *cD1, const int *lD1,
                            const double *a2,  const double *r2,  const int *nx2,    const int *ny2, const int *nz2,
                                  double *value);
      bool   D_Coulomb_Y_shell (const double *aD1, const double *rD1, const double *cD1, const int *lD1,
                                const double *a2,  const double *r2,  const int *l2, double *value);
      bool   D_Coulomb_C_shell (const double *aD1, const double *rD1, const double *cD1, const int *lD1,
                                const double *a2,  const double *r2,  const int *l2, double *value);
      bool   D_Coulomb_Y_shells(const double *aD1, const double *rD1, const double *cD1, const int *lD1,
                                const double *a2,  const double *x2,  const double *y2,  const double *z2,
                                const int *l2   ,  const int *nshell,       double *value);
      bool   D_Coulomb_C_shells(const double *aD1, const double *rD1, const double *cD1, const int *lD1,
                                const double *a2,  const double *x2,  const double *y2,  const double *z2,
                                const int *l2   ,  const int *nshell,       double *value);
      void   Y_Rotation_Matrix(const int *l, const double *R, double *Rl, int *Rl_dim);
      void   R_Rotation_Matrix(const int *l, const double *R, double *Rl, int *Rl_dim);
      double Potential_from_C    (const double *r_test, const double *alpha, const double *r, const int *nx, const int *ny, const int *nz);
      double Potential_from_Y    (const double *r_test, const double *alpha, const double *r, const int *l , const int *m );
      void   Field_from_C        (const double *r_test, const double *alpha, const double *r, const int *nx, const int *ny, const int *nz, double *E_test);
      void   Field_from_Y        (const double *r_test, const double *alpha, const double *r, const int *l , const int *m , double *E_test);
      void   Electrostatics_from_C(const double *r_test, const double *alpha, const double *r, const int *nx, const int *ny, const int *nz, double *V_test, double *E_test);
      void   Electrostatics_from_Y(const double *r_test, const double *alpha, const double *r, const int *l , const int *m , double *V_test, double *E_test);
      void   Cumulative_electrostatics_on_D(const double *r_source,  //< source charge and dipole location
                                            const double *q_source,  //< source charge value
                                            const double *dx_source, //< source x dipole value
                                            const double *dy_source, //< source y dipole value
                                            const double *dz_source, //< source z dipole value
                                            const int    *n_source,  //< number of source charges and dipoles
                                            const double *alpha,     //< exponent of hermit harmonic test charge distributions
                                            const double *r,         //< center of hermit harmonic test charge distributions
                                            const int    *l,         //< max l for hermit harmonic test charge distributions
                                                  double *V_test,    //< hermit harmonic test charge distributions cumulated interaction with source charge
                                                  double *Ex_test,   //< hermit harmonic test charge distributions cumulated interaction with source x dipoles
                                                  double *Ey_test,   //< hermit harmonic test charge distributions cumulated interaction with source y dipoles
                                                  double *Ez_test);  //< hermit harmonic test charge distributions cumulated interaction with source z dipoles
      void   Cumulative_electrostatics_on_C(const double *r_source,  //< source charge and dipole location
                                            const double *q_source,  //< source charge value
                                            const double *dx_source, //< source x dipole value
                                            const double *dy_source, //< source y dipole value
                                            const double *dz_source, //< source z dipole value
                                            const int    *n_source,  //< number of source charges and dipoles
                                            const double *alpha,     //< exponent of cubic harmonic test charge distribution
                                            const double *r,         //< center of cubic harmonic test charge distribution
                                            const int    *nx,        //< x power for cubic harmonic test charge distribution
                                            const int    *ny,        //< y power for cubic harmonic test charge distribution
                                            const int    *nz,        //< z power for cubic harmonic test charge distribution
                                                  double *V_test,    //< cubic harmonic test charge distribution cumulated interaction with source charge
                                                  double *Ex_test,   //< cubic harmonic test charge distribution cumulated interaction with source x dipoles
                                                  double *Ey_test,   //< cubic harmonic test charge distribution cumulated interaction with source y dipoles
                                                  double *Ez_test);  //< cubic harmonic test charge distribution cumulated interaction with source z dipoles
      void   Cumulative_electrostatics_on_Y(const double *r_source,  //< source charge and dipole location
                                            const double *q_source,  //< source charge value
                                            const double *dx_source, //< source x dipole value
                                            const double *dy_source, //< source y dipole value
                                            const double *dz_source, //< source z dipole value
                                            const int    *n_source,  //< number of source charges and dipoles
                                            const double *alpha,     //< exponent of solid harmonic test charge distribution
                                            const double *r,         //< center of solid harmonic test charge distribution
                                            const int    *l,         //< l momentum for solid harmonic test charge distribution
                                            const int    *m,         //< m momentum for solid harmonic test charge distribution
                                                  double *V_test,    //< solid harmonic test charge distribution cumulated interaction with source charge
                                                  double *Ex_test,   //< solid harmonic test charge distribution cumulated interaction with source x dipoles
                                                  double *Ey_test,   //< solid harmonic test charge distribution cumulated interaction with source y dipoles
                                                  double *Ez_test);  //< solid harmonic test charge distribution cumulated interaction with source z dipoles
      void   Cumulative_electrostatics_on_YY(const double *r_source,  //< source charge and dipole location
                                             const double *q_source,  //< source charge value
                                             const double *dx_source, //< source x dipole value
                                             const double *dy_source, //< source y dipole value
                                             const double *dz_source, //< source z dipole value
                                             const int    *n_source,  //< number of source charges and dipoles
                                             const double *alpha1,    //< exponent of first solid harmonic test charge distribution
                                             const double *r1,        //< center of first solid harmonic test charge distribution
                                             const int    *l1,        //< l momentum for first solid harmonic test charge distribution
                                             const int    *m1,        //< m momentum for first solid harmonic test charge distribution
                                             const double *alpha2,    //< exponent of second solid harmonic test charge distribution
                                             const double *r2,        //< center of second solid harmonic test charge distribution
                                             const int    *l2,        //< l momentum for second solid harmonic test charge distribution
                                             const int    *m2,        //< m momentum for second solid harmonic test charge distribution
                                                   double *V_test,    //< solid harmonic test charge distribution cumulated interaction with source charge
                                                   double *Ex_test,   //< solid harmonic test charge distribution cumulated interaction with source x dipoles
                                                   double *Ey_test,   //< solid harmonic test charge distribution cumulated interaction with source y dipoles
                                                   double *Ez_test);  //< solid harmonic test charge distribution cumulated interaction with source z dipoles
      
#ifdef __cplusplus
  }
#endif

#endif
