/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <math.h>
#include <stdint.h>

double my_erf_ultimate_ref(double a)
{
  // fit parameters
  const double p1=0.131506884694059871829665553377708429e0;
  const double p2=0.504665528409229025008642111060908103e0;
  const double C[16]={
   4.29463428976381402415959464760932200e0,
   0.01865143915153509079509433526017075e0,
 -14.20868456394042349279792943707128640e0,
  -0.19756010201375526523019673104134380e0,
  23.62028718318021596636366125371244120e0,
   0.89263080906490788552186786030257620e0,
 -22.24394322043539279890633770070421540e0,
  -2.05914903824048358200942904768149371e0,
  13.53812364637492398439908583302899160e0,
   1.79205237867671881372831820553640887e0,
  -4.95872327050368255234677868091306920e0,
   2.32435844794304634100765629026017495e0,
   0.64814731776613766124188486514352200e0,
  -5.68817140081665065814356854440259275e0,
  -1.16492141026395367274102441919508674e0,
   4.39226749429304471743993461394588706e0};
  
  // calculation variabes
  double t1;
  double t2;
  double s1;
  double s2;
  
  // set t's
  t1 = 1.0e0 / ( 1.0e0 + p1 * a );
  t2 = 1.0e0 / ( 1.0e0 + p2 * a );
  
  // compute sum
  s1 = C[0];
  s2 = C[1];
  int j;
  for ( j=1 ; j<8 ; j++ )
  {
    s1 = C[2*j+0] + t1 * s1;
    s2 = C[2*j+1] + t2 * s2;
  }
  
  // last comp
  return 1.0e0 - exp(-a*a) * ( t1 * s1 + t2 * s2 );
}


double my_erf_ultimate(double a)
{
  // calculation variabes
  double t1;
  double t2;
  double s1;
  double s2;
  
  // set t's
  t1 = 1.0e0 / ( 1.0e0 + 0.131506884694059871829665553377708429e0 * a );
  t2 = 1.0e0 / ( 1.0e0 + 0.504665528409229025008642111060908103e0 * a );
  
  // compute sum
  s1 = 4.29463428976381402415959464760932200e0;
  s2 = 0.01865143915153509079509433526017075e0;
  s1 = t1 * s1 - 14.20868456394042349279792943707128640e0;
  s2 = t2 * s2 -  0.19756010201375526523019673104134380e0;
  s1 = t1 * s1 + 23.62028718318021596636366125371244120e0;
  s2 = t2 * s2 +  0.89263080906490788552186786030257620e0;
  s1 = t1 * s1 - 22.24394322043539279890633770070421540e0;
  s2 = t2 * s2 -  2.05914903824048358200942904768149371e0;
  s1 = t1 * s1 + 13.53812364637492398439908583302899160e0;
  s2 = t2 * s2 +  1.79205237867671881372831820553640887e0;
  s1 = t1 * s1 -  4.95872327050368255234677868091306920e0;
  s2 = t2 * s2 +  2.32435844794304634100765629026017495e0;
  s1 = t1 * s1 +  0.64814731776613766124188486514352200e0;
  s2 = t2 * s2 -  5.68817140081665065814356854440259275e0;
  s1 = t1 * s1 -  1.16492141026395367274102441919508674e0;
  s2 = t2 * s2 +  4.39226749429304471743993461394588706e0;
  
  // last comp
  return copysign( 1.0e0 - exp(-a*a) * ( t1 * s1 + t2 * s2 ), a);
}

#ifdef __AVX2__

#include <immintrin.h>

void my_erf_v4_(const double * a, double * c, double * b){
  
//  for ( int i=0 ; i<4  ; i++ ) c[i]=exp(-a[i]*a[i]);
//  for ( int i=0 ; i<4  ; i++ ) b[i]=erf(a[i]);
//  return;
  
  int32_t k;
  __m256d s10;
  __m256d s20;
  __m256d t10;
  __m256d t20;
  __m256d a20;
  t10 = _mm256_loadu_pd( &a[ + (4) * (0)] );
  t20 = _mm256_loadu_pd( &a[ + (4) * (0)] );
  a20 = _mm256_mul_pd( t10, t10 );
  t10 = _mm256_div_pd( _mm256_set1_pd( 1.0 ), _mm256_fmadd_pd( _mm256_set1_pd( 0.13150688469405988 ), t10, _mm256_set1_pd( 1.0 ) ) );
  t20 = _mm256_div_pd( _mm256_set1_pd( 1.0 ), _mm256_fmadd_pd( _mm256_set1_pd( 0.5046655284092291 ), t20, _mm256_set1_pd( 1.0 ) ) );
  s10 = _mm256_set1_pd( 4.294634289763814 );
  s20 = _mm256_set1_pd( 0.018651439151535092 );
  s10 = _mm256_fmadd_pd( t10, s10, _mm256_set1_pd( -14.208684563940423 ) );
  s20 = _mm256_fmadd_pd( t20, s20, _mm256_set1_pd( -0.19756010201375526 ) );
  s10 = _mm256_fmadd_pd( t10, s10, _mm256_set1_pd( 23.620287183180217 ) );
  s20 = _mm256_fmadd_pd( t20, s20, _mm256_set1_pd( 0.8926308090649079 ) );
  s10 = _mm256_fmadd_pd( t10, s10, _mm256_set1_pd( -22.243943220435394 ) );
  s20 = _mm256_fmadd_pd( t20, s20, _mm256_set1_pd( -2.059149038240484 ) );
  s10 = _mm256_fmadd_pd( t10, s10, _mm256_set1_pd( 13.538123646374924 ) );
  s20 = _mm256_fmadd_pd( t20, s20, _mm256_set1_pd( 1.7920523786767188 ) );
  s10 = _mm256_fmadd_pd( t10, s10, _mm256_set1_pd( -4.958723270503683 ) );
  s20 = _mm256_fmadd_pd( t20, s20, _mm256_set1_pd( 2.3243584479430464 ) );
  s10 = _mm256_fmadd_pd( t10, s10, _mm256_set1_pd( 0.6481473177661377 ) );
  s20 = _mm256_fmadd_pd( t20, s20, _mm256_set1_pd( -5.688171400816651 ) );
  s10 = _mm256_fmadd_pd( t10, s10, _mm256_set1_pd( -1.1649214102639536 ) );
  s20 = _mm256_fmadd_pd( t20, s20, _mm256_set1_pd( 4.392267494293045 ) );
  for (k = 0; k <= 3; k += 1) {
    a20[k] = exp( -(a20[k]));
  }
  _mm256_storeu_pd( (double * ) &c[ + (4) * (0)], a20 );
  _mm256_storeu_pd( (double * ) &b[ + (4) * (0)], _mm256_fmadd_pd(  -(a20), _mm256_fmadd_pd( t10, s10, _mm256_mul_pd( t20, s20 ) ), _mm256_set1_pd( 1.0 ) ) );
}

#else

void my_erf_v4_(const double * a, double * c, double * b)
{
  int i;
  for ( i=0 ; i<4  ; i++ ) c[i]=exp(-a[i]*a[i]);
  for ( i=0 ; i<4  ; i++ ) b[i]=erf(a[i]);
  return;
}


#endif

