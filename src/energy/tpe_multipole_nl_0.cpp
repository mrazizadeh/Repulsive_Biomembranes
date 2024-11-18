
#include "energy/tpe_multipole_nl_0.h"

namespace biorsurfaces
{
    
    template<typename T1, typename T2>
    mreal TPEnergyMultipole_Normals0::FarField( T1 alpha, T2 betahalf)
    {
        T2 minus_betahalf = -betahalf;
        
        auto S = bct->S;
        auto T = bct->T;
        mint b_m = bct->far->b_m;
        mint  const * restrict const b_outer = bct->far->b_outer;
        mint  const * restrict const b_inner = bct->far->b_inner;
        mint nthreads = std::min( S->thread_count, T->thread_count);
        
        // Dunno why "restrict" helps with P_near. It is actually a lie here when S = T.
        // Well, it isn't in the far field, because no cluster may interact with itself...
        mreal const * restrict const A  = S->C_far[0];
        mreal const * restrict const X1 = S->C_far[1];
        mreal const * restrict const X2 = S->C_far[2];
        mreal const * restrict const X3 = S->C_far[3];
        mreal const * restrict const N1 = S->C_far[4];
        mreal const * restrict const N2 = S->C_far[5];
        mreal const * restrict const N3 = S->C_far[6];
        
        mreal const * restrict const B  = T->C_far[0];
        mreal const * restrict const Y1 = T->C_far[1];
        mreal const * restrict const Y2 = T->C_far[2];
        mreal const * restrict const Y3 = T->C_far[3];
        mreal const * restrict const M1 = T->C_far[4];
        mreal const * restrict const M2 = T->C_far[5];
        mreal const * restrict const M3 = T->C_far[6];
        
        mreal sum = 0.;
        
    #pragma omp parallel for num_threads( nthreads ) reduction( + : sum )
        for( mint i = 0; i < b_m; ++i )
        {
            mreal x1 = X1[i];
            mreal x2 = X2[i];
            mreal x3 = X3[i];
            mreal n1 = N1[i];
            mreal n2 = N2[i];
            mreal n3 = N3[i];
            
            mreal block_sum = 0.;
            
            // This loop can be SIMDized straight-forwardly (horizontal SIMDization).
            
            mint k_begin = b_outer[i];
            mint k_end = b_outer[i+1];
            #pragma omp simd aligned( B, Y1, Y2, Y3, M1, M2, M3 : ALIGN )
            for( mint k = k_begin; k < k_end; ++k )
            {
                mint j = b_inner[k];
                
                if(i <= j)
                {
                    mreal v1 = Y1[j] - x1;
                    mreal v2 = Y2[j] - x2;
                    mreal v3 = Y3[j] - x3;
                    mreal m1 = M1[j];
                    mreal m2 = M2[j];
                    mreal m3 = M3[j];
                    
                    mreal rCosPhi = v1 * n1 + v2 * n2 + v3 * n3;
                    mreal rCosPsi = v1 * m1 + v2 * m2 + v3 * m3;
                    mreal r2 = v1 * v1 + v2 * v2 + v3 * v3 ;
                    
                    block_sum += ( mypow( fabs(rCosPhi), alpha ) + mypow( fabs(rCosPsi), alpha) ) * mypow( r2, minus_betahalf ) * B[j];
                }
            }
            
            sum += A[i] * block_sum;
        }
        return sum;
    } // FarField


    template<typename T1, typename T2>
    mreal TPEnergyMultipole_Normals0::NearField(T1 alpha, T2 betahalf)
    {
        T2 minus_betahalf = -betahalf;
        
        auto S = bct->S;
        auto T = bct->T;
        mint b_m = bct->near->b_m;
        mint nthreads = std::min( S->thread_count, T->thread_count);
        
        mint  const * restrict const b_row_ptr = S->leaf_cluster_ptr;
        mint  const * restrict const b_col_ptr = T->leaf_cluster_ptr;
        
        mint  const * restrict const b_outer   = bct->near->b_outer;
        mint  const * restrict const b_inner   = bct->near->b_inner;
        
        // Dunno why "restrict" helps with P_near. It is actually a lie here.
        mreal const * restrict const A  = S->P_near[0];
        mreal const * restrict const X1 = S->P_near[1];
        mreal const * restrict const X2 = S->P_near[2];
        mreal const * restrict const X3 = S->P_near[3];
        mreal const * restrict const N1 = S->P_near[4];
        mreal const * restrict const N2 = S->P_near[5];
        mreal const * restrict const N3 = S->P_near[6];
        
        mreal const * restrict const B  = T->P_near[0];
        mreal const * restrict const Y1 = T->P_near[1];
        mreal const * restrict const Y2 = T->P_near[2];
        mreal const * restrict const Y3 = T->P_near[3];
        mreal const * restrict const M1 = T->P_near[4];
        mreal const * restrict const M2 = T->P_near[5];
        mreal const * restrict const M3 = T->P_near[6];
        
        mreal sum = 0.;
    #pragma omp parallel for num_threads( nthreads ) reduction( + : sum)
        for( mint b_i = 0; b_i < b_m; ++b_i )
        {
            
            mint i_begin = b_row_ptr[b_i];
            mint i_end   = b_row_ptr[b_i+1];
            
            for( mint k = b_outer[b_i]; k < b_outer[b_i+1]; ++k )
            {
                mint b_j = b_inner[k];
                if(b_i <= b_j)
                {
                    mint j_begin = b_col_ptr[b_j];
                    mint j_end   = b_col_ptr[b_j+1];
                    mreal block_sum = 0.;
                    
                    for( mint i = i_begin; i < i_end ; ++i )
                    {
                        mreal x1 = X1[i];
                        mreal x2 = X2[i];
                        mreal x3 = X3[i];
                        mreal n1 = N1[i];
                        mreal n2 = N2[i];
                        mreal n3 = N3[i];
                        
                        mreal i_sum = 0.;
                        
                        // if b_i == b_j, we loop only over the upper triangular block, diagonal excluded
                        mint begin = (b_i != b_j) ? j_begin : i + 1;
                        // Here, one could do a bit of horizontal vectorization. However, the number of js an x interacts with varies greatly..
                        #pragma omp simd aligned( B, Y1, Y2, Y3, M1, M2, M3 : ALIGN ) reduction( + : block_sum )
                        for( mint j = begin; j < j_end; ++j )
                        {
                            if( i != j )
                            {
                                mreal v1 = Y1[j] - x1;
                                mreal v2 = Y2[j] - x2;
                                mreal v3 = Y3[j] - x3;
                                mreal m1 = M1[j];
                                mreal m2 = M2[j];
                                mreal m3 = M3[j];
                                
                                mreal rCosPhi = v1 * n1 + v2 * n2 + v3 * n3;
                                mreal rCosPsi = v1 * m1 + v2 * m2 + v3 * m3;
                                mreal r2 = v1 * v1 + v2 * v2 + v3 * v3 ;
                                
                                mreal en = ( mypow( fabs(rCosPhi), alpha ) + mypow( fabs(rCosPsi), alpha) ) * mypow( r2, minus_betahalf );
                                
                                
                                i_sum += en * B[j];
                            }
                        }
                        block_sum += A[i] * i_sum;
                    }
                    sum += block_sum;
                }
            }
        }
        return sum;
    }; //NearField

    template<typename T1, typename T2>
    mreal TPEnergyMultipole_Normals0::DFarField(T1 alpha, T2 betahalf)
    {
        
        T1 alpha_minus_2 = alpha - 2;
        T2 minus_betahalf_minus_1 = -betahalf - 1;
        
        mreal beta = 2. * betahalf;
        
        mreal sum = 0.;
        
        auto S = bct->S;
        auto T = bct->T;

        mint b_m = bct->far->b_m;
        mint nthreads = std::min( S->thread_count, T->thread_count);
        
        mint  const * restrict const b_outer = bct->far->b_outer;
        mint  const * restrict const b_inner = bct->far->b_inner;
        
        
        // Dunno why "restrict" helps with C_far. It is actually a lie here.
        mreal const * restrict const A  = S->C_far[0];
        mreal const * restrict const X1 = S->C_far[1];
        mreal const * restrict const X2 = S->C_far[2];
        mreal const * restrict const X3 = S->C_far[3];
        mreal const * restrict const N1 = S->C_far[4];
        mreal const * restrict const N2 = S->C_far[5];
        mreal const * restrict const N3 = S->C_far[6];
        
        mreal const * restrict const B  = T->C_far[0];
        mreal const * restrict const Y1 = T->C_far[1];
        mreal const * restrict const Y2 = T->C_far[2];
        mreal const * restrict const Y3 = T->C_far[3];
        mreal const * restrict const M1 = T->C_far[4];
        mreal const * restrict const M2 = T->C_far[5];
        mreal const * restrict const M3 = T->C_far[6];
        
    #pragma omp parallel for num_threads( nthreads ) reduction( + : sum )
        for( mint i = 0; i < b_m; ++i )
        {
            mint thread = omp_get_thread_num();
            
            mreal * restrict const U = &S->C_D_far[thread][0];
            mreal * restrict const V = &T->C_D_far[thread][0];
            
            mreal a  =  A[i];
            mreal x1 = X1[i];
            mreal x2 = X2[i];
            mreal x3 = X3[i];
            mreal n1 = N1[i];
            mreal n2 = N2[i];
            mreal n3 = N3[i];
            
            mreal  da = 0.;
            mreal dx1 = 0.;
            mreal dx2 = 0.;
            mreal dx3 = 0.;
            mreal dn1 = 0.;
            mreal dn2 = 0.;
            mreal dn3 = 0.;
            
            mreal block_sum = 0.;
            
            mint k_begin = b_outer[i];
            mint k_end = b_outer[i+1];
            
            // This loop can be SIMDized straight-forwardly (horizontal SIMDization).
            #pragma omp simd aligned( B, Y1, Y2, Y3, M1, M2, M3 : ALIGN ) reduction( + : block_sum )
            for( mint k = k_begin; k < k_end; ++k )
            {
                mint j = b_inner[k];
                
                if( i <= j )
                {
                    mreal  b = B [j];
                    mreal y1 = Y1[j];
                    mreal y2 = Y2[j];
                    mreal y3 = Y3[j];
                    mreal m1 = M1[j];
                    mreal m2 = M2[j];
                    mreal m3 = M3[j];
                    
                    mreal v1 = y1 - x1;
                    mreal v2 = y2 - x2;
                    mreal v3 = y3 - x3;
                    
                    mreal rCosPhi = v1 * n1 + v2 * n2 + v3 * n3;
                    mreal rCosPsi = v1 * m1 + v2 * m2 + v3 * m3;
                    mreal r2      = v1 * v1 + v2 * v2 + v3 * v3;
                    
                    mreal rBetaMinus2 = mypow( r2, minus_betahalf_minus_1 );
                    mreal rBeta = rBetaMinus2 * r2;
                    
                    mreal rCosPhiAlphaMinus1 = mypow( fabs(rCosPhi), alpha_minus_2 ) * rCosPhi;
                    mreal rCosPhiAlpha = rCosPhiAlphaMinus1 * rCosPhi;
                    
                    mreal rCosPsiAlphaMinus1 = mypow( fabs(rCosPsi), alpha_minus_2 ) * rCosPsi;
                    mreal rCosPsiAlpha = rCosPsiAlphaMinus1 * rCosPsi;
                    
                    
                    mreal Num = rCosPhiAlpha + rCosPsiAlpha;
                    mreal factor0 = rBeta * alpha;
                    mreal density = rBeta * Num;
                    block_sum += a * b * density;
                    
                    mreal F = factor0 * rCosPhiAlphaMinus1;
                    mreal G = factor0 * rCosPsiAlphaMinus1;
                    mreal H = beta * rBetaMinus2 * Num;
                    
                    mreal bF = b * F;
                    mreal aG = a * G;
                    
                    mreal Z1 = ( - n1 * F - m1 * G + v1 * H );
                    mreal Z2 = ( - n2 * F - m2 * G + v2 * H );
                    mreal Z3 = ( - n3 * F - m3 * G + v3 * H );
                    
                    da += b * (
                               density
                               +
                               F * ( n1 * (x1 - v1) + n2 * (x2 - v2) + n3 * (x3 - v3) )
                               +
                               G * ( m1 * x1 + m2 * x2 + m3 * x3 )
                               -
                               H * ( v1 * x1 + v2 * x2 + v3 * x3 )
                               );
                    
                    V[ 7 * j ] += a * (
                                              density
                                              -
                                              F * ( n1 * y1 + n2 * y2 + n3 * y3 )
                                              -
                                              G * ( m1 * (y1 + v1) + m2 * (y2 + v2) + m3 * (y3 + v3) )
                                              +
                                              H * ( v1 * y1 + v2 * y2 + v3 * y3 )
                                              );
                    
                    dx1 += b  * Z1;
                    dx2 += b  * Z2;
                    dx3 += b  * Z3;
                    dn1 += bF * v1;
                    dn2 += bF * v2;
                    dn3 += bF * v3;
                    
                    V[ 7 * j + 1 ] -= a  * Z1;
                    V[ 7 * j + 2 ] -= a  * Z2;
                    V[ 7 * j + 3 ] -= a  * Z3;
                    V[ 7 * j + 4 ] += aG * v1;
                    V[ 7 * j + 5 ] += aG * v2;
                    V[ 7 * j + 6 ] += aG * v3;
                    
                } // if( i<= j )
            } // for( mint k = b_outer[i]; k < b_outer[i+1]; ++k )
            
            U[ 7 * i     ] +=  da;
            U[ 7 * i + 1 ] += dx1;
            U[ 7 * i + 2 ] += dx2;
            U[ 7 * i + 3 ] += dx3;
            U[ 7 * i + 4 ] += dn1;
            U[ 7 * i + 5 ] += dn2;
            U[ 7 * i + 6 ] += dn3;
            
            sum += block_sum;
        } // for( mint i = 0; i < b_m; ++i )
        
        return sum;
    }; //DFarField

    template<typename T1, typename T2>
    mreal TPEnergyMultipole_Normals0::DNearField(T1 alpha, T2 betahalf)
    {
        T1 alpha_minus_2 = alpha - 2;
        T2 minus_betahalf_minus_1 = -betahalf - 1;
        
        mreal beta = 2. * betahalf;
        
        mreal sum = 0.;
        
        auto S = bct->S;
        auto T = bct->T;
        
        mint b_m = bct->near->b_m;
        mint nthreads = std::min( S->thread_count, T->thread_count);
        
        mint  const * restrict const b_row_ptr = S->leaf_cluster_ptr;
        mint  const * restrict const b_col_ptr = T->leaf_cluster_ptr;
        
        mint  const * restrict const b_outer = &bct->near->b_outer[0];
        mint  const * restrict const b_inner = &bct->near->b_inner[0];
        
        // Dunno why "restrict" helps with P_near. It is actually a lie here.
        mreal const * restrict const A  = S->P_near[0];
        mreal const * restrict const X1 = S->P_near[1];
        mreal const * restrict const X2 = S->P_near[2];
        mreal const * restrict const X3 = S->P_near[3];
        mreal const * restrict const N1 = S->P_near[4];
        mreal const * restrict const N2 = S->P_near[5];
        mreal const * restrict const N3 = S->P_near[6];
        
        mreal const * restrict const B  = T->P_near[0];
        mreal const * restrict const Y1 = T->P_near[1];
        mreal const * restrict const Y2 = T->P_near[2];
        mreal const * restrict const Y3 = T->P_near[3];
        mreal const * restrict const M1 = T->P_near[4];
        mreal const * restrict const M2 = T->P_near[5];
        mreal const * restrict const M3 = T->P_near[6];
        
        
        #pragma omp parallel for num_threads( nthreads ) reduction( +: sum )
        for( mint b_i = 0; b_i < b_m; ++b_i )
        {
            mint thread = omp_get_thread_num();
            
            mreal * restrict const U = &S->P_D_near[thread][0];
            mreal * restrict const V = &T->P_D_near[thread][0];
            
            mint i_begin = b_row_ptr[b_i];
            mint i_end   = b_row_ptr[b_i+1];
            
            for( mint k = b_outer[b_i]; k < b_outer[b_i+1]; ++k )
            {
                mint b_j = b_inner[k];
                if( b_i <= b_j )
                {
                    mint j_begin = b_col_ptr[b_j];
                    mint j_end   = b_col_ptr[b_j+1];
                    
                    for( mint i = i_begin; i < i_end ; ++i )
                    {
                        mreal  a = A [i];
                        mreal x1 = X1[i];
                        mreal x2 = X2[i];
                        mreal x3 = X3[i];
                        mreal n1 = N1[i];
                        mreal n2 = N2[i];
                        mreal n3 = N3[i];
                        
                        mreal  da = 0.;
                        mreal dx1 = 0.;
                        mreal dx2 = 0.;
                        mreal dx3 = 0.;
                        mreal dn1 = 0.;
                        mreal dn2 = 0.;
                        mreal dn3 = 0.;
                        
                        
                        // if i == j, we loop only over the upper triangular block, diagonal excluded
                        mint begin = (b_i != b_j) ? j_begin : i + 1;
                        
                        // Here, one could do a bit of horizontal vectorization. However, the number of js an x interacts with is small and varies greatly..
                        #pragma omp simd aligned( B, Y1, Y2, Y3, M1, M2, M3 : ALIGN ) reduction( + : sum)
                        for( mint j = begin; j < j_end; ++j )
                        {
                            if( i != j )
                            {
                                mreal  b = B [j];
                                mreal y1 = Y1[j];
                                mreal y2 = Y2[j];
                                mreal y3 = Y3[j];
                                mreal m1 = M1[j];
                                mreal m2 = M2[j];
                                mreal m3 = M3[j];
                                
                                mreal v1 = y1 - x1;
                                mreal v2 = y2 - x2;
                                mreal v3 = y3 - x3;
                                
                                mreal rCosPhi = v1 * n1 + v2 * n2 + v3 * n3;
                                mreal rCosPsi = v1 * m1 + v2 * m2 + v3 * m3;
                                mreal r2      = v1 * v1 + v2 * v2 + v3 * v3;
                                
                                mreal rBetaMinus2 = mypow( r2, minus_betahalf_minus_1 );
                                mreal rBeta = rBetaMinus2 * r2;
                                
                                mreal rCosPhiAlphaMinus1 = mypow( fabs(rCosPhi), alpha_minus_2 ) * rCosPhi;
                                mreal rCosPhiAlpha = rCosPhiAlphaMinus1 * rCosPhi;
                                
                                mreal rCosPsiAlphaMinus1 = mypow( fabs(rCosPsi), alpha_minus_2 ) * rCosPsi;
                                mreal rCosPsiAlpha = rCosPsiAlphaMinus1 * rCosPsi;
                                
                                
                                mreal Num = rCosPhiAlpha + rCosPsiAlpha;
                                mreal factor0 = rBeta * alpha;
                                mreal density = rBeta * Num;
                                sum += a * b * density;
                                
                                mreal F = factor0 * rCosPhiAlphaMinus1;
                                mreal G = factor0 * rCosPsiAlphaMinus1;
                                mreal H = beta * rBetaMinus2 * Num;
                                
                                mreal bF = b * F;
                                mreal aG = a * G;
                                
                                mreal Z1 = ( - n1 * F - m1 * G + v1 * H );
                                mreal Z2 = ( - n2 * F - m2 * G + v2 * H );
                                mreal Z3 = ( - n3 * F - m3 * G + v3 * H );
                                
                                da += b * (
                                           density
                                           +
                                           F * ( n1 * (x1 - v1) + n2 * (x2 - v2) + n3 * (x3 - v3) )
                                           +
                                           G * ( m1 * x1 + m2 * x2 + m3 * x3 )
                                           -
                                           H * ( v1 * x1 + v2 * x2 + v3 * x3 )
                                           );
                                
                                V[ 7 * j ] += a * (
                                                          density
                                                          -
                                                          F * ( n1 * y1 + n2 * y2 + n3 * y3 )
                                                          -
                                                          G * ( m1 * (y1 + v1) + m2 * (y2 + v2) + m3 * (y3 + v3) )
                                                          +
                                                          H * ( v1 * y1 + v2 * y2 + v3 * y3 )
                                                          );
                                
                                dx1 += b  * Z1;
                                dx2 += b  * Z2;
                                dx3 += b  * Z3;
                                dn1 += bF * v1;
                                dn2 += bF * v2;
                                dn3 += bF * v3;
                                
                                V[ 7 * j + 1 ] -= a  * Z1;
                                V[ 7 * j + 2 ] -= a  * Z2;
                                V[ 7 * j + 3 ] -= a  * Z3;
                                V[ 7 * j + 4 ] += aG * v1;
                                V[ 7 * j + 5 ] += aG * v2;
                                V[ 7 * j + 6 ] += aG * v3;
                            }
                        } // for( mint j = ( b_i != b_j ? j_begin : i + 1 ); j < j_end; ++j )
                        
                        
                        U[ 7 * i     ] +=  da;
                        U[ 7 * i + 1 ] += dx1;
                        U[ 7 * i + 2 ] += dx2;
                        U[ 7 * i + 3 ] += dx3;
                        U[ 7 * i + 4 ] += dn1;
                        U[ 7 * i + 5 ] += dn2;
                        U[ 7 * i + 6 ] += dn3;
                        
                    }// for( mint i = i_begin; i < i_end ; ++i )
                }// if( b_i <= b_j )
            }// for( mint k = b_outer[b_i]; k < b_outer[b_i+1]; ++k )
        } // for( mint b_i = 0; b_i < b_m; ++b_i )
        
        return sum;
    }; //DNearField


    // Returns the current value of the energy.
    double TPEnergyMultipole_Normals0::Value()
    {
        
        if( use_int )
        {
            mint int_alpha = std::round(alpha);
            mint int_betahalf = std::round(beta/2);
            return weight * (FarField( int_alpha, int_betahalf ) + NearField (int_alpha, int_betahalf ));
        }
        else
        {
            mreal real_alpha = alpha;
            mreal real_betahalf = beta/2;
            return weight * (FarField( real_alpha, real_betahalf ) + NearField( real_alpha, real_betahalf ));
        }
    } // Value

    // Returns the current differential of the energy, stored in the given
    // V x 3 matrix, where each row holds the differential (a 3-vector) with
    // respect to the corresponding vertex.
    void TPEnergyMultipole_Normals0::Differential(Eigen::MatrixXd &output)
    {
        if( bct->S->near_dim != 7)
        {
            eprint("in TPEnergyBarnesHut_Projectors0::Differential: near_dim != 7");
        }
        if( bct->S->far_dim != 7)
        {
            eprint("in TPEnergyBarnesHut_Projectors0::Differential: far_dim != 7");
        }
                
        EigenMatrixRM P_D_near( bct->S->primitive_count, bct->S->near_dim );
        EigenMatrixRM P_D_far ( bct->S->primitive_count, bct->S->far_dim );
        
        bct->S->CleanseD();
//        bct->T->CleanseD();
        
        if( use_int )
        {
            mint int_alpha = std::round(alpha);
            mint int_betahalf = std::round(beta/2);
            DNearField( int_alpha, int_betahalf );
            DFarField ( int_alpha, int_betahalf );
            
        }
        else
        {
            mreal real_alpha = alpha;
            mreal real_betahalf = beta/2;
            DNearField( real_alpha, real_betahalf );
            DFarField ( real_alpha, real_betahalf );
        }
        
        bct->S->CollectDerivatives( P_D_near.data(), P_D_far.data() );
                
        AssembleDerivativeFromACNData( mesh, geom, P_D_near, output, weight );
        AssembleDerivativeFromACNData( mesh, geom, P_D_far, output, weight );
        
    } // Differential


    // Update the energy to reflect the current state of the mesh. This could
    // involve building a new BVH for Barnes-Hut energies, for instance.
    void TPEnergyMultipole_Normals0::Update()
    {
        throw std::runtime_error("Multipole energy not supported for flow");
    }

    // Get the mesh associated with this energy.
    MeshPtr TPEnergyMultipole_Normals0::GetMesh()
    {
        return mesh;
    }

    // Get the geometry associated with this geometry.
    GeomPtr TPEnergyMultipole_Normals0::GetGeom()
    {
        return geom;
    }

    // Get the exponents of this energy; only applies to tangent-point energies.
    Vector2 TPEnergyMultipole_Normals0::GetExponents()
    {
        return Vector2{alpha, beta};
    }

    // Get a pointer to the current BVH for this energy.
    // Return 0 if the energy doesn't use a BVH.
    OptimizedClusterTree * TPEnergyMultipole_Normals0::GetBVH()
    {
        return 0;
    }

    // Return the separation parameter for this energy.
    // Return 0 if this energy doesn't do hierarchical approximation.
    double TPEnergyMultipole_Normals0::GetTheta()
    {
        return sqrt(bct->theta2);
    }
    
} // namespace biorsurfaces
