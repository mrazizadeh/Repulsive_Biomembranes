
#include "energy/tp_obstacle_barnes_hut_pr_0.h"

namespace biorsurfaces
{
    
    template <typename T1, typename T2>
    mreal TPObstacleBarnesHut_Projectors0::Energy(T1 alphahalf, T2 betahalf)
    {
        T2 minus_betahalf = -betahalf;
        mreal theta2 = theta * theta;

        mint nthreads = bvh->thread_count;

        mreal sum = 0.;

        {
            auto S = bvh;
            auto T = o_bvh;
            mreal const * restrict const  C_xmin1 = S->C_min[0];
            mreal const * restrict const  C_xmin2 = S->C_min[1];
            mreal const * restrict const  C_xmin3 = S->C_min[2];
            mreal const * restrict const  C_xmax1 = S->C_max[0];
            mreal const * restrict const  C_xmax2 = S->C_max[1];
            mreal const * restrict const  C_xmax3 = S->C_max[2];

            mreal const * restrict const  C_xr2 = S->C_squared_radius;
            
            mreal const * restrict const P_A  = S->P_near[0];
            mreal const * restrict const P_X1 = S->P_near[1];
            mreal const * restrict const P_X2 = S->P_near[2];
            mreal const * restrict const P_X3 = S->P_near[3];
            mreal const * restrict const P_P11 = S->P_near[4];
            mreal const * restrict const P_P12 = S->P_near[5];
            mreal const * restrict const P_P13 = S->P_near[6];
            mreal const * restrict const P_P22 = S->P_near[7];
            mreal const * restrict const P_P23 = S->P_near[8];
            mreal const * restrict const P_P33 = S->P_near[9];
            

            mint  const * restrict const  C_xbegin = S->C_begin;
            mint  const * restrict const  C_xend = S->C_end;

            mint  const * restrict const  leaf = S->leaf_clusters;

            mreal const * restrict const  C_ymin1 = T->C_min[0];
            mreal const * restrict const  C_ymin2 = T->C_min[1];
            mreal const * restrict const  C_ymin3 = T->C_min[2];
            mreal const * restrict const  C_ymax1 = T->C_max[0];
            mreal const * restrict const  C_ymax2 = T->C_max[1];
            mreal const * restrict const  C_ymax3 = T->C_max[2];

            mreal const * restrict const  C_yr2 = T->C_squared_radius;

            mreal const * restrict const  P_B = T->P_near[0];
            mreal const * restrict const  P_Y1 = T->P_near[1];
            mreal const * restrict const  P_Y2 = T->P_near[2];
            mreal const * restrict const  P_Y3 = T->P_near[3];

            mreal const * restrict const  C_B = T->C_far[0];
            mreal const * restrict const  C_Y1 = T->C_far[1];
            mreal const * restrict const  C_Y2 = T->C_far[2];
            mreal const * restrict const  C_Y3 = T->C_far[3];

            mint  const * restrict const  C_ybegin = T->C_begin;
            mint  const * restrict const  C_yend = T->C_end;

            mint  const * restrict const  C_left = T->C_left;
            mint  const * restrict const  C_right = T->C_right;

            A_Vector<A_Vector<mint>> thread_stack(nthreads);

            #pragma omp parallel for num_threads(nthreads) reduction(+ : sum)
            for (mint k = 0; k < S->leaf_cluster_count; ++k)
            {
                mint thread = omp_get_thread_num();

                A_Vector<mint> *stack = &thread_stack[thread];

                stack->clear();
                stack->push_back(0);

                mint l = leaf[k];
                mint i_begin = C_xbegin[l];
                mint i_end = C_xend[l];

                mreal xmin1 = C_xmin1[l];
                mreal xmin2 = C_xmin2[l];
                mreal xmin3 = C_xmin3[l];

                mreal xmax1 = C_xmax1[l];
                mreal xmax2 = C_xmax2[l];
                mreal xmax3 = C_xmax3[l];

                mreal r2l = C_xr2[l];

                mreal local_sum = 0.;

                while (!stack->empty())
                {
                    mint C = stack->back();
                    stack->pop_back();

                    mreal h2 = std::max(r2l, C_yr2[C]);

                    // Compute squared distance between bounding boxes.
                    // Inpired by https://gamedev.stackexchange.com/questions/154036/efficient-minimum-distance-between-two-axis-aligned-squares

                    mreal ymin1 = C_ymin1[C];
                    mreal ymin2 = C_ymin2[C];
                    mreal ymin3 = C_ymin3[C];

                    mreal ymax1 = C_ymax1[C];
                    mreal ymax2 = C_ymax2[C];
                    mreal ymax3 = C_ymax3[C];

                    mreal d1 = mymax(0., mymax(xmin1, ymin1) - mymin(xmax1, ymax1));
                    mreal d2 = mymax(0., mymax(xmin2, ymin2) - mymin(xmax2, ymax2));
                    mreal d3 = mymax(0., mymax(xmin3, ymin3) - mymin(xmax3, ymax3));

                    mreal R2 = d1 * d1 + d2 * d2 + d3 * d3;

                    if (h2 < theta2 * R2)
                    {
                        mreal b = C_B[C];
                        mreal y1 = C_Y1[C];
                        mreal y2 = C_Y2[C];
                        mreal y3 = C_Y3[C];

                        mreal local_local_sum = 0.;

                        for (mint i = i_begin; i < i_end; ++i)
                        {
                            mreal a = P_A[i];
                            mreal x1 = P_X1[i];
                            mreal x2 = P_X2[i];
                            mreal x3 = P_X3[i];
                            mreal p11 = P_P11[i];
                            mreal p12 = P_P12[i];
                            mreal p13 = P_P13[i];
                            mreal p22 = P_P22[i];
                            mreal p23 = P_P23[i];
                            mreal p33 = P_P33[i];

                            mreal v1 = y1 - x1;
                            mreal v2 = y2 - x2;
                            mreal v3 = y3 - x3;

                            mreal rCosPhi2 = v1*(p11*v1 + p12*v2 + p13*v3) + v2*(p12*v1 + p22*v2 + p23*v3) + v3*(p13*v1 + p23*v2 + p33*v3);
                            mreal r2 = v1 * v1 + v2 * v2 + v3 * v3;
                            local_local_sum += a * mypow( fabs(rCosPhi2), alphahalf ) * mypow(r2, minus_betahalf);
                        }
                        local_sum += local_local_sum * b;
                    }
                    else
                    {
                        mint left = C_left[C];
                        mint right = C_right[C];
                        if (left >= 0 && right >= 0)
                        {
                            stack->push_back(right);
                            stack->push_back(left);
                        }
                        else
                        {
                            // near field loop
                            mint j_begin = C_ybegin[C];
                            mint j_end = C_yend[C];

                            for (mint i = i_begin; i < i_end; ++i)
                            {
                                mreal a = P_A[i];
                                mreal x1 = P_X1[i];
                                mreal x2 = P_X2[i];
                                mreal x3 = P_X3[i];
                                mreal p11 = P_P11[i];
                                mreal p12 = P_P12[i];
                                mreal p13 = P_P13[i];
                                mreal p22 = P_P22[i];
                                mreal p23 = P_P23[i];
                                mreal p33 = P_P33[i];

                                mreal local_local_sum = 0.;

                                //                        #pragma omp simd aligned( P_A, P_X1, P_X3 : ALIGN )
                                for (mint j = j_begin; j < j_end; ++j)
                                {
                                    mreal b = P_B[j];
                                    mreal v1 = P_Y1[j] - x1;
                                    mreal v2 = P_Y2[j] - x2;
                                    mreal v3 = P_Y3[j] - x3;

                                    mreal rCosPhi2 = v1*(p11*v1 + p12*v2 + p13*v3) + v2*(p12*v1 + p22*v2 + p23*v3) + v3*(p13*v1 + p23*v2 + p33*v3);
                                    mreal r2 = v1 * v1 + v2 * v2 + v3 * v3;

                                    local_local_sum += mypow( fabs(rCosPhi2), alphahalf )  * mypow(r2, minus_betahalf) * b;
                                }

                                local_sum += a * local_local_sum;
                            }
                        }
                    }
                }

                sum += local_sum;
            }
        }

        {
            auto S = o_bvh;
            auto T = bvh;
            mreal const * restrict const  C_xmin1 = S->C_min[0];
            mreal const * restrict const  C_xmin2 = S->C_min[1];
            mreal const * restrict const  C_xmin3 = S->C_min[2];
            mreal const * restrict const  C_xmax1 = S->C_max[0];
            mreal const * restrict const  C_xmax2 = S->C_max[1];
            mreal const * restrict const  C_xmax3 = S->C_max[2];

            mreal const * restrict const  C_xr2 = S->C_squared_radius;

            mreal const * restrict const P_A  = S->P_near[0];
            mreal const * restrict const P_X1 = S->P_near[1];
            mreal const * restrict const P_X2 = S->P_near[2];
            mreal const * restrict const P_X3 = S->P_near[3];
            mreal const * restrict const P_P11 = S->P_near[4];
            mreal const * restrict const P_P12 = S->P_near[5];
            mreal const * restrict const P_P13 = S->P_near[6];
            mreal const * restrict const P_P22 = S->P_near[7];
            mreal const * restrict const P_P23 = S->P_near[8];
            mreal const * restrict const P_P33 = S->P_near[9];

            mint  const * restrict const  C_xbegin = S->C_begin;
            mint  const * restrict const  C_xend = S->C_end;

            mint  const * restrict const  leaf = S->leaf_clusters;

            mreal const * restrict const  C_ymin1 = T->C_min[0];
            mreal const * restrict const  C_ymin2 = T->C_min[1];
            mreal const * restrict const  C_ymin3 = T->C_min[2];
            mreal const * restrict const  C_ymax1 = T->C_max[0];
            mreal const * restrict const  C_ymax2 = T->C_max[1];
            mreal const * restrict const  C_ymax3 = T->C_max[2];

            mreal const * restrict const  C_yr2 = T->C_squared_radius;

            mreal const * restrict const  P_B = T->P_near[0];
            mreal const * restrict const  P_Y1 = T->P_near[1];
            mreal const * restrict const  P_Y2 = T->P_near[2];
            mreal const * restrict const  P_Y3 = T->P_near[3];

            mreal const * restrict const  C_B = T->C_far[0];
            mreal const * restrict const  C_Y1 = T->C_far[1];
            mreal const * restrict const  C_Y2 = T->C_far[2];
            mreal const * restrict const  C_Y3 = T->C_far[3];

            mint  const * restrict const  C_ybegin = T->C_begin;
            mint  const * restrict const  C_yend = T->C_end;

            mint  const * restrict const  C_left = T->C_left;
            mint  const * restrict const  C_right = T->C_right;

            A_Vector<A_Vector<mint>> thread_stack(nthreads);

            #pragma omp parallel for num_threads(nthreads) reduction(+ : sum)
            for (mint k = 0; k < S->leaf_cluster_count; ++k)
            {
                mint thread = omp_get_thread_num();

                A_Vector<mint> *stack = &thread_stack[thread];

                stack->clear();
                stack->push_back(0);

                mint l = leaf[k];
                mint i_begin = C_xbegin[l];
                mint i_end = C_xend[l];

                mreal xmin1 = C_xmin1[l];
                mreal xmin2 = C_xmin2[l];
                mreal xmin3 = C_xmin3[l];

                mreal xmax1 = C_xmax1[l];
                mreal xmax2 = C_xmax2[l];
                mreal xmax3 = C_xmax3[l];

                mreal r2l = C_xr2[l];

                mreal local_sum = 0.;

                while (!stack->empty())
                {
                    mint C = stack->back();
                    stack->pop_back();

                    mreal h2 = std::max(r2l, C_yr2[C]);

                    // Compute squared distance between bounding boxes.
                    // Inpired by https://gamedev.stackexchange.com/questions/154036/efficient-minimum-distance-between-two-axis-aligned-squares

                    mreal ymin1 = C_ymin1[C];
                    mreal ymin2 = C_ymin2[C];
                    mreal ymin3 = C_ymin3[C];

                    mreal ymax1 = C_ymax1[C];
                    mreal ymax2 = C_ymax2[C];
                    mreal ymax3 = C_ymax3[C];

                    mreal d1 = mymax(0., mymax(xmin1, ymin1) - mymin(xmax1, ymax1));
                    mreal d2 = mymax(0., mymax(xmin2, ymin2) - mymin(xmax2, ymax2));
                    mreal d3 = mymax(0., mymax(xmin3, ymin3) - mymin(xmax3, ymax3));

                    mreal R2 = d1 * d1 + d2 * d2 + d3 * d3;

                    if (h2 < theta2 * R2)
                    {
                        mreal b = C_B[C];
                        mreal y1 = C_Y1[C];
                        mreal y2 = C_Y2[C];
                        mreal y3 = C_Y3[C];

                        mreal local_local_sum = 0.;

                        for (mint i = i_begin; i < i_end; ++i)
                        {
                            mreal a = P_A[i];
                            mreal x1 = P_X1[i];
                            mreal x2 = P_X2[i];
                            mreal x3 = P_X3[i];
                            mreal p11 = P_P11[i];
                            mreal p12 = P_P12[i];
                            mreal p13 = P_P13[i];
                            mreal p22 = P_P22[i];
                            mreal p23 = P_P23[i];
                            mreal p33 = P_P33[i];

                            mreal v1 = y1 - x1;
                            mreal v2 = y2 - x2;
                            mreal v3 = y3 - x3;

                            mreal rCosPhi2 = v1*(p11*v1 + p12*v2 + p13*v3) + v2*(p12*v1 + p22*v2 + p23*v3) + v3*(p13*v1 + p23*v2 + p33*v3);
                            mreal r2 = v1 * v1 + v2 * v2 + v3 * v3;
                            local_local_sum += a * mypow( fabs(rCosPhi2), alphahalf ) * mypow(r2, minus_betahalf);
                        }
                        local_sum += local_local_sum * b;
                    }
                    else
                    {
                        mint left = C_left[C];
                        mint right = C_right[C];
                        if (left >= 0 && right >= 0)
                        {
                            stack->push_back(right);
                            stack->push_back(left);
                        }
                        else
                        {
                            // near field loop
                            mint j_begin = C_ybegin[C];
                            mint j_end = C_yend[C];

                            for (mint i = i_begin; i < i_end; ++i)
                            {
                                mreal a = P_A[i];
                                mreal x1 = P_X1[i];
                                mreal x2 = P_X2[i];
                                mreal x3 = P_X3[i];
                                mreal p11 = P_P11[i];
                                mreal p12 = P_P12[i];
                                mreal p13 = P_P13[i];
                                mreal p22 = P_P22[i];
                                mreal p23 = P_P23[i];
                                mreal p33 = P_P33[i];

                                mreal local_local_sum = 0.;

                                //                        #pragma omp simd aligned( P_A, P_X1, P_X3 : ALIGN )
                                for (mint j = j_begin; j < j_end; ++j)
                                {
                                    mreal b = P_B[j];
                                    mreal v1 = P_Y1[j] - x1;
                                    mreal v2 = P_Y2[j] - x2;
                                    mreal v3 = P_Y3[j] - x3;

                                    mreal rCosPhi2 = v1*(p11*v1 + p12*v2 + p13*v3) + v2*(p12*v1 + p22*v2 + p23*v3) + v3*(p13*v1 + p23*v2 + p33*v3);
                                    mreal r2 = v1 * v1 + v2 * v2 + v3 * v3;

                                    local_local_sum += mypow( fabs(rCosPhi2), alphahalf ) * mypow(r2, minus_betahalf) * b;
                                }

                                local_sum += a * local_local_sum;
                            }
                        }
                    }
                }

                sum += local_sum;
            }
        }

        return sum;
    }; //Energy
    
    template <typename T1, typename T2>
    mreal TPObstacleBarnesHut_Projectors0::DEnergy(T1 alphahalf, T2 betahalf)
    {

        T1 alphahalf_minus_1 = alphahalf - 1;
        T2 minus_betahalf_minus_1 = -betahalf - 1;

        mreal beta = 2. * betahalf;
        mreal theta2 = theta * theta;
        mreal sum = 0.;

        mint nthreads = bvh->thread_count;

        {
            auto S = bvh;
            auto T = o_bvh;
            mreal const * restrict const  C_xmin1 = S->C_min[0];
            mreal const * restrict const  C_xmin2 = S->C_min[1];
            mreal const * restrict const  C_xmin3 = S->C_min[2];
            mreal const * restrict const  C_xmax1 = S->C_max[0];
            mreal const * restrict const  C_xmax2 = S->C_max[1];
            mreal const * restrict const  C_xmax3 = S->C_max[2];

            mreal const * restrict const  C_xr2 = S->C_squared_radius;

            mreal const * restrict const  P_A = S->P_near[0];
            mreal const * restrict const  P_X1 = S->P_near[1];
            mreal const * restrict const  P_X2 = S->P_near[2];
            mreal const * restrict const  P_X3 = S->P_near[3];
            mreal const * restrict const  P_P11 = S->P_near[4];
            mreal const * restrict const  P_P12 = S->P_near[5];
            mreal const * restrict const  P_P13 = S->P_near[6];
            mreal const * restrict const  P_P22 = S->P_near[7];
            mreal const * restrict const  P_P23 = S->P_near[8];
            mreal const * restrict const  P_P33 = S->P_near[9];

            mint  const * restrict const  C_xbegin = S->C_begin;
            mint  const * restrict const  C_xend = S->C_end;

            mint  const * restrict const  leaf = S->leaf_clusters;

            mreal const * restrict const  C_ymin1 = T->C_min[0];
            mreal const * restrict const  C_ymin2 = T->C_min[1];
            mreal const * restrict const  C_ymin3 = T->C_min[2];
            mreal const * restrict const  C_ymax1 = T->C_max[0];
            mreal const * restrict const  C_ymax2 = T->C_max[1];
            mreal const * restrict const  C_ymax3 = T->C_max[2];

            mreal const * restrict const  C_yr2 = T->C_squared_radius;

            mreal const * restrict const  P_B = T->P_near[0];
            mreal const * restrict const  P_Y1 = T->P_near[1];
            mreal const * restrict const  P_Y2 = T->P_near[2];
            mreal const * restrict const  P_Y3 = T->P_near[3];

            mreal const * restrict const  C_B = T->C_far[0];
            mreal const * restrict const  C_Y1 = T->C_far[1];
            mreal const * restrict const  C_Y2 = T->C_far[2];
            mreal const * restrict const  C_Y3 = T->C_far[3];

            mint  const * restrict const  C_ybegin = T->C_begin;
            mint  const * restrict const  C_yend = T->C_end;

            mint  const * restrict const  C_left = T->C_left;
            mint  const * restrict const  C_right = T->C_right;

            A_Vector<A_Vector<mint>> thread_stack(nthreads);

            #pragma omp parallel for num_threads(nthreads) reduction(+ : sum)
            for (mint k = 0; k < S->leaf_cluster_count; ++k)
            {
                mint thread = omp_get_thread_num();

                A_Vector<mint> *stack = &thread_stack[thread];

                mreal *restrict const P_U = &S->P_D_near[thread][0];

                stack->clear();
                stack->push_back(0);

                mint l = leaf[k];
                mint i_begin = C_xbegin[l];
                mint i_end = C_xend[l];

                mreal xmin1 = C_xmin1[l];
                mreal xmin2 = C_xmin2[l];
                mreal xmin3 = C_xmin3[l];

                mreal xmax1 = C_xmax1[l];
                mreal xmax2 = C_xmax2[l];
                mreal xmax3 = C_xmax3[l];

                mreal r2l = C_xr2[l];

                while (!stack->empty())
                {
                    mint C = stack->back();
                    stack->pop_back();

                    mreal h2 = std::max(r2l, C_yr2[C]);

                    // Compute squared distance between bounding boxes.
                    // Inpired by https://gamedev.stackexchange.com/questions/154036/efficient-minimum-distance-between-two-axis-aligned-squares

                    mreal ymin1 = C_ymin1[C];
                    mreal ymin2 = C_ymin2[C];
                    mreal ymin3 = C_ymin3[C];

                    mreal ymax1 = C_ymax1[C];
                    mreal ymax2 = C_ymax2[C];
                    mreal ymax3 = C_ymax3[C];

                    mreal d1 = mymax(0., mymax(xmin1, ymin1) - mymin(xmax1, ymax1));
                    mreal d2 = mymax(0., mymax(xmin2, ymin2) - mymin(xmax2, ymax2));
                    mreal d3 = mymax(0., mymax(xmin3, ymin3) - mymin(xmax3, ymax3));

                    mreal R2 = d1 * d1 + d2 * d2 + d3 * d3;

                    if (h2 < theta2 * R2)
                    {
                        mreal b = C_B[C];
                        mreal y1 = C_Y1[C];
                        mreal y2 = C_Y2[C];
                        mreal y3 = C_Y3[C];

                        for (mint i = i_begin; i < i_end; ++i)
                        {
                            mreal a = P_A[i];
                            mreal x1 = P_X1[i];
                            mreal x2 = P_X2[i];
                            mreal x3 = P_X3[i];
                            mreal p11 = P_P11[i];
                            mreal p12 = P_P12[i];
                            mreal p13 = P_P13[i];
                            mreal p22 = P_P22[i];
                            mreal p23 = P_P23[i];
                            mreal p33 = P_P33[i];

                            mreal v1 = y1 - x1;
                            mreal v2 = y2 - x2;
                            mreal v3 = y3 - x3;

                            mreal v11 = v1 * v1;
                            mreal v22 = v2 * v2;
                            mreal v33 = v3 * v3;
                            
                            mreal v12 = 2. * v1 * v2;
                            mreal v13 = 2. * v1 * v3;
                            mreal v23 = 2. * v2 * v3;
                            mreal r2 = v11 + v22 + v33;
                            
                            mreal Pv1 = p11*v1 + p12*v2 + p13*v3;
                            mreal Pv2 = p12*v1 + p22*v2 + p23*v3;
                            mreal Pv3 = p13*v1 + p23*v2 + p33*v3;
                            mreal rCosPhi2 = v1*Pv1 + v2*Pv2 + v3*Pv3;
                            
                            mreal rCosPhiAlphaMinus2 = mypow( fabs(rCosPhi2), alphahalf_minus_1);
                            mreal rMinusBetaMinus2 = mypow( r2, minus_betahalf_minus_1 );
                            
                            mreal rMinusBeta = rMinusBetaMinus2 * r2;
                            mreal rCosPhiAlpha = rCosPhiAlphaMinus2 * rCosPhi2;
                            mreal Num = rCosPhiAlpha;
                        
                            mreal E = Num * rMinusBeta;
                            sum += a * b * E;
                            
                            mreal factor = alphahalf * rMinusBeta;
                            mreal F = factor * rCosPhiAlphaMinus2;
                            mreal H = - beta * rMinusBetaMinus2 * Num;
                            
                            mreal bF = b * F;
                            
                            mreal dEdv1 = 2. * F * Pv1 + H * v1;
                            mreal dEdv2 = 2. * F * Pv2 + H * v2;
                            mreal dEdv3 = 2. * F * Pv3 + H * v3;
                            
                            P_U[10 * i + 0] += b * ( E + dEdv1 * x1 + dEdv2 * x2 + dEdv3 * x3 - factor * rCosPhiAlpha );
                            P_U[10 * i + 1] -= b * dEdv1;
                            P_U[10 * i + 2] -= b * dEdv2;
                            P_U[10 * i + 3] -= b * dEdv3;
                            P_U[10 * i + 4] += bF * v11;
                            P_U[10 * i + 5] += bF * v12;
                            P_U[10 * i + 6] += bF * v13;
                            P_U[10 * i + 7] += bF * v22;
                            P_U[10 * i + 8] += bF * v23;
                            P_U[10 * i + 9] += bF * v33;
                        }
                    }
                    else
                    {
                        mint left = C_left[C];
                        mint right = C_right[C];
                        if (left >= 0 && right >= 0)
                        {
                            stack->push_back(right);
                            stack->push_back(left);
                        }
                        else
                        {
                            // near field loop
                            mint j_begin = C_ybegin[C];
                            mint j_end = C_yend[C];

                            for (mint i = i_begin; i < i_end; ++i)
                            {
                                mreal a = P_A[i];
                                mreal x1 = P_X1[i];
                                mreal x2 = P_X2[i];
                                mreal x3 = P_X3[i];
                                mreal p11 = P_P11[i];
                                mreal p12 = P_P12[i];
                                mreal p13 = P_P13[i];
                                mreal p22 = P_P22[i];
                                mreal p23 = P_P23[i];
                                mreal p33 = P_P33[i];

                                mreal da = 0.;
                                mreal dx1 = 0.;
                                mreal dx2 = 0.;
                                mreal dx3 = 0.;
                                mreal dp11 = 0.;
                                mreal dp12 = 0.;
                                mreal dp13 = 0.;
                                mreal dp22 = 0.;
                                mreal dp23 = 0.;
                                mreal dp33 = 0.;

//                                #pragma omp simd aligned(P_B, P_Y1, P_Y2, P_Y3 : ALIGN) reduction(+ : sum)
                                for (mint j = j_begin; j < j_end; ++j)
                                {
                                    mreal b = P_B[j];
                                    mreal y1 = P_Y1[j];
                                    mreal y2 = P_Y2[j];
                                    mreal y3 = P_Y3[j];

                                    mreal v1 = y1 - x1;
                                    mreal v2 = y2 - x2;
                                    mreal v3 = y3 - x3;

                                    mreal v11 = v1 * v1;
                                    mreal v22 = v2 * v2;
                                    mreal v33 = v3 * v3;
                                    
                                    mreal v12 = 2. * v1 * v2;
                                    mreal v13 = 2. * v1 * v3;
                                    mreal v23 = 2. * v2 * v3;
                                    mreal r2 = v11 + v22 + v33;
                                    
                                    mreal Pv1 = p11*v1 + p12*v2 + p13*v3;
                                    mreal Pv2 = p12*v1 + p22*v2 + p23*v3;
                                    mreal Pv3 = p13*v1 + p23*v2 + p33*v3;
                                    mreal rCosPhi2 = v1*Pv1 + v2*Pv2 + v3*Pv3;
                                    
                                    mreal rCosPhiAlphaMinus2 = mypow( fabs(rCosPhi2), alphahalf_minus_1);
                                    mreal rMinusBetaMinus2 = mypow( r2, minus_betahalf_minus_1 );
                                    
                                    mreal rMinusBeta = rMinusBetaMinus2 * r2;
                                    mreal rCosPhiAlpha = rCosPhiAlphaMinus2 * rCosPhi2;
                                    mreal Num = rCosPhiAlpha;
                                
                                    mreal E = Num * rMinusBeta;
                                    sum += a * b * E;
                                    
                                    mreal factor = alphahalf * rMinusBeta;
                                    mreal F = factor * rCosPhiAlphaMinus2;
                                    mreal H = - beta * rMinusBetaMinus2 * Num;
                                    
                                    mreal bF = b * F;
                                    
                                    mreal dEdv1 = 2. * F * Pv1 + H * v1;
                                    mreal dEdv2 = 2. * F * Pv2 + H * v2;
                                    mreal dEdv3 = 2. * F * Pv3 + H * v3;
                                    
                                    da += b * ( E + dEdv1 * x1 + dEdv2 * x2 + dEdv3 * x3 - factor * rCosPhiAlpha );
                                    dx1 -= b * dEdv1;
                                    dx2 -= b * dEdv2;
                                    dx3 -= b * dEdv3;
                                    dp11 += bF * v11;
                                    dp12 += bF * v12;
                                    dp13 += bF * v13;
                                    dp22 += bF * v22;
                                    dp23 += bF * v23;
                                    dp33 += bF * v33;
                                }
                                P_U[ 10 * i + 0 ] +=  da;
                                P_U[ 10 * i + 1 ] += dx1;
                                P_U[ 10 * i + 2 ] += dx2;
                                P_U[ 10 * i + 3 ] += dx3;
                                P_U[ 10 * i + 4 ] += dp11;
                                P_U[ 10 * i + 5 ] += dp12;
                                P_U[ 10 * i + 6 ] += dp13;
                                P_U[ 10 * i + 7 ] += dp22;
                                P_U[ 10 * i + 8 ] += dp23;
                                P_U[ 10 * i + 9 ] += dp33;

                            }
                        }
                    }
                }
            }
        }

        {
            auto S = o_bvh;
            auto T = bvh;
            mreal const * restrict const  C_xmin1 = S->C_min[0];
            mreal const * restrict const  C_xmin2 = S->C_min[1];
            mreal const * restrict const  C_xmin3 = S->C_min[2];
            mreal const * restrict const  C_xmax1 = S->C_max[0];
            mreal const * restrict const  C_xmax2 = S->C_max[1];
            mreal const * restrict const  C_xmax3 = S->C_max[2];
            
            mreal const * restrict const  C_xr2 = S->C_squared_radius;
            
            mreal const * restrict const  P_A = S->P_near[0];
            mreal const * restrict const  P_X1 = S->P_near[1];
            mreal const * restrict const  P_X2 = S->P_near[2];
            mreal const * restrict const  P_X3 = S->P_near[3];
            mreal const * restrict const  P_P11 = S->P_near[4];
            mreal const * restrict const  P_P12 = S->P_near[5];
            mreal const * restrict const  P_P13 = S->P_near[6];
            mreal const * restrict const  P_P22 = S->P_near[7];
            mreal const * restrict const  P_P23 = S->P_near[8];
            mreal const * restrict const  P_P33 = S->P_near[9];
            
            mint  const * restrict const  C_xbegin = S->C_begin;
            mint  const * restrict const  C_xend = S->C_end;
            
            mint  const * restrict const  leaf = S->leaf_clusters;
            
            mreal const * restrict const  C_ymin1 = T->C_min[0];
            mreal const * restrict const  C_ymin2 = T->C_min[1];
            mreal const * restrict const  C_ymin3 = T->C_min[2];
            mreal const * restrict const  C_ymax1 = T->C_max[0];
            mreal const * restrict const  C_ymax2 = T->C_max[1];
            mreal const * restrict const  C_ymax3 = T->C_max[2];
            
            mreal const * restrict const  C_yr2 = T->C_squared_radius;
            
            mreal const * restrict const  P_B = T->P_near[0];
            mreal const * restrict const  P_Y1 = T->P_near[1];
            mreal const * restrict const  P_Y2 = T->P_near[2];
            mreal const * restrict const  P_Y3 = T->P_near[3];
            
            mreal const * restrict const  C_B = T->C_far[0];
            mreal const * restrict const  C_Y1 = T->C_far[1];
            mreal const * restrict const  C_Y2 = T->C_far[2];
            mreal const * restrict const  C_Y3 = T->C_far[3];
            
            mint  const * restrict const  C_ybegin = T->C_begin;
            mint  const * restrict const  C_yend = T->C_end;
            
            mint  const * restrict const  C_left = T->C_left;
            mint  const * restrict const  C_right = T->C_right;
            
            A_Vector<A_Vector<mint>> thread_stack(nthreads);
            
            #pragma omp parallel for num_threads(nthreads) reduction(+ : sum)
            for (mint k = 0; k < S->leaf_cluster_count; ++k)
            {
                mint thread = omp_get_thread_num();
                
                A_Vector<mint> *stack = &thread_stack[thread];
                
                mreal * restrict const P_V = &T->P_D_near[thread][0];
                mreal * restrict const C_V = &T->C_D_far[thread][0];
                
                stack->clear();
                stack->push_back(0);
                
                mint l = leaf[k];
                mint i_begin = C_xbegin[l];
                mint i_end = C_xend[l];
                
                mreal xmin1 = C_xmin1[l];
                mreal xmin2 = C_xmin2[l];
                mreal xmin3 = C_xmin3[l];
                
                mreal xmax1 = C_xmax1[l];
                mreal xmax2 = C_xmax2[l];
                mreal xmax3 = C_xmax3[l];
                
                mreal r2l = C_xr2[l];
                
                while (!stack->empty())
                {
                    mint C = stack->back();
                    stack->pop_back();
                    
                    mreal h2 = std::max(r2l, C_yr2[C]);
                    
                    // Compute squared distance between bounding boxes.
                    // Inpired by https://gamedev.stackexchange.com/questions/154036/efficient-minimum-distance-between-two-axis-aligned-squares
                    
                    mreal ymin1 = C_ymin1[C];
                    mreal ymin2 = C_ymin2[C];
                    mreal ymin3 = C_ymin3[C];
                    
                    mreal ymax1 = C_ymax1[C];
                    mreal ymax2 = C_ymax2[C];
                    mreal ymax3 = C_ymax3[C];
                    
                    mreal d1 = mymax(0., mymax(xmin1, ymin1) - mymin(xmax1, ymax1));
                    mreal d2 = mymax(0., mymax(xmin2, ymin2) - mymin(xmax2, ymax2));
                    mreal d3 = mymax(0., mymax(xmin3, ymin3) - mymin(xmax3, ymax3));
                    
                    mreal R2 = d1 * d1 + d2 * d2 + d3 * d3;
                    
                    if (h2 < theta2 * R2)
                    {
                        mreal b = C_B[C];
                        mreal y1 = C_Y1[C];
                        mreal y2 = C_Y2[C];
                        mreal y3 = C_Y3[C];
                        
                        mreal db = 0.;
                        mreal dy1 = 0.;
                        mreal dy2 = 0.;
                        mreal dy3 = 0.;
                        
//                        #pragma omp simd aligned(P_A, P_X1, P_X2, P_X3, P_P11, P_P12, P_P13, P_P22, P_P23, P_P33 : ALIGN) reduction(+ : sum)
                        for (mint i = i_begin; i < i_end; ++i)
                        {
                            mreal a = P_A[i];
                            mreal x1 = P_X1[i];
                            mreal x2 = P_X2[i];
                            mreal x3 = P_X3[i];
                            mreal p11 = P_P11[i];
                            mreal p12 = P_P12[i];
                            mreal p13 = P_P13[i];
                            mreal p22 = P_P22[i];
                            mreal p23 = P_P23[i];
                            mreal p33 = P_P33[i];
                            
                            mreal v1 = y1 - x1;
                            mreal v2 = y2 - x2;
                            mreal v3 = y3 - x3;
                            
                            mreal v11 = v1 * v1;
                            mreal v22 = v2 * v2;
                            mreal v33 = v3 * v3;
                            
                            mreal v12 = 2. * v1 * v2;
                            mreal v13 = 2. * v1 * v3;
                            mreal v23 = 2. * v2 * v3;
                            mreal r2 = v11 + v22 + v33;
                            
                            mreal Pv1 = p11*v1 + p12*v2 + p13*v3;
                            mreal Pv2 = p12*v1 + p22*v2 + p23*v3;
                            mreal Pv3 = p13*v1 + p23*v2 + p33*v3;
                            mreal rCosPhi2 = v1*Pv1 + v2*Pv2 + v3*Pv3;
                            
                            mreal rCosPhiAlphaMinus2 = mypow( fabs(rCosPhi2), alphahalf_minus_1);
                            mreal rMinusBetaMinus2 = mypow( r2, minus_betahalf_minus_1 );
                            
                            mreal rMinusBeta = rMinusBetaMinus2 * r2;
                            mreal rCosPhiAlpha = rCosPhiAlphaMinus2 * rCosPhi2;
                            mreal Num = rCosPhiAlpha;
                        
                            mreal E = Num * rMinusBeta;
                            sum += a * b * E;
                            
                            mreal factor = alphahalf * rMinusBeta;
                            mreal F = factor * rCosPhiAlphaMinus2;
                            mreal H = - beta * rMinusBetaMinus2 * Num;
                            
                            mreal bF = b * F;
                            
                            mreal dEdv1 = 2. * F * Pv1 + H * v1;
                            mreal dEdv2 = 2. * F * Pv2 + H * v2;
                            mreal dEdv3 = 2. * F * Pv3 + H * v3;
                            
                            db += a * ( E - dEdv1 * y1 - dEdv2 * y2 - dEdv3 * y3);
                            dy1 += a * dEdv1;
                            dy2 += a * dEdv2;
                            dy3 += a * dEdv3;
                        }
                        C_V[10 * C + 0] += db;
                        C_V[10 * C + 1] += dy1;
                        C_V[10 * C + 2] += dy2;
                        C_V[10 * C + 3] += dy3;
                    }
                    else
                    {
                        mint left = C_left[C];
                        mint right = C_right[C];
                        if (left >= 0 && right >= 0)
                        {
                            stack->push_back(right);
                            stack->push_back(left);
                        }
                        else
                        {
                            // near field loop
                            mint j_begin = C_ybegin[C];
                            mint j_end = C_yend[C];
                            
                            for (mint i = i_begin; i < i_end; ++i)
                            {
                                mreal a = P_A[i];
                                mreal x1 = P_X1[i];
                                mreal x2 = P_X2[i];
                                mreal x3 = P_X3[i];
                                mreal p11 = P_P11[i];
                                mreal p12 = P_P12[i];
                                mreal p13 = P_P13[i];
                                mreal p22 = P_P22[i];
                                mreal p23 = P_P23[i];
                                mreal p33 = P_P33[i];
                                
//                                #pragma omp simd aligned(P_B, P_Y1, P_Y2, P_Y3, P_V : ALIGN) reduction(+ : sum)
                                for (mint j = j_begin; j < j_end; ++j)
                                {
                                    mreal b = P_B[j];
                                    mreal y1 = P_Y1[j];
                                    mreal y2 = P_Y2[j];
                                    mreal y3 = P_Y3[j];
                                    
                                    mreal v1 = y1 - x1;
                                    mreal v2 = y2 - x2;
                                    mreal v3 = y3 - x3;
                                    
                                    mreal v11 = v1 * v1;
                                    mreal v22 = v2 * v2;
                                    mreal v33 = v3 * v3;
                                    
                                    mreal v12 = 2. * v1 * v2;
                                    mreal v13 = 2. * v1 * v3;
                                    mreal v23 = 2. * v2 * v3;
                                    mreal r2 = v11 + v22 + v33;
                                    
                                    mreal Pv1 = p11*v1 + p12*v2 + p13*v3;
                                    mreal Pv2 = p12*v1 + p22*v2 + p23*v3;
                                    mreal Pv3 = p13*v1 + p23*v2 + p33*v3;
                                    mreal rCosPhi2 = v1*Pv1 + v2*Pv2 + v3*Pv3;
                                    
                                    mreal rCosPhiAlphaMinus2 = mypow( fabs(rCosPhi2), alphahalf_minus_1);
                                    mreal rMinusBetaMinus2 = mypow( r2, minus_betahalf_minus_1 );
                                    
                                    mreal rMinusBeta = rMinusBetaMinus2 * r2;
                                    mreal rCosPhiAlpha = rCosPhiAlphaMinus2 * rCosPhi2;
                                    mreal Num = rCosPhiAlpha;
                                
                                    mreal E = Num * rMinusBeta;
                                    sum += a * b * E;
                                    
                                    mreal factor = alphahalf * rMinusBeta;
                                    mreal F = factor * rCosPhiAlphaMinus2;
                                    mreal H = - beta * rMinusBetaMinus2 * Num;
                                    
                                    mreal bF = b * F;
                                    
                                    mreal dEdv1 = 2. * F * Pv1 + H * v1;
                                    mreal dEdv2 = 2. * F * Pv2 + H * v2;
                                    mreal dEdv3 = 2. * F * Pv3 + H * v3;
                                    
                                    P_V[10 * j + 0] += a * ( E - dEdv1 * y1 - dEdv2 * y2 - dEdv3 * y3);
                                    P_V[10 * j + 1] += a * dEdv1;
                                    P_V[10 * j + 2] += a * dEdv2;
                                    P_V[10 * j + 3] += a * dEdv3;
                                    
                                }
                            }
                        }
                    }
                }
            }
        }
        
        return sum;
    }; // DEnergy
    
    double TPObstacleBarnesHut_Projectors0::Value()
    {
        bvh = bvhSharedFrom->GetBVH();
        if (!bvh)
        {
            throw std::runtime_error("Obstacle energy is sharing BVH from an energy that has no BVH.");
        }
        if (use_int)
        {
            mint int_alphahalf = std::round(alpha/2);
            mint int_betahalf = std::round(beta/2);
            return weight * Energy(int_alphahalf, int_betahalf);
        }
        else
        {
            mreal real_alphahalf = alpha/2;
            mreal real_betahalf = beta/2;
            return weight * Energy(real_alphahalf, real_betahalf);
        }
    } // Value

    void TPObstacleBarnesHut_Projectors0::Differential(Eigen::MatrixXd &output)
    {
        bvh = bvhSharedFrom->GetBVH();
        if (!bvh)
        {
            throw std::runtime_error("Obstacle energy is sharing BVH from an energy that has no BVH.");
        }
        
        if( bvh->near_dim != 10)
        {
            eprint("in TPObstacleBarnesHut_Projectors0::Differential: bvh->near_dim != 10");
        }
        if( bvh->far_dim != 10)
        {
            eprint("in TPObstacleBarnesHut_Projectors0::Differential: bvh->far_dim != 10");
        }
        
        if( o_bvh->near_dim != 10)
        {
            eprint("in TPObstacleBarnesHut_Projectors0::Differential: o_bvh->near_dim != 10");
        }
        if( o_bvh->far_dim != 10)
        {
            eprint("in TPObstacleBarnesHut_Projectors0::Differential: o_bvh->far_dim != 10");
        }
        
        EigenMatrixRM P_D_near_( bvh->primitive_count, 10 );
        EigenMatrixRM P_D_far_ ( bvh->primitive_count, 10 );
        
        bvh->CleanseD();
        
        if( use_int )
        {
            mint int_alphahalf = std::round(alpha/2);
            mint int_betahalf = std::round(beta/2);
            DEnergy( int_alphahalf, int_betahalf );
            
        }
        else
        {
            mreal real_alphahalf = alpha/2;
            mreal real_betahalf = beta/2;
            DEnergy( real_alphahalf, real_betahalf );
        }
        
        bvh->CollectDerivatives( P_D_near_.data(), P_D_far_.data() );
    
        AssembleDerivativeFromACPData( mesh, geom, P_D_near_, output, weight );
        AssembleDerivativeFromACPData( mesh, geom, P_D_far_ , output, weight );
  
    } // Differential
    
    // Update the energy to reflect the current state of the mesh. This could
    // involve building a new BVH for Barnes-Hut energies, for instance.
    void TPObstacleBarnesHut_Projectors0::Update()
    {
        // Invalidate the old BVH pointer
        bvh = 0;
        // bvhSharedFrom is responsible for reallocating it in its Update() function
        bvh = bvhSharedFrom->GetBVH();
        if (!bvh)
        {
            throw std::runtime_error("Obstacle energy is sharing BVH from an energy that has no BVH.");
        }
    }

    // Get the mesh associated with this energy.
    MeshPtr TPObstacleBarnesHut_Projectors0::GetMesh()
    {
        return mesh;
    }

    // Get the geometry associated with this geometry.
    GeomPtr TPObstacleBarnesHut_Projectors0::GetGeom()
    {
        return geom;
    }

    // Get the exponents of this energy; only applies to tangent-point energies.
    Vector2 TPObstacleBarnesHut_Projectors0::GetExponents()
    {
        return Vector2{alpha, beta};
    }

    // Get a pointer to the current BVH for this energy.
    // Return 0 if the energy doesn't use a BVH.
    OptimizedClusterTree *TPObstacleBarnesHut_Projectors0::GetBVH()
    {
        return o_bvh;
    }

    // Return the separation parameter for this energy.
    // Return 0 if this energy doesn't do hierarchical approximation.
    double TPObstacleBarnesHut_Projectors0::GetTheta()
    {
        return theta;
    }
} // namespace biorsurfaces
