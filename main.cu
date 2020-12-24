#include<iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <cuda.h>
#include <time.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <iomanip>

#include "include/params.cuh"
#include "include/DataStructure.cuh"
#include "include/MCMPC.cuh"
#include "include/init.cuh"
#include "include/cuSolverForMCMPC.cuh"

#define Linear

void printMatrix(int m, int n, float*A, int lda, const char* name)
{
    for(int row = 0 ; row < m ; row++){
        for(int col = 0 ; col < n ; col++){
            float Areg = A[row + col*lda];
            printf("%s(%d,%d) = %f\n", name, row+1, col+1, Areg);
            //printf("%s[%d] = %f\n", name, row + col*lda, Areg);
        }
    }
}
int main(int argc, char **argv)
{
    cusolverDnHandle_t cusolverH = NULL;
    cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;
    cudaError_t cudaStat1 = cudaSuccess;
    cudaError_t cudaStat2 = cudaSuccess;
    cudaError_t cudaStat3 = cudaSuccess;

    cusolver_status = cusolverDnCreate(&cusolverH);
    assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
    /*データ書き込みファイルの定義*/
    FILE *fp;
    time_t timeValue;
    struct tm *timeObject;
    time( &timeValue );
    timeObject = localtime( &timeValue );
    char filename1[35];
    sprintf(filename1,"data_system_%d%d_%d%d.txt",timeObject->tm_mon + 1, timeObject->tm_mday, timeObject->tm_hour,timeObject->tm_min);
    fp = fopen(filename1,"w");


    float params[dim_param], state[dim_state], /*h_constraint[NUM_CONST],*/ h_matrix[dim_weight_matrix];
    float *device_param, *device_matrix;
    Mat_sys_A( params );
    init_state( state );
    // init_constraint( h_constraint );
    init_Weight_matrix( h_matrix );
    //cudaMemcpyToSymbol(d_param, &params, dim_param * sizeof(float));
    cudaMalloc(&device_param, sizeof(float)*dim_param);
    cudaMalloc(&device_matrix, sizeof(float)*dim_weight_matrix);
    cudaMemcpy(device_param, params, sizeof(float)*dim_param, cudaMemcpyHostToDevice);
    cudaMemcpy(device_matrix, h_matrix, sizeof(float)*dim_weight_matrix, cudaMemcpyHostToDevice);
    

#ifdef Linear
    float opt[HORIZON], Error[HORIZON];
    init_opt( opt );
#endif


    /* GPUの設定 */
    unsigned int numBlocks, randomBlocks, randomNums/*, minId_cpu*/;
    int Blocks;
    randomNums = N_OF_SAMPLES * (dim_U+1) * HORIZON;
    randomBlocks = countBlocks(randomNums, THREAD_PER_BLOCKS);
    numBlocks = countBlocks(N_OF_SAMPLES, THREAD_PER_BLOCKS);
    printf("#NumBlocks = %d\n", numBlocks);
    Blocks = numBlocks;

    /* CPU to GPU dataExchanger */
    Data1 *h_dataFromBlocks;
    Data1 *d_dataFromBlocks;

#ifdef USING_THRUST
    thrust::device_vector<int> indices_device_vec( N_OF_SAMPLES );
    //thrust::device_vector<int> indices_vec_dev_temp( N_OF_SAMPLES );
    //indices_device_vec = indices_vec_dev_temp;
    thrust::device_vector<float> cost_device_vec_for_sorting( N_OF_SAMPLES );
    //thrust::device_vector<float> cost_vec_dev_temp( N_OF_SAMPLES );
    //cost_device_vec_for_sorting = cost_vec_dev_temp;
    
    Input_vec *d_Input_vec;
    //Input_vec *h_Input_vec;
    //h_Input_vec = (Input_vec *)malloc(sizeof(Input_vec) * N_OF_SAMPLES);
    cudaMalloc(&d_Input_vec, sizeof(Input_vec) * N_OF_SAMPLES);
    set_Input_vec<<<N_OF_SAMPLES,1>>>(d_Input_vec, 0.0f);

    /*for(int i = 0; i < N_OF_SAMPLES; i++){
        for(int k = 0; k < HORIZON; k++){
            h_Input_vec[i].Input[k] = 0.0f;
        }
    }*/
#endif
    /*Data1 *h_dataFromBlocks;
    Data1 *d_dataFromBlocks;*/
    h_dataFromBlocks = (Data1 *)malloc(sizeof(Data1)*numBlocks);
    cudaMalloc(&d_dataFromBlocks, sizeof(Data1) * numBlocks);



    /* curand の設定 */
    curandState *devStates;
    cudaMalloc((void **)&devStates, randomNums * sizeof(curandState));
    setup_kernel<<<N_OF_SAMPLES * (dim_U+1), HORIZON>>>(devStates,rand());
    cudaDeviceSynchronize();

    /* Covariance の定義 */
    float *h_hat_Q, *Diag_D;
    float *device_cov;
    float *device_diag_eig = NULL;
    float *d_hat_Q;
    h_hat_Q = (float *)malloc(sizeof(float)*dim_hat_Q);
    Diag_D = (float *)malloc(sizeof(float)*dim_hat_Q);
    cudaMalloc(&device_cov, sizeof(float)*dim_hat_Q);
    cudaMalloc(&device_diag_eig, sizeof(float)*dim_hat_Q);
    cudaMalloc(&d_hat_Q, sizeof(float)*dim_hat_Q);

    setup_init_Covariance<<<HORIZON, HORIZON>>>(d_hat_Q);

    /* 準最適制御入力列 */
    float *Us_host, *Us_device;
    Us_host = (float *)malloc(sizeof(float) * HORIZON);
    for(int i = 0; i < HORIZON; i++){
        Us_host[i] = 0.0f;
    }
    cudaMalloc(&Us_device, sizeof(float) * HORIZON);


    float var;
    float now_u;
    for(int i = 0; i < Blocks; i++){
        for(int k = 0; k < HORIZON; k++){
            h_dataFromBlocks[i].Input[k] = 0.0f;
        }
    }

    /* 固有値の取得 */
    

    const int m = HORIZON;
    const int lda = m;

    float eig_vec[m] = { };

    float *d_A;
    float *d_W;
    int *devInfo;
    float *d_work;
    int lwork = 0;

    // int work_size;
    // float *work_space;

    int info_gpu = 0;

    cudaStat1 = cudaMalloc ((void**)&d_A, sizeof(float) * lda * m);
    cudaStat2 = cudaMalloc ((void**)&d_W, sizeof(float) * m);
    cudaStat3 = cudaMalloc ((void**)&devInfo, sizeof(int));
    assert(cudaSuccess == cudaStat1);
    assert(cudaSuccess == cudaStat2);
    assert(cudaSuccess == cudaStat3);
    cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; // compute eigenvalues and eigenvectors.
    cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;



    for(int time = 0; time < TIME; time++){
        for(int repeat = 0; repeat < Recalc; repeat++){
            var = Variavility * pow(0.8,repeat);
            //var = Variavility;
            cudaMemcpy(d_dataFromBlocks, h_dataFromBlocks, sizeof(Data1)*numBlocks, cudaMemcpyHostToDevice);
            cudaDeviceSynchronize();
            // MCMPC_GPU<<<numBlocks, THREAD_PER_BLOCKS>>>(state, devStates, d_dataFromBlocks, var, Blocks, d_hat_Q);
#ifdef USING_THRUST
            Using_Thrust_MCMPC_Linear<<<numBlocks, THREAD_PER_BLOCKS>>>(state[0],state[1],state[2],devStates, d_Input_vec, var, Blocks, d_hat_Q, device_param, device_matrix, thrust::raw_pointer_cast( cost_device_vec_for_sorting.data() ));
            thrust::sequence( indices_device_vec.begin(), indices_device_vec.end() );
            thrust::sort_by_key( cost_device_vec_for_sorting.begin(), cost_device_vec_for_sorting.end(), indices_device_vec.begin() );
            callback_elite_sample<<<Blocks,1>>>(d_dataFromBlocks, d_Input_vec, thrust::raw_pointer_cast( indices_device_vec.data() ));
#else
            MCMPC_GPU_Linear_Example<<<numBlocks, THREAD_PER_BLOCKS>>>(state[0],state[1],state[2], devStates, d_dataFromBlocks, var, Blocks, d_hat_Q, device_param, device_matrix);
            cudaDeviceSynchronize();
            //cudaMemcpy(h_dataFromBlocks, d_dataFromBlocks, sizeof(Data1) * numBlocks, cudaMemcpyDeviceToHost);
#endif
            cudaMemcpy(h_dataFromBlocks, d_dataFromBlocks, sizeof(Data1) * numBlocks, cudaMemcpyDeviceToHost);
            printf("TOP  W == %f WORST W == %f\n",h_dataFromBlocks[0].W,  h_dataFromBlocks[Blocks-1].W);
            weighted_mean(h_dataFromBlocks, Blocks, Us_host);
            //printMatrix(m,1,Us_host, m, "u");
            cudaMemcpy(Us_device, Us_host, sizeof(float) * HORIZON, cudaMemcpyHostToDevice);
            //printf("hoge\n");
            calc_Var_Cov_matrix<<<HORIZON, HORIZON>>>(device_cov, d_dataFromBlocks, Us_device, Blocks);
            cudaDeviceSynchronize();
            cudaMemcpy(h_hat_Q, device_cov, sizeof(float)*dim_hat_Q, cudaMemcpyDeviceToHost);
            //printMatrix(m,m,h_hat_Q, lda, "DA");
            
            //cudaStat1 = cudaMemcpy(d_A, h_hat_Q, sizeof(float) * lda * m, cudaMemcpyHostToDevice);
            cudaStat1 = cudaMemcpy(d_A, h_hat_Q, sizeof(float) * lda * m, cudaMemcpyHostToDevice);
            assert(cudaSuccess == cudaStat1);
            cusolver_status = cusolverDnSsyevd_bufferSize(
                cusolverH,
                jobz,
                uplo,
                m,
                d_A,
                lda,
                d_W,
                &lwork);
            assert (cusolver_status == CUSOLVER_STATUS_SUCCESS);
            //cudaMemcpy(Diag_D, d_A, sizeof(float)*lda*m, cudaMemcpyDeviceToHost);
            //printMatrix(m,m,Diag_D, m, "V");
            

            cudaStat1 = cudaMalloc((void**)&d_work, sizeof(float)*lwork);
            //assert(cudaSuccess == cudaStat1);

            cusolver_status = cusolverDnSsyevd(
                cusolverH,
                jobz,
                uplo,
                m,
                d_A,
                lda,
                d_W,
                d_work,
                lwork,
                devInfo);

            cudaDeviceSynchronize();
            //assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
            //assert(cudaSuccess == cudaStat1);

            cudaStat1 = cudaMemcpy(eig_vec, d_W, sizeof(float)*m, cudaMemcpyDeviceToHost);
            cudaStat2 = cudaMemcpy(Diag_D, d_A, sizeof(float)*lda*m, cudaMemcpyDeviceToHost);
            cudaStat3 = cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost);
            //assert(cudaSuccess == cudaStat1);
            //assert(cudaSuccess == cudaStat2);
            //assert(cudaSuccess == cudaStat3);
            // printMatrix(m,1,eig_vec, lda, "C");
            // printf("=====Upper is eigen value====");
            //printMatrix(m,m,Diag_D, lda, "C");
            make_Diagonalization<<<HORIZON,HORIZON>>>(d_W, d_A);
            cudaMemcpy(h_hat_Q, d_A, sizeof(float)*lda*m, cudaMemcpyDeviceToHost);
            //printMatrix(m,m,h_hat_Q, lda, "C");
            cudaMemcpy(device_diag_eig, h_hat_Q, sizeof(float)*dim_hat_Q, cudaMemcpyHostToDevice);
            cudaMemcpy(d_hat_Q, Diag_D, sizeof(float)*dim_hat_Q, cudaMemcpyHostToDevice);
            tanspose<<<HORIZON,HORIZON>>>(device_cov, d_hat_Q);
            pwr_matrix_answerB<<<HORIZON,HORIZON>>>(device_cov, device_diag_eig);
            cudaDeviceSynchronize();

            /*cusolver_status = cusolverDnSpotrf_bufferSize(cusolverH, uplo, m, device_cov, m, &work_size);
            assert( cusolver_status == CUSOLVER_STATUS_SUCCESS );
            // float* workspace ;
            cudaMalloc((void**)&work_space, sizeof(float)*work_size);
            cusolver_status = cusolverDnSpotrf(cusolverH, uplo, m, device_cov, m , work_space, work_size, devInfo);
            assert( cusolver_status == CUSOLVER_STATUS_SUCCESS );
            setup_init_Covariance<<<HORIZON, HORIZON>>>(d_hat_Q);
            cusolver_status = cusolverDnSpotrs(cusolverH, uplo, m, m , device_cov, m, d_hat_Q, m, devInfo);
            assert( cusolver_status == CUSOLVER_STATUS_SUCCESS );*/
            tanspose<<<HORIZON,HORIZON>>>(d_hat_Q, device_cov);
            // pwr_matrix_answerA<<<HORIZON,HORIZON>>>(device_diag_eig, device_cov);
            pwr_matrix_answerA<<<HORIZON,HORIZON>>>(device_diag_eig, d_hat_Q);
            cudaDeviceSynchronize();
            tanspose<<<HORIZON,HORIZON>>>(d_hat_Q, device_diag_eig);
            //cudaMemcpy(h_hat_Q, device_diag_eig, sizeof(float)*dim_hat_Q, cudaMemcpyDeviceToHost);
            cudaMemcpy(h_hat_Q, d_hat_Q, sizeof(float)*dim_hat_Q, cudaMemcpyDeviceToHost);
            //cudaMemcpy(d_hat_Q, h_hat_Q, sizeof(float)*dim_hat_Q, cudaMemcpyHostToDevice);
            printMatrix(m,m,h_hat_Q, lda, "C");

            fprintf(fp,"%f %f %f %f %f %f %f %f %f %f\n",Us_host[0], Us_host[1],
                    Us_host[2], Us_host[3], Us_host[4], Us_host[5], Us_host[6], Us_host[7], Us_host[8], Us_host[9]);

            for(int count = 0; count < HORIZON; count++){
                h_dataFromBlocks[0].Input[count] = Us_host[count];
            }
#ifdef USING_THRUST
                reset_Input_vec<<<numBlocks,THREAD_PER_BLOCKS>>>(d_Input_vec, Us_device);
                cudaDeviceSynchronize();
#endif

#ifdef Linear
            float RSME=0.0f;
            for(int d = 0; d < HORIZON; d++){
                Error[d] = Us_host[d] - opt[d];
                RSME += powf(Error[d],2);
            }
            printf("RSME == %f\n", RSME / HORIZON);
#endif
        }
        //printMatrix(m,m,h_hat_Q, lda, "C");
        now_u = Us_host[0];
        calc_Linear_example(state, now_u, params, state);
        for(int i = 0; i < Blocks; i++){
            for(int k = 0; k < HORIZON - 1; k++){
                h_dataFromBlocks[i].Input[k] = Us_host[k+1];
            }
            h_dataFromBlocks[i].Input[HORIZON-1] = Us_host[HORIZON - 1];
        }
    }
    if (d_A    ) cudaFree(d_A);
    if (d_W    ) cudaFree(d_W);
    if (devInfo) cudaFree(devInfo);
    if (d_work ) cudaFree(d_work);
    //cudaFree(indices_device_vec);
    //cudaFree(cost_device_vec_for_sorting);

    if (cusolverH) cusolverDnDestroy(cusolverH);
    fclose(fp);
    // fclose(hp);
    //thrust::reduce(indices_device_vec,cost_device_vec_for_sorting);
    cudaDeviceReset();
    printf("%s\n", cudaGetErrorString(cudaGetLastError()));
    return 0;
}
