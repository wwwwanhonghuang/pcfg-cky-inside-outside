// main.cu
#include <iostream>
#include <cuda_runtime.h>

#define N 1000
#define MAX_SEQUENCE_SIZE 512  

// CUDA kernel for vector addition
__global__ void dp_task(int i, int j) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    return;
}

int main() {

    // host data (data preparation)
    float *h_data = (float*)malloc(N * MAX_SEQUENCE_SIZE * sizeof(float));



    float *d_data;
    cudaMalloc((void**)&d_data, N * MAX_SEQUENCE_SIZE * sizeof(float));
    cudaMemcpy(d_data, h_data, N * MAX_SEQUENCE_SIZE * sizeof(float), cudaMemcpyHostToDevice);




    return 0;
}
