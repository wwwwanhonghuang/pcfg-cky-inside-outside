// main.cu
#include <iostream>
#include <cuda_runtime.h>

#define N 1000
#define MAX_SEQUENCE_SIZE 512  


__global__ void dp_task(float* d_data, int max_size) {
    int i = blockIdx.x;  // row index
    int j = threadIdx.x + blockIdx.y * blockDim.x;  // column index

    if (i < max_size && j < max_size) {
        for (int k = i; k < j; k++) {
            // DP relation: dp[i, j] += dp[i, k] + dp[k + 1, j]
            d_data[i * MAX_SEQUENCE_SIZE + j] += 
                d_data[i * MAX_SEQUENCE_SIZE + k] + d_data[(k + 1) * MAX_SEQUENCE_SIZE + j];
        }
    }
}

int main() {

    // host data (data preparation)
    float *h_data = (float*)malloc(N * MAX_SEQUENCE_SIZE * sizeof(float));

    float *d_data;
    cudaMalloc((void**)&d_data, N * MAX_SEQUENCE_SIZE * sizeof(float));
    cudaMemcpy(d_data, h_data, N * MAX_SEQUENCE_SIZE * sizeof(float), cudaMemcpyHostToDevice);

    int blockSize = 512;
    int numBlocks = (int)((MAX_SEQUENCE_SIZE + blockSize - 1) / blockSize);
    
    dp_task<<<numBlocks, blockSize>>>(d_data, MAX_SEQUENCE_SIZE);
  


    return 0;
}
