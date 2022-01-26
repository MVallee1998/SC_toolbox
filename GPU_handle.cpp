#include <iostream>
#include <cuda.h>
using namespace std;

__shared__ int r[2688];
__device__ unsigned long ai[6][256];
__device__ unsigned short mi[6][5][256];
__global__ void kernel(){
    unsigned long a[6];
    unsigned short m[6][5];
    bool stop;
    for (int i=0;i<2688;i++){
        unsigned long x=X[i];
        for (int j=0; j<6; j++){
        if __popc(x & v[i])
        }
    }
}
