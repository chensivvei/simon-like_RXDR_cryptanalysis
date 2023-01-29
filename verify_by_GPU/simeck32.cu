#include <iostream>
#include <cstdint>
#include <random>
#include "common.h"

#define KEYSIZE 1024

#define ROUND 15

#if ROUND == 14 ||ROUND == 15

#define ROUND_CONST 0x0004

#endif

#if ROUND == 16

#define ROUND_CONST 0x0006

#endif




#define ROTL(x, n) ((((uint16_t)(x)<<(n)) | (uint16_t)(x)>>(16-(n))))
#define ROTR(x, n) ((((uint16_t)(x)>>(n)) | (uint16_t)(x)<<(16-(n))))


__global__ void encrypt(uint32_t *data1) {
    uint32_t tid = blockDim.x * blockIdx.x + threadIdx.x;
    uint16_t right = tid;
    uint16_t left = tid >> 16;
    uint16_t tmp;
    uint32_t key[11] = {0x94f6, 0x3b02, 0xb740, 0x1a89, 0xff71, 0x806, 0x7ef0, 0xed8c, 0x5fe0, 0xb870, 0x6ce};
    for (int r = 0; r < 11; r++) {
        tmp = left;
        left = (ROTL(left, 8) & ROTL(left, 1)) ^ ROTL(left, 2) ^ right ^ key[r];
        right = tmp;
    }
    data1[tid] = ((uint32_t) left << 16) ^ right;
}

__global__ void encrypt2(uint32_t *data2) {
    uint32_t tid = blockDim.x * blockIdx.x + threadIdx.x;
    uint16_t right;
    uint16_t left;
    uint16_t tmp;

    uint32_t key2[11] = {0x29eb, 0x7604, 0x6e81, 0x3512, 0xfee3, 0x100a, 0x5de6, 0x651f, 0x3e3, 0xcca1, 0x50d5};
    right = tid;
    right = ROTL(right, 1) ^ 0x0006;
    left = tid >> 16;
    left = ROTL(left, 1);
    for (int r = 0; r < 11; r++) {
        tmp = left;
        left = (ROTL(left, 8) & ROTL(left, 1)) ^ ROTL(left, 2) ^ right ^ key2[r];
        right = tmp;
    }
    data2[tid] = ((uint32_t) left << 16) ^ right;
}

__global__ void define(const uint32_t *data1, const uint32_t *data2, uint32_t *count) {
    uint64_t tid = blockDim.x * blockIdx.x + threadIdx.x;
    uint32_t diff;
    for (uint64_t i = 0; i < (1 << 7); i++) {
        if (tid > i) {
            diff = data1[tid] ^ data1[i];
            if (diff == 0x20220008 && ((data2[tid] ^ data2[i]) == 0x40440010)) {
                atomicAdd(count, 1);
            }
        }
    }
}

__global__ void define11(const uint32_t *data1, const uint32_t *data2, uint32_t *count) {
    uint64_t tid = blockDim.x * blockIdx.x + threadIdx.x;
    uint32_t diff;
    for (uint64_t i = 0; i < (1 << 21); i++) {
        if (tid > i) {
            diff = data1[tid] ^ data1[i];
            if (diff == 0x20220008 && ((data2[tid] ^ data2[i]) == 0x40440010)) {
                atomicAdd(count, 1);
            }
        }
    }
}

__global__ void define12(const uint32_t *data1, const uint32_t *data2, uint32_t *count) {
    uint64_t tid = blockDim.x * blockIdx.x + threadIdx.x;
    uint32_t diff;
    for (uint64_t i = 0; i < (1 << 10); i++) {
        if (tid > i) {
            diff = data1[tid] ^ data1[i];
            if (diff == 0x20220008 && ((data2[tid] ^ data2[i]) == 0x40440010)) {
                atomicAdd(count, 1);
            }
        }
    }
}



__global__ void define3(uint32_t *count, uint16_t *keylist) {
    uint32_t tid = blockDim.x * blockIdx.x + threadIdx.x;
    uint16_t right;
    uint16_t left;
    uint16_t tmp;

    //decrypt c1
    for (int i = 0; i < KEYSIZE; i++) {
        int key_index = i * ROUND * 2;
        uint16_t *key = keylist + key_index;
        uint16_t *key2 = keylist + key_index + ROUND;

        right = tid >> 16;
        left = tid;
        for (int r = 0; r < ROUND; r++) {
            tmp = left;
            left = (ROTL(left, 5) & left) ^ ROTL(left, 1) ^ right ^ key[ROUND - 1 - r];
            right = tmp;
        }
        uint32_t p1 = ((uint32_t) right << 16) ^ left;

        //decrypt c3

        uint32_t c3 = tid ^ 0x00150008;

        right = c3 >> 16;
        left = c3;
        for (int r = 0; r < ROUND; r++) {
            tmp = left;
            left = (ROTL(left, 5) & left) ^ ROTL(left, 1) ^ right ^ key[ROUND - 1 - r];
            right = tmp;
        }
        uint32_t p3 = ((uint32_t) right << 16) ^ left;

        right = p1;
        left = p1 >> 16;

        //P2 = P1 <<< 1 ^ 0x00000004


        uint32_t p2 = ((uint32_t) ROTL(left, 1) << 16) ^ ROTL(right, 1) ^ ROUND_CONST;




        right = p3;
        left = p3 >> 16;



        uint32_t p4 = ((uint32_t) ROTL(left, 1) << 16) ^ ROTL(right, 1) ^ ROUND_CONST;





        //encrypt p2
        right = p2;
        left = p2 >> 16;
        for (int r = 0; r < ROUND; r++) {
            tmp = left;
            left = (ROTL(left, 5) & left) ^ ROTL(left, 1) ^ right ^ key2[r];
            right = tmp;
        }
        uint32_t c2 = ((uint32_t) left << 16) ^ right;

        //encrypt p4
        right = p4;
        left = p4 >> 16;
        for (int r = 0; r < ROUND; r++) {
            tmp = left;
            left = (ROTL(left, 5) & left) ^ ROTL(left, 1) ^ right ^ key2[r];
            right = tmp;
        }
        uint32_t c4 = ((uint32_t) left << 16) ^ right;

        if ((c4 ^ c2) == 0x002a0010) {
            atomicAdd(count + i, 1);
//            if (tid < 0xffffff){
//                printf("c1 %#x c2 %#x c3 %#x c4 %#x p1 %#x p2 %#x p3 %#x p4 %#x\n",tid, c2,c3,c4,p1,p2,p3,p4);
//            }

        }

    }
}

void keyschedule(uint16_t *key, int R) {
    uint16_t c;
    for (int i = 0; i < R - 4; i++) {

        int r = i + 1;

        if ((r == 1) || (r == 2) || (r == 3) || (r == 4)
            || (r == 5) || (r == 9) || (r == 10) ||
            (r == 12) || (r == 13) || (r == 14) ||
            (r == 16) || (r == 18) || (r == 23) ||
            (r == 26) || (r == 28) || (r == 29) || (r == 32)) {
            c = 0xfffd;
        } else {
            c = 0xfffc;
        }

        key[i + 4] = (ROTL(key[i + 1], 5) & key[i + 1]) ^ ROTL(key[i + 1], 1) ^ c;
    }
}

void init_key(uint16_t *key) {
    using std::default_random_engine;
    default_random_engine e;
    e.seed(time(NULL));
    for (int i = 0; i < KEYSIZE; i++) {
        int index = i * ROUND * 2;
        key[index] = e();
        key[index + 1] = e();
        key[index + 2] = e();
        key[index + 3] = e();
        keyschedule(key + index, ROUND);
        key[index + ROUND] = ROTL(key[index], 1) ^ ROUND_CONST;
        key[index + ROUND + 1] = ROTL(key[index + 1], 1);
        key[index + ROUND + 2] = ROTL(key[index + 2], 1);
        key[index + ROUND + 3] = ROTL(key[index + 3], 1);
        keyschedule(key + index + ROUND, ROUND);
    }
}


int main() {
//    printf("%#x \n %#x\n", ROTL(0x2022,1), ROTL(0x0008,1));//0x4044   0x10

    cudaSetDevice(3);

    size_t nWords = 1ULL << 32;

    uint32_t *count;
    uint16_t *key;

    CHECK(cudaMallocManaged((uint32_t **) &count, KEYSIZE * sizeof(uint32_t)));
    CHECK(cudaMallocManaged((uint32_t **) &key, KEYSIZE * sizeof(uint16_t) * ROUND * 2));

    memset(count, 0, KEYSIZE * sizeof(uint32_t));

    init_key(key);

//    printf("key: ");
//    for (int j = 0; j < 11; j++) {
//        printf("%#x ", key[j]);
//    }
//    printf("\n");
//    printf("key2: ");
//    for (int j = 0; j < 11; j++) {
//        printf("%#x ", key[11 + j]);
//    }
//    printf("\n");

    int blockx = 1024;
    dim3 block = blockx;
    dim3 grid = nWords / blockx;

    define3<<<grid, block>>>(count, key);

    CHECK(cudaDeviceSynchronize());

//    for (int i = 0; i < KEYSIZE; i++) {
//        for (int j = 0; j < 11; j++) {
//            printf("%#x ",key[j + i * 22]);
//        }
//        printf("\n");
//        for (int j = 0; j < 11; j++) {
//            printf("%#x ",key[j + 11 + i * 22]);
//        }
//        printf("\n");
//        printf("\n");
//    }


    for (int i = 0; i < KEYSIZE; i++) {
        printf("count %d is %d\n", i, count[i]);
    }

    uint64_t problity = 0;
    for (int i = 0; i < KEYSIZE; i++) {
        problity += count[i];
    }
    printf("\nThe probability is %f\n", (double) problity / KEYSIZE);

//    cudaFree(d_data1);
//    cudaFree(d_data2);
//    free(h_data1);
//    free(h_data2);
    return 0;
}
