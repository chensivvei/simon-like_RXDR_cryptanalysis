# simon-like_RXDR_cryptanalysis
This project is related to the paper "rotational-XOR differential ractangle cryptanalysis on Simon-like ciphers".
There are four files in total. The file "multiple_differential.cpp" is used to calculate the probability of multiple differential corresponding to the paper's Algorithm 2. The files "simon32.cu" and "simeck32.cu" are CUDA files and correspond to the paper's Algorithm 1, which are used to verify the RXDR distinguishers of Simon32/64 and Simeck32.
