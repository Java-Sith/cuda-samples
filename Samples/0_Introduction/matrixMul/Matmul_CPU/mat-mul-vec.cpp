#include <iostream>
#include <cstdlib>
#include <ctime>
#include <chrono>

// Function to generate a matrix with given rows and columns
void generateMatrix(float* matrix, int rows, int cols) {
    srand(static_cast<unsigned>(time(0))); // Seed for random number generation

    for (int i = 0; i < rows * cols; ++i) {
        matrix[i] = static_cast<float>(rand() % 1000) / 100.0f; // Generate random float values between 0.0 and 9.99
    }
}

static void vec_dot(const int n, float * s, float * x, float * y) {
    float sumf = 0.0;
    for (int i = 0; i < n; ++i) {
        sumf += x[i] * y[i];
    }
    *s = sumf;
}

// Function to multiply two matrices
static void gemm(int m, int n, int k,
                 float * A,
                 float * B,
                 float * C,
                 const int ith, const int nth) {
    int m0, m1, n0, n1;
    if (m > n) {
        n0 = 0;
        n1 = n;
        const int np = m;
        const int dp = (np + nth - 1)/nth;
        m0 = dp*ith;
        m1 = std::min(m0 + dp, np);
    } else {
        m0 = 0;
        m1 = m;
        const int np = n;
        const int dp = (np + nth - 1)/nth;
        n0 = dp*ith;
        n1 = std::min(n0 + dp, np);
    }

    for (int i = m0; i < m1; i++) {
        for (int j = n0; j < n1; j++) {
            vec_dot(k, C + i * n + j, A + i * k, B + j * k);
        }
    }
}

void perform_gemm_test(float* a, float* b, int M, int N, int K) {
    printf("\nPerforming GEMM test:\n");

    auto start = std::chrono::high_resolution_clock::now();

    float* gemm_out = new float[M * N];
    gemm(M, N, K, a, b, gemm_out, 0, 1);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Time in seconds: " << elapsed_seconds.count() << " s\n";
    std::cout << "Time in milliseconds: " << elapsed_milliseconds.count() << " ms\n";
}


int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <M> <N> <K>" << std::endl;
        return 1;
    }

    int M = std::atoi(argv[1]);
    int N = std::atoi(argv[2]);
    int K = std::atoi(argv[3]);

    // Allocate memory for the matrices
    float* matrixA = new float[M * K];
    float* matrixB = new float[K * N];
    float* matrixC = new float[M * N];

    // Generate random matrices
    generateMatrix(matrixA, M, K);
    generateMatrix(matrixB, K, N);

    perform_gemm_test(matrixA, matrixB, M, N, K);

    // Free allocated memory
    delete[] matrixA;
    delete[] matrixB;
    delete[] matrixC;

    return 0;
}
