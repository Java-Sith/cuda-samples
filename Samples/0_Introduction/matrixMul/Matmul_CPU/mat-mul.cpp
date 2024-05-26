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

// Function to multiply two matrices
void multiplyMatrices(const float* matrixA, const float* matrixB, float* matrixC, int M, int N, int K) {
    // Initialize the result matrix with zeros
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            matrixC[i * N + j] = 0.0f;
        }
    }

    // Perform matrix multiplication
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < K; ++k) {
                matrixC[i * N + j] += matrixA[i * K + k] * matrixB[k * N + j];
            }
        }
    }
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

    // Measure the time taken for matrix multiplication
    auto start = std::chrono::high_resolution_clock::now();
    multiplyMatrices(matrixA, matrixB, matrixC, M, N, K);
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration in milliseconds
    std::chrono::duration<float, std::milli> duration = end - start;

    // Print the time taken for multiplication
    std::cout << "Time taken for matrix multiplication: " << duration.count() << " milliseconds" << std::endl;

    // Free allocated memory
    delete[] matrixA;
    delete[] matrixB;
    delete[] matrixC;

    return 0;
}

