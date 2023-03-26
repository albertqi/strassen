#include <stdio.h>
#include <unordered_set>
#include <string>
#include <fstream>
#include <random>
#include <chrono>
#include <thread>

#define N0 68
#define NUM_THREADS 10

long zero = 0;
thread_local std::mt19937 generator;
std::mutex m;
std::unordered_set<long *> pointers;

class Matrix
{
private:
    const int original_dim, dim, row_start, col_start;
    long *arr;

    Matrix(const int &original_dim, const int &dim, const int &row_start, const int &col_start, long *arr)
        : original_dim(original_dim), dim(dim), row_start(row_start), col_start(col_start), arr(arr)
    {
        m.lock();
        pointers.insert(arr);
        m.unlock();
    }

    const Matrix quarter(const int &q) const
    {
        const int new_dim = (dim + 1) / 2;
        const int new_row = q >= 2 ? row_start + new_dim : row_start,
                  new_col = q % 2 == 1 ? col_start + new_dim : col_start;
        return Matrix(original_dim, new_dim, new_row, new_col, arr);
    }

    static const Matrix shift(const Matrix &X, const Matrix &Y, const int &delta)
    {
        // Handle mismatched dimensions
        if (X.original_dim != Y.original_dim || X.dim != Y.dim)
        {
            throw std::invalid_argument("Mismatched dimensions");
        }

        const int dim = X.dim;
        const Matrix Z(dim, new long[dim * dim]);

        for (int i = 0; i < dim; ++i)
        {
            for (int j = 0; j < dim; ++j)
            {
                Z(i, j) = X(i, j) + delta * Y(i, j);
            }
        }

        return Z;
    }

    static const Matrix add(const Matrix &X, const Matrix &Y)
    {
        return shift(X, Y, 1);
    }

    static const Matrix subtract(const Matrix &X, const Matrix &Y)
    {
        return shift(X, Y, -1);
    }

    long &operator()(const int &row, const int &col) const
    {
        if (row_start + row >= original_dim || col_start + col >= original_dim)
        {
            return zero;
        }

        return arr[(row_start + row) * original_dim + col_start + col];
    }

public:
    Matrix(const int &dim, long *arr)
        : original_dim(dim), dim(dim), row_start(0), col_start(0), arr(arr)
    {
        m.lock();
        pointers.insert(arr);
        m.unlock();
    }

    void print_diagonal() const
    {
        for (int i = 0; i < dim; ++i)
        {
            printf("%ld\n", (*this)(i, i));
        }
    }

    const long sum_diagonal() const
    {
        long x = 0;
        for (int i = 0; i < dim; ++i)
        {
            x += (*this)(i, i);
        }
        return x;
    }

    static const Matrix multiply_conventional(const Matrix &X, const Matrix &Y)
    {
        // Handle mismatched dimensions
        if (X.dim != Y.dim)
        {
            throw std::invalid_argument("Mismatched dimensions");
        }

        const int dim = X.dim;
        const Matrix Z(dim, new long[dim * dim]);

        for (int i = 0; i < dim; ++i)
        {
            for (int j = 0; j < dim; ++j)
            {
                for (int k = 0; k < dim; ++k)
                {
                    Z(i, j) += X(i, k) * Y(k, j);
                }
            }
        }

        return Z;
    }

    static const Matrix multiply_strassen_threshold(const Matrix &X,
                                                    const Matrix &Y,
                                                    const int &threshold)
    {
        // Handle mismatched dimensions
        if (X.dim != Y.dim)
        {
            throw std::invalid_argument("Mismatched dimensions");
        }

        const int dim = X.dim;

        if (dim <= threshold)
        {
            return multiply_conventional(X, Y);
        }

        const Matrix A = X.quarter(0),
                     B = X.quarter(1),
                     C = X.quarter(2),
                     D = X.quarter(3),
                     E = Y.quarter(0),
                     F = Y.quarter(1),
                     G = Y.quarter(2),
                     H = Y.quarter(3);

        const Matrix P1 = multiply_strassen(A, subtract(F, H)),
                     P2 = multiply_strassen(add(A, B), H),
                     P3 = multiply_strassen(add(C, D), E),
                     P4 = multiply_strassen(D, subtract(G, E)),
                     P5 = multiply_strassen(add(A, D), add(E, H)),
                     P6 = multiply_strassen(subtract(B, D), add(G, H)),
                     P7 = multiply_strassen(subtract(C, A), add(E, F));

        const Matrix Q1 = add(P4, add(P5, subtract(P6, P2))),
                     Q2 = add(P1, P2),
                     Q3 = add(P3, P4),
                     Q4 = add(P1, add(P5, subtract(P7, P3)));

        const Matrix Z(dim, new long[dim * dim]);
        const bool should_trim = Q1.dim * 2 != dim;

        for (int i = 0; i < (dim + 1) / 2; ++i)
        {
            for (int j = 0; j < (dim + 1) / 2; ++j)
            {
                Z(i, j) = Q1(i, j);
                if (!should_trim || j != (dim - 1) / 2)
                {
                    Z(i, j + (dim + 1) / 2) = Q2(i, j);
                }
                if (!should_trim || i != (dim - 1) / 2)
                {
                    Z(i + (dim + 1) / 2, j) = Q3(i, j);
                }
                if (!should_trim || (i != (dim - 1) / 2) && j != (dim - 1) / 2)
                {
                    Z(i + (dim + 1) / 2, j + (dim + 1) / 2) = Q4(i, j);
                }
            }
        }

        return Z;
    }

    static const Matrix multiply_strassen(const Matrix &X, const Matrix &Y)
    {
        return multiply_strassen_threshold(X, Y, N0);
    }
};

int main(int argc, char *argv[])
{
    // Check for valid command line arguments
    if (argc != 4)
    {
        printf("Usage: ./strassen [flag] [dimension] [input_file]\n");
        return 1;
    }

    // Parse command line arguments
    const double flag = atof(argv[1]);
    const int dim = atoi(argv[2]);
    const std::string input_file = argv[3];

    // Count triangles in random graph
    if (flag > 0 && flag < 1)
    {
        long total_triangles = 0;

        std::vector<std::thread> threads;
        for (int i = 0; i < NUM_THREADS; ++i)
        {
            threads.push_back(std::thread([&flag, &total_triangles]()
            {
                // Set up random generator
                std::random_device random_dev;
                generator.seed(random_dev());
                std::uniform_real_distribution<double> unif(0.0, 1.0);

                long *x = new long[1024 * 1024];
                for (int i = 0; i < 1024; ++i)
                {
                    for (int j = 0; j < i; ++j)
                    {
                        const bool include_edge = unif(generator) <= flag;

                        x[i * 1024 + j] = include_edge;
                        x[j * 1024 + i] = include_edge;
                    }
                }

                // Create matrix `X`
                const Matrix X(1024, x);

                // Calculate product `Y`
                const Matrix Y = Matrix::multiply_strassen(X, Matrix::multiply_strassen(X, X));

                m.lock();
                total_triangles += Y.sum_diagonal() / 6;
                m.unlock();
            }));
        }
        for (std::thread &t : threads)
        {
            t.join();
        }

        // Print results
        printf("%f\n", (double)total_triangles / NUM_THREADS);

        // Free dynamically allocated memory
        for (const auto &p : pointers)
        {
            delete[] p;
        }

        return 0;
    }

    // Open input file
    std::ifstream file(input_file);

    // Parse input file
    long *a = new long[dim * dim], *b = new long[dim * dim];
    for (int i = 0; i < dim * dim; ++i)
    {
        file >> a[i];
    }
    for (int i = 0; i < dim * dim; ++i)
    {
        file >> b[i];
    }

    // Create matrices `A` and `B`
    const Matrix A(dim, a), B(dim, b);

    // Find best threshold for `multiply_strassen_threshold`
    if (flag == 1)
    {
        // Keep track of `min_time` and `best_threshold`
        double min_time = INT_MAX;
        int best_threshold;

        // Calculate product `C` for various values of `threshold`
        for (int threshold = 1; threshold <= 100; ++threshold)
        {
            double total_runtime = 0;

            std::vector<std::thread> threads;
            for (int i = 0; i < NUM_THREADS; ++i)
            {
                threads.push_back(std::thread([&A, &B, &threshold, &total_runtime]()
                {
                    const auto start = std::chrono::high_resolution_clock::now();
                    const Matrix C = Matrix::multiply_strassen_threshold(A, B, threshold);
                    const auto end = std::chrono::high_resolution_clock::now();
                    const std::chrono::duration<double> diff = end - start;
                    
                    m.lock();
                    total_runtime += diff.count();
                    m.unlock();
                }));
            }
            for (std::thread &t : threads)
            {
                t.join();
            }

            const double runtime = total_runtime / NUM_THREADS;

            // Print results
            printf("%f seconds with threshold of %d\n", runtime, threshold);

            // Update `min_time` and `best_threshold`
            if (runtime < min_time)
            {
                min_time = runtime;
                best_threshold = threshold;
            }

            // Free dynamically allocated memory
            for (const auto &p : pointers)
            {
                if (p != a && p != b)
                {
                    delete[] p;
                }
            }
            pointers = std::unordered_set<long *>({a, b});
        }

        // Print `best_threshold`
        printf("%d\n", best_threshold);

        // Free dynamically allocated memory
        for (const auto &p : pointers)
        {
            delete[] p;
        }

        return 0;
    }

    // Calculate product `C`
    const Matrix C = Matrix::multiply_strassen(A, B);

    // Print results
    C.print_diagonal();

    // Free dynamically allocated memory
    for (const auto &p : pointers)
    {
        delete[] p;
    }

    return 0;
}