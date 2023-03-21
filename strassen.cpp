#include <stdio.h>
#include <unordered_set>
#include <string>
#include <fstream>
#include <random>
#include <chrono>

#define N0 68

long zero = 0;
std::unordered_set<long *> pointers;
std::mt19937 generator;

class Matrix
{
private:
    const int original_dim, dim, row_start, col_start;
    long *arr;

    Matrix(const int &original_dim, const int &dim, const int &row_start, const int &col_start, long *arr)
        : original_dim(original_dim), dim(dim), row_start(row_start), col_start(col_start), arr(arr)
    {
        pointers.insert(arr);
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
        pointers.insert(arr);
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

    static const Matrix multiply_strassen(const Matrix &X, const Matrix &Y)
    {
        // Handle mismatched dimensions
        if (X.dim != Y.dim)
        {
            throw std::invalid_argument("Mismatched dimensions");
        }

        const int dim = X.dim;

        if (dim <= N0)
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

        // Print results
        printf("%ld\n", Y.sum_diagonal() / 6);

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

    // Find runtimes for `multiply_conventional` and `multiply_strassen`
    if (flag == 1)
    {
        // Calculate product `C`
        auto start = std::chrono::high_resolution_clock::now();
        const Matrix C = Matrix::multiply_conventional(A, B);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        const double conventional_time = diff.count();

        // Calculate product `D`
        start = std::chrono::high_resolution_clock::now();
        const Matrix D = Matrix::multiply_strassen(A, B);
        end = std::chrono::high_resolution_clock::now();
        diff = end - start;
        const double strassen_time = diff.count();

        // Print results
        printf("%f seconds for conventional\n", conventional_time);
        printf("%f seconds for Strassen\n", strassen_time);

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