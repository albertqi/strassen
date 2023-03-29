#include <stdio.h>
#include <unordered_set>
#include <string>
#include <fstream>
#include <random>
#include <chrono>
#include <thread>
#include <mutex>
#include <limits>

#define N0 65
#define NUM_TRIALS 10

long zero = 0;
thread_local std::mt19937 generator;
std::mutex m;
std::unordered_set<long *> pointers;

/**
 * Matrix with support for both conventional and Strassen multiplication.
 */
class Matrix
{
private:
    const int original_dim, dim, row_start, col_start;
    long *arr;

    /**
     * Constructs a new `Matrix` object.
     *
     * @param original_dim Original dimension when first created.
     * @param dim Current dimension of matrix.
     * @param row_start Starting row of matrix.
     * @param col_start Starting column of matrix.
     * @param arr Pointer to array.
     */
    Matrix(const int &original_dim, const int &dim, const int &row_start, const int &col_start, long *arr)
        : original_dim(original_dim), dim(dim), row_start(row_start), col_start(col_start), arr(arr)
    {
        m.lock();
        pointers.insert(arr);
        m.unlock();
    }

    /**
     * Finds the `q`-th quarter.
     *
     * @param q The quarter to be returned.
     * @return Quarter `q` of matrix.
     */
    const Matrix quarter(const int &q) const
    {
        const int new_dim = (dim + 1) / 2;
        const int new_row = q >= 2 ? row_start + new_dim : row_start,
                  new_col = q % 2 == 1 ? col_start + new_dim : col_start;
        return Matrix(original_dim, new_dim, new_row, new_col, arr);
    }

    /**
     * Shifts matrix `X` by `delta` * `Y`.
     *
     * @param X First matrix.
     * @param Y Second matrix.
     * @param delta Shift constant.
     * @return Sum `X` + `delta` * `Y`.
     */
    static const Matrix shift(const Matrix &X, const Matrix &Y, const int &delta)
    {
        // Handle mismatched dimensions.
        if (X.original_dim != Y.original_dim || X.dim != Y.dim)
        {
            throw std::invalid_argument("Mismatched dimensions");
        }

        // Create new matrix `Z`.
        const int dim = X.dim;
        const Matrix Z(dim, new long[dim * dim]);

        // Update entries of `Z`.
        for (int i = 0; i < dim; ++i)
        {
            for (int j = 0; j < dim; ++j)
            {
                Z(i, j) = X(i, j) + delta * Y(i, j);
            }
        }

        return Z;
    }

    /**
     * Adds matrices `X` and `Y`.
     *
     * @param X First matrix.
     * @param Y Second matrix.
     * @return Sum `X` + `Y`.
     */
    static const Matrix add(const Matrix &X, const Matrix &Y)
    {
        return shift(X, Y, 1);
    }

    /**
     * Subtracts matrices `X` and `Y`.
     *
     * @param X First matrix.
     * @param Y Second matrix.
     * @return Difference `X` - `Y`.
     */
    static const Matrix subtract(const Matrix &X, const Matrix &Y)
    {
        return shift(X, Y, -1);
    }

    /**
     * Indexes into current matrix at position `row` and `col`,
     * or `zero` if position is invalid.
     *
     * @param row Row number.
     * @param col Column number.
     * @return Reference to position `row` and `col` in matrix.
     */
    long &operator()(const int &row, const int &col) const
    {
        if (row_start + row >= original_dim || col_start + col >= original_dim)
        {
            return zero;
        }

        return arr[(row_start + row) * original_dim + col_start + col];
    }

public:
    /**
     * Constructs a new `Matrix` object.
     *
     * @param dim Current dimension of matrix.
     * @param arr Pointer to array.
     */
    Matrix(const int &dim, long *arr)
        : original_dim(dim), dim(dim), row_start(0), col_start(0), arr(arr)
    {
        m.lock();
        pointers.insert(arr);
        m.unlock();
    }

    /**
     * Conventional matrix multiplication.
     *
     * @param X First matrix.
     * @param Y Second matrix.
     * @return Product `X` * `Y`.
     */
    static const Matrix multiply_conventional(const Matrix &X, const Matrix &Y)
    {
        // Handle mismatched dimensions.
        if (X.dim != Y.dim)
        {
            throw std::invalid_argument("Mismatched dimensions");
        }

        // Create new matrix `Z`.
        const int dim = X.dim;
        const Matrix Z(dim, new long[dim * dim]);

        // Update entries of `Z`.
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

    /**
     * Strassen multiplication with an input `threshold`.
     *
     * @param X First matrix.
     * @param Y Second matrix.
     * @param threshold Dimension at which we switch to conventional multiplication.
     * @return Product `X` * `Y`.
     */
    static const Matrix multiply_strassen_threshold(const Matrix &X,
                                                    const Matrix &Y,
                                                    const int &threshold)
    {
        // Handle mismatched dimensions.
        if (X.dim != Y.dim)
        {
            throw std::invalid_argument("Mismatched dimensions");
        }

        // Store dimension `dim`.
        const int dim = X.dim;

        // Use conventional multiplication if `dim <= threshold`.
        if (dim <= threshold)
        {
            return multiply_conventional(X, Y);
        }

        // Partition `X` and `Y` into quarters.
        const Matrix A = X.quarter(0),
                     B = X.quarter(1),
                     C = X.quarter(2),
                     D = X.quarter(3),
                     E = Y.quarter(0),
                     F = Y.quarter(1),
                     G = Y.quarter(2),
                     H = Y.quarter(3);

        // Calculate products `P1` through `P7`.
        const Matrix P1 = multiply_strassen_threshold(A, subtract(F, H), threshold),
                     P2 = multiply_strassen_threshold(add(A, B), H, threshold),
                     P3 = multiply_strassen_threshold(add(C, D), E, threshold),
                     P4 = multiply_strassen_threshold(D, subtract(G, E), threshold),
                     P5 = multiply_strassen_threshold(add(A, D), add(E, H), threshold),
                     P6 = multiply_strassen_threshold(subtract(B, D), add(G, H), threshold),
                     P7 = multiply_strassen_threshold(subtract(C, A), add(E, F), threshold);

        // Calculate new quarters for product.
        const Matrix Q1 = add(P4, add(P5, subtract(P6, P2))),
                     Q2 = add(P1, P2),
                     Q3 = add(P3, P4),
                     Q4 = add(P1, add(P5, subtract(P7, P3)));

        // Create new matrix `Z`.
        const Matrix Z(dim, new long[dim * dim]);
        const bool should_trim = Q1.dim * 2 != dim;

        // Update entries of `Z`.
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

    /**
     * Strassen multiplication with threshold `N0`.
     *
     * @param X First matrix.
     * @param Y Second matrix.
     * @return Product `X` * `Y`.
     */
    static const Matrix multiply_strassen(const Matrix &X, const Matrix &Y)
    {
        return multiply_strassen_threshold(X, Y, N0);
    }

    /**
     * Prints diagonal entries.
     */
    void print_diagonal() const
    {
        for (int i = 0; i < dim; ++i)
        {
            printf("%ld\n", (*this)(i, i));
        }
    }

    /**
     * Finds sum of diagonal entries (i.e., the trace).
     *
     * @return Sum of diagonal entries.
     */
    const long sum_diagonal() const
    {
        long x = 0;
        for (int i = 0; i < dim; ++i)
        {
            x += (*this)(i, i);
        }
        return x;
    }
};

/**
 * Finds best threshold for Strassen multiplication.
 *
 * @param A First matrix.
 * @param B Second matrix.
 * @param a Pointer to array of `A`.
 * @param b Pointer to array of `B`.
 * @param dim Dimension of `A` and `B`.
 */
void find_best_threshold(const Matrix &A, const Matrix &B, long *a, long *b, const int &dim)
{
    // Keep track of `min_time` and `best_threshold`.
    double min_time = std::numeric_limits<double>::max();
    int best_threshold;

    // Calculate product `C` for various values of `threshold`.
    for (int threshold = 10; threshold <= dim; ++threshold)
    {
        // Initialize `total_runtime`.
        double total_runtime = 0;

        for (int i = 0; i < NUM_TRIALS; ++i)
        {
            // Time `multiply_strassen_threshold`.
            const auto start = std::chrono::high_resolution_clock::now();
            const Matrix C = Matrix::multiply_strassen_threshold(A, B, threshold);
            const auto end = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double> diff = end - start;

            // Update `total_runtime`.
            total_runtime += diff.count();

            // Free dynamically allocated memory.
            for (const auto &p : pointers)
            {
                if (p != a && p != b)
                {
                    delete[] p;
                }
            }

            pointers = std::unordered_set<long *>({a, b});
        }

        // Calculate average runtime.
        const double runtime = total_runtime / NUM_TRIALS;

        // Print results.
        printf("%f seconds with threshold of %d\n", runtime, threshold);

        // Update `min_time` and `best_threshold`.
        if (runtime < min_time)
        {
            min_time = runtime;
            best_threshold = threshold;
        }
    }

    // Print `best_threshold`.
    printf("%d\n", best_threshold);

    // Free dynamically allocated memory.
    for (const auto &p : pointers)
    {
        delete[] p;
    }
}

/**
 * Counts number of triangles in a random graph.
 *
 * @param flag Probability of including an edge.
 */
void count_triangles(const double &flag)
{
    // Initialize `total_triangles`.
    long total_triangles = 0;

    std::vector<std::thread> threads;
    for (int i = 0; i < NUM_TRIALS; ++i)
    {
        threads.push_back(std::thread([&flag, &total_triangles]()
        {
            // Set up random generator.
            std::random_device random_dev;
            generator.seed(random_dev());
            std::uniform_real_distribution<double> unif(0.0, 1.0);

            // Initialize adjacency matrix.
            long *x = new long[1024 * 1024];
            for (int i = 0; i < 1024; ++i)
            {
                for (int j = 0; j < i; ++j)
                {
                    // Add edge with probability `flag`.
                    const bool include_edge = unif(generator) <= flag;
                    x[i * 1024 + j] = include_edge;
                    x[j * 1024 + i] = include_edge;
                }
            }

            // Create matrix `X`.
            const Matrix X(1024, x);

            // Calculate product `Y`.
            const Matrix Y = Matrix::multiply_strassen(X, Matrix::multiply_strassen(X, X));

            // Update `total_triangles`.
            m.lock();
            total_triangles += Y.sum_diagonal() / 6;
            m.unlock();
        }));
    }
    for (std::thread &t : threads)
    {
        t.join();
    }

    // Print results.
    printf("%f\n", (double)total_triangles / NUM_TRIALS);

    // Free dynamically allocated memory.
    for (const auto &p : pointers)
    {
        delete[] p;
    }
}

int main(int argc, char *argv[])
{
    // Check for valid command line arguments.
    if (argc != 4)
    {
        printf("Usage: ./strassen [flag] [dimension] [input_file]\n");
        return 1;
    }

    // Parse command line arguments.
    const double flag = atof(argv[1]);
    const int dim = atoi(argv[2]);
    const std::string input_file = argv[3];

    // Count triangles in random graph.
    if (flag > 0 && flag < 1)
    {
        count_triangles(flag);
        return 0;
    }

    // Open input file.
    std::ifstream file(input_file);

    // Parse input file.
    long *a = new long[dim * dim], *b = new long[dim * dim];
    for (int i = 0; i < dim * dim; ++i)
    {
        file >> a[i];
    }
    for (int i = 0; i < dim * dim; ++i)
    {
        file >> b[i];
    }

    // Create matrices `A` and `B`.
    const Matrix A(dim, a), B(dim, b);

    // Find best threshold for Strassen multiplication.
    if (flag == 1)
    {
        find_best_threshold(A, B, a, b, dim);
        return 0;
    }

    // Calculate product `C`.
    const Matrix C = Matrix::multiply_strassen(A, B);

    // Print results.
    C.print_diagonal();

    // Free dynamically allocated memory.
    for (const auto &p : pointers)
    {
        delete[] p;
    }

    return 0;
}