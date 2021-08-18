#include <iostream>
#include <random>
#include <iomanip>
#include <omp.h>
#include <fstream>

using namespace std;

enum MatrixType {
    IDENTITY,
    RANDOM,
    ZERO
};

enum errorCode {
    ZERO_DET
};

double** initMatrix(size_t n, MatrixType type = RANDOM) {
    double **m = new double*[n];
    for (int i = 0; i < n; i++) {
        m[i] = new double[n];
        for (int j = 0; j < n; j++) {
            switch (type) {
            case RANDOM:
                m[i][j] = (double)(rand() % 100);
                break;
            case IDENTITY:
                m[i][j] = (i == j) ? 1 : 0;
                break;
            default:
                m[i][j] = 0;
                break;
            }
        }
    }
    return m;
}

double **copyMatrix(double **m, size_t n) {
    double **m1 = new double*[n];
    for (int i = 0; i < n; i++) {
        m1[i] = new double[n];
        for (int j = 0; j < n; j++) {
            m1[i][j] = m[i][j];
        }
    }
    return m1;
}

void printMatrix(double **m, size_t n) {
    cout << "{\n";
    for (int i = 0; i < n; i++) {
        cout << "    {";
        for (int j = 0; j < n; j++) {
            cout << ((m[i][j] >= 0.0 && m[i][j] < 10.0) ? " " : "") << fixed << setprecision(3) << m[i][j] << (j != n - 1 ? ", " : "");
        }
        cout << "}" << (i != n - 1 ? "," : "") << endl;
    }
    cout << '}';
}



void multRow(double **m, int row, double k, size_t n) {
    for (int i = 0; i < n; i++) {
        m[row][i] *= k;
    }
}

void addRow(double **m, int rowPut, int rowGet, size_t n, double k = 1.0) {
    for (int i = 0; i < n; i++) {
        m[rowPut][i] += k * m[rowGet][i];
    }
}

void swapRows(double **m, int row1, int row2, size_t n) {
    double *tmp = m[row1];
    m[row1] = m[row2];
    m[row2] = tmp;
}

double **inverseMatrix(double **m_orig, size_t n) {
    double **m = copyMatrix(m_orig, n);
    double **I = initMatrix(n , IDENTITY);
    int i, j, k;
    for (j = 0; j < n; j++) {
        // swap lines if a_jj == 0
        if (m[j][j] == 0) {
            bool flag = true;
            for (k = j; k < n; k++) {
                if (m[k][j] != 0) {
                    swapRows(m, j, k, n);
                    flag = false;
                    break;
                }
            }
            if (flag) {
                throw ZERO_DET;
            }
        }

        // делим строку на a_jj
        double div = m[j][j];
        for (int i = 0; i < n; i++) {
            m[j][i] /= div;
            I[j][i] /= div;
        }
        // из всех строк, кроме этой, вычитаем эту с коэф-ом
        for (int i = 0; i < n; i++) {
            if (i == j)
                continue;
                double c = m[i][j];
            for (int k = 0; k < n; k++) {
                m[i][k] -= c * m[j][k];
                I[i][k] -= c * I[j][k];
            }
        }

        // old
        // for (i = 0; i < n; i++) {
        //     if (i == j)
        //         continue;
        //     addRow(I, i, j, n, -m[i][j]/m[j][j]);
        //     addRow(m, i, j, n, -m[i][j]/m[j][j]);
        // }
        // multRow(I, j, 1.0/m[j][j], n);
        // multRow(m, j, 1.0/m[j][j], n);
    }
    return I;
}

double **inverseMatrixOpenMP(double **m, size_t n, int threadsNum = 8) {
    double **I = initMatrix(n , IDENTITY);

    omp_set_num_threads(threadsNum);
    int i, j, k;

#pragma omp parallel for shared(m, I) private(i, j, k)
    for (j = 0; j < n; j++) {
        // swap lines if a_jj == 0
        if (m[j][j] == 0) {
            bool flag = true;
            for (k = j; k < n; k++) {
                if (m[k][j] != 0) {
                    swapRows(m, j, k, n);
                    flag = false;
                    break;
                }
            }
            if (flag) {
                throw ZERO_DET;
            }
        }
        // делим строку на a_jj
        double div = m[j][j];
        for (int i = 0; i < n; i++) {
            m[j][i] /= div;
            I[j][i] /= div;
        }
        // вот здесь надо будет распараллелить: из всех строк, кроме этой, вычитаем эту с коэф-ом
        for (int i = 0; i < n; i++) {
            if (i == j)
                continue;
                double c = m[i][j];
            
            for (int k = 0; k < n; k++) {
                m[i][k] -= c * m[j][k];
                I[i][k] -= c * I[j][k];
            }
        }
    }
    return I;
}

int main(int argc, char **argv) {

    bool test = false;
    int matrix_size = 3;
    int threads = 8;
    
    cout<<"Enter the size of matrix = ";
    cin>> matrix_size;

    if (argc > 1) {
        if ((matrix_size = stoi(argv[1])) > 0) {
            test = false;
        }
        if (argc > 2) {
            if (stoi(argv[2]) > 0)
                threads = stoi(argv[2]);
        }
    }

    double **T = initMatrix(matrix_size, ZERO);

    if (test) {
        T[0][0] = 1, T[0][1] = 2, T[0][2] = 1,
        T[1][0] = 3, T[1][1] = 5, T[1][2] = 2,
        T[2][0] = 2, T[2][1] = 3, T[2][2] = 3;
    } else {
        T = initMatrix(matrix_size, RANDOM);
    }

    if (matrix_size <= 1000) {
        cout << "\nInput matrix (generated randomly):\n\n";
        printMatrix(T, matrix_size);
    }

    double timerOpenMp;

    double **N;

    // обычный алгоритм

    try {
        timerOpenMp = omp_get_wtime();
        N = inverseMatrix(T, matrix_size);
        timerOpenMp = omp_get_wtime() - timerOpenMp;
    }
    catch(errorCode err) {
        if (err == ZERO_DET) {
            cout << "Inverse matrix doesn't exist\n";
            return 0;
        }
    }
    if (matrix_size <= 10000) {
        cout << "\n\nResult:\n\n";
        printMatrix(N, matrix_size);
    }
    
    cout << "\nRegular algorithm execution time: " << timerOpenMp << endl;

    // OpenMP

    cout << "Matrix size: " << matrix_size << ", threads: " << threads << endl;

    try {
        timerOpenMp = omp_get_wtime();
        N = inverseMatrixOpenMP(T, matrix_size, threads);
        timerOpenMp = omp_get_wtime() - timerOpenMp;
    }
    catch(errorCode err) {
        if (err == ZERO_DET) {
            cout << "Inverse matrix doesn't exist\n";
            return 0;
        }
    }
    if (matrix_size <= 10000) {
        cout << "\n\nResult:\n\n";
        printMatrix(N, matrix_size);
    }
    
    cout << "\nOpenMP execution time: " << timerOpenMp << endl;


    return 0;
}
