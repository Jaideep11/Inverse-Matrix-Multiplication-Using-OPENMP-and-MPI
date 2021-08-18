#include <iostream>
#include <random>
#include <iomanip>
#include <mpi.h>
#include <stdlib.h>
#include <malloc.h>

using namespace std;

enum MatrixType {
    IDENTITY,
    RANDOM,
    ZERO
};

enum errorCode {
    ZERO_DET
};

double **alloc_2d(int rows, int cols) {
    double *data = (double *)malloc(rows*cols*sizeof(double));
    double **array= (double **)malloc(rows*sizeof(double*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}

double** initMatrix(size_t n, MatrixType type = RANDOM) {
    double **m = alloc_2d(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            switch (type) {
            case RANDOM:
                m[i][j] = (double)(rand() % 1000 + 1);
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

double** initMatrix1(size_t n) {
    double **m = alloc_2d(n, n * 2);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n * 2; j++) {
            m[i][j] = (j < n) ? rand() % 899 + 100 : ((j - n == i) ? 1 : 0);
        }
    }
    return m;
}

void printMatrix(double **m, size_t n) {
    cout << "{\n";
    for (int i = 0; i < n; i++) {
        cout << "    {";
        for (int j = 0; j < n; j++) {
            cout << ((m[i][j] >= 10.0 && m[i][j] < 100.0) ? " " : "") << fixed << setprecision(3) << m[i][j] << (j != n - 1 ? ", " : "");
        }
        cout << "}" << (i != n - 1 ? "," : "") << endl;
    }
    cout << '}';
}

int main(int argc, char **argv) {

    int id, p;
    double a, c;
    
    /*
        id - номер текущего процесса
        p - количесво процессов
        c, a - коэф-ты, которые нужно будет умножать элементы
    */
    double time_start, time_end;
    double **m, **I;
    int request = 0;

    bool test = true;
    int matrix_size = 3;
    
    cout<<"Enter the size of matrix = ";
    // cin>> matrix_size;
    
    if (argc > 1) {
        if ((matrix_size = stoi(argv[1])) > 0) {
            test = false;
        }
    }

    MPI_Status status;

    MPI_Init(&argc, &argv);
    // id данного процесса
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // количество процессов
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    m = initMatrix1(matrix_size);
    I = initMatrix(matrix_size, IDENTITY);
    
    // если это главный процесс, то он инициализирует матрицу
    if (id == 0) {
        cout << matrix_size << ", " << p << ", ";

        if (test) {
            m[0][0] = 1, m[0][1] = 2, m[0][2] = 1, m[0][3] = 1, m[0][4] = 0, m[0][5] = 0,
            m[1][0] = 3, m[1][1] = 5, m[1][2] = 2, m[1][3] = 0, m[1][4] = 1, m[1][5] = 0,
            m[2][0] = 2, m[2][1] = 3, m[2][2] = 3, m[2][3] = 0, m[2][4] = 0, m[2][5] = 1;
        } else {
            m = initMatrix1(matrix_size);
            
        }

        if (matrix_size <= 1000) {
            cout << "\nInput matrix (generated randomly):\n\n";
            // printMatrix(m, matrix_size);
        }

        time_start = MPI_Wtime();
    }
    
    // все процессы получают исходную матрицу и единичную матрицу
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&(m[0][0]), matrix_size*matrix_size*2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double lead[matrix_size * 2];
    int block, start_row, end_row;
    block = matrix_size / (p + 1);
    // cout<<endl<<block<<matrix_size<<"Helllell        aldfldjflksdhalkfhdslkfaldhfladshflasdlfdslflds\n";
    start_row = block * id; 
    end_row = block * (id + 1);
    if (id == p - 1)
        end_row = matrix_size;
    
    for (int i = 0; i < matrix_size; i++) {
    // cout<<"Hello\n"<< matrix_size<<(p + 1)<<"\n";

        if (id == i / block) {
            for (int j = 0; j < matrix_size * 2; j++) {
                lead[j] = m[i][j];
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(lead, matrix_size * 2, MPI_DOUBLE, i / block, MPI_COMM_WORLD);

        for (int j = start_row; j < end_row; j++) {
                // cout<<m[][j];

            if (j == i) {
                double d = m[i][i];
                for (int k = 0; k < matrix_size * 2; k++) 
                    m[j][k] /= d;
                continue;
            }
            double d = m[j][i] / lead[i];
            for (int k = 0; k < matrix_size * 2; k++) {
                m[j][k] -= d * lead[k];
                // cout<<m[j][k];
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);


    cout<<endl<<"Helllellaldfldjflksdhalkfhdslkfaldhfladshflasdlfdslflds\n";
    if (matrix_size <= 10) {
        for (int i = 0; i < matrix_size; i++) {
            if (i / block == id) {
                for (int j = matrix_size; j < matrix_size * 2; j++) {
                    cout << ((m[i][j] >= 0.0 && m[i][j] < 10.0) ? " " : "") << fixed << setprecision(3) << m[i][j] << (j != matrix_size - 1 ? ", " : "");
                }
                cout << endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    printMatrix(I, matrix_size);
    
    if (id == 0) {
        if (matrix_size <= 10 ) {
            cout << "\n\nResult:\n\n";
            printMatrix(I, matrix_size);
        }
        if (id == 0) {
            time_end = MPI_Wtime();
            cout << time_end - time_start << endl;
        }
    } 
    MPI_Finalize();
    
    return 0;
}

