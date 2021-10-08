#include <cmath>
#include <fstream>
#include <iostream>

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;

typedef double MyType;
const MyType eps0 = (typeid(MyType).name()[0] == 'd') ? 1e-12 : 1e-6;
const MyType EPS = 1e-7;
const char *const settingsFileName = "settings.dat";

int readSettings(const char *fileName, size_t &n, string &matrixName, string &vectorName,
                 string &outputName);

void memoryAllocation(MyType **&A, MyType **&L, MyType **&D, MyType **&U, MyType *&b,
                      MyType *&x, int n);

void memoryDeallocation(MyType **&A, MyType **&L, MyType **&D, MyType **&U, MyType *&b,
                        MyType *&x, int n);

int readData(MyType *const *A, MyType *b, const string &matrixName, const string &vectorName);

void printSystem(const MyType *const *A, const MyType *b, int n);

void printVector(const MyType *x, int n);

void printMatrix(const MyType *const *A, int n);

void backward(const MyType *const *A, const MyType *b, MyType *x, int n);

void printResidual(const MyType *const *A, const MyType *b, const MyType *x, int n,
                   MyType normVector(const MyType *const, const int));

MyType computeNorm1(const MyType *x, int n);

MyType computeMatrixNorm1(const MyType *const *A, int n);

MyType computeNormInf(const MyType *x, int n);

MyType computeMatrixNormInf(const MyType *const *A, int n);

int gaussMethod(MyType **A, MyType *b, MyType *x, int n);

int computeInverseMatrix(const MyType *const *A, MyType *const *invA, int n);

void correction(MyType *const *A, MyType *b, int n);

void LDU(const MyType *const *A, MyType *const *L, MyType *const *D, MyType *const *U, int n);

int fixedPointIterMethod(const MyType *const *A, const MyType *b, MyType *x, MyType tau, int n);

int jacobiMethod(const MyType *const *A, const MyType *b, MyType *x, int n);

int seidelMethod(const MyType *const *A, MyType *const *L, MyType *const *D, MyType *const *U,
                 MyType *b, MyType *x, MyType omega, int n);

int seidelThreeDiag(const MyType *a, const MyType *b, const MyType *c, const MyType *d,
                    MyType *x, MyType omega, int n);


int main()  // Точка входа в приложение
{
    setlocale(LC_ALL, "Russian");

    size_t n;  // Размер матрицы
    string matrixFileName, vectorFileName, outputFileName;
    int result = readSettings(settingsFileName, n, matrixFileName,
                              vectorFileName, outputFileName);
    if (result != 0) {
        return result;
    }

    MyType **A,   // Матрица коэффициентов
            **L,  // Нижний треугольник матрицы А
            **D,  // Диагональ матрицы А
            **U,  // Верхний треугольник матрицы A
            *q,   // Вектор правой части
            *x;   // Вектор неизвестных

    memoryAllocation(A, L, D, U, q, x, n);

    result = readData(A, q, matrixFileName, vectorFileName);
    if (result != 0) {
        return result;
    }

    cout << "\nВариант 19, система " << matrixFileName[11] << ", тип "
         << ((typeid(MyType).name()[0] == 'd') ? "double" : "float")
         << ", кубическая норма, точность " << EPS << "\n\n";

    correction(A, q, n);
    printSystem(A, q, n);

    // Решение системы методом простой итерации
    cout << "\nМетод простой итерации:" << '\n';
    MyType tau = 0.01384;  // Оптимальный параметр метода
    cout << "tau = " << tau << " (оптимальный для данной системы)\n";
    fixedPointIterMethod(A, q, x, tau, n);

    // Решение системы методом Якоби
    cout << "\nМетод Якоби:\n";
    jacobiMethod(A, q, x, n);

    // Решение системы методом Зейделя
    cout << "\nМетод Зейделя:\n";
    LDU(A, L, D, U, n);
    seidelMethod(A, L, D, U, q, x, 1, n);

    // Решение системы методом релаксации
    cout << "\nМетод релаксации:\n";
    MyType omega = 1.0;  // Оптимальный параметр метода
    cout << "omega = " << omega << " (оптимальный для данной системы)\n";
    // seidelMethod(A, L, D, U, b, x, omega, n);
    cout << "Совпадает с решением, полученным методом Зейделя\n";

    memoryDeallocation(A, L, D, U, q, x, n);

    // Решение методом релаксации для трехдиагональных матриц
    const int N_VARIANT = 19;
    const int m = 200 + N_VARIANT;

    MyType *a,   // Поддиагональ большой матрицы
            *b,  // Диагональ большой матрицы
            *c,  // Наддиагональ большой матрицы
            *d;  // Правая часть системы

    a = new MyType[m];
    b = new MyType[m];
    c = new MyType[m];
    d = new MyType[m];
    x = new MyType[m];

    cout << "\nМетод релаксации для трехдиагональной матрицы:\n"
            "a[i] = 6,\tb[i] = 8,\tc[i] = 1,\td[i] = i\n";
    for (int i = 0; i < m; ++i) {
        a[i] = 6;
        b[i] = 8;
        c[i] = 1;
        d[i] = i;
    }
    seidelThreeDiag(a, b, c, d, x, 1, m);

    for (int i = 0; i < m; ++i) {
        a[i] = 6;
        b[i] = 8;
        c[i] = 1;
        d[i] = i;
    }

    cout << "\nМетод релаксации для трехдиагональных матриц:\n"
            "a[i] = 6,\tb[i] = 8,\tc[i] = 1,\td[i] = i\n";
    for (int i = 0; i < m; ++i) {
        a[i] = 1;
        b[i] = 8;
        c[i] = 6;
        d[i] = i;
    }
    seidelThreeDiag(a, b, c, d, x, 1, m);

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
    delete[] x;

    return 0;
}


int readSettings(const char *const fileName, size_t &n, string &matrixName,
                 string &vectorName, string &outputName)  // Чтение параметров из файла
{
    ifstream file;
    file.open(fileName);
    if (!file.is_open()) {
        cerr << "Error: file with settings is not open!\n";
        return 1;
    }

    file >> n;
    string str;
    getline(file, str);
    file >> matrixName;
    getline(file, str);
    file >> vectorName;
    getline(file, str);
    file >> outputName;

    file.close();
    return 0;
}


void memoryAllocation(MyType **&A, MyType **&L, MyType **&D, MyType **&U, MyType *&b,
                      MyType *&x, const int n)  // Выделение памяти
{
    A = new MyType *[n];
    L = new MyType *[n];
    D = new MyType *[n];
    U = new MyType *[n];
    for (int i = 0; i < n; ++i) {
        A[i] = new MyType[n];
        L[i] = new MyType[n];
        D[i] = new MyType[n];
        U[i] = new MyType[n];
    }
    b = new MyType[n];
    x = new MyType[n];
}


void memoryDeallocation(MyType **&A, MyType **&L, MyType **&D, MyType **&U, MyType *&b,
                        MyType *&x, const int n)  // Освобождение памяти
{
    for (int i = 0; i < n; ++i) {
        delete[] A[i];
        delete[] L[i];
        delete[] D[i];
        delete[] U[i];
    }
    delete[] A;
    delete[] L;
    delete[] D;
    delete[] U;
    delete[] b;
    delete[] x;
}


int readData(MyType *const *const A, MyType *const b, const string &matrixName,
             const string &vectorName)  // Чтение данных из файлов
{
    ifstream file;
    file.open(matrixName);
    if (!file.is_open()) {
        cerr << "Error: file with matrix is not open!\n";
        return 2;
    }

    int n;
    file >> n;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file >> A[i][j];
        }
    }

    file.close();
    file.open(vectorName);
    if (!file.is_open()) {
        cerr << "Error: file with vector is not open!\n";
        return 3;
    }

    file >> n;
    for (int i = 0; i < n; ++i) {
        file >> b[i];
    }

    file.close();
    return 0;
}


void printSystem(const MyType *const *const A, const MyType *const b, const int n)  // Вывод системы уравнений
{
    cout << "Система уравнений:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (j > 0) {
                cout << ((A[i][j] < 0) ? " - " : " + ");
            } else {
                cout << ((A[i][j] < 0) ? "-" : " ");
            }
            cout << std::setiosflags(std::ios::left)
                 << std::setiosflags(std::ios::fixed)
                 << std::setprecision(3)
                 << std::setw(8) << fabs(A[i][j])
                 << std::setprecision(6)
                 << std::resetiosflags(std::ios::fixed)
                 << std::resetiosflags(std::ios::left)
                 << " * x" << j + 1 << "\t";
        }
        cout << " = " << b[i] << '\n';
    }
}


void printVector(const MyType *const x, const int n)  // Вывод вектора решения
{
    cout << "{ " << x[0];
    for (int i = 1; i < n; ++i) {
        cout << ", " << x[i];
    }
    cout << " }\n";
}


void printMatrix(const MyType *const *const A, const int n)  // Вывод матрицы
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << std::setw(15) << A[i][j];
        }
        cout << '\n';
    }
}


void backward(const MyType *const *const A, const MyType *const b, MyType *const x, const int n)  // Обратный ход
{
    for (int i = 1; i < n; ++i) {
        x[i] = 0.0;
    }

    MyType t;
    for (int i = n - 1; i >= 0; --i) {
        t = b[i];
        for (int j = i + 1; j < n; ++j) {
            t -= A[i][j] * x[j];
        }

        x[i] = t / A[i][i];
        if (fabs(x[i]) < eps0) {
            x[i] = 0.0;
        }
    }
}


void printResidual(const MyType *const *const A, const MyType *const b, const MyType *const x, const int n,
                   MyType normVector(const MyType *const, const int))  // Подсчет невязки
{
    MyType *residual = new MyType[n];

    for (int i = 0; i < n; ++i) {
        MyType bi = 0;
        for (int j = 0; j < n; ++j) {
            bi += A[i][j] * x[j];
        }
        residual[i] = fabs(bi - b[i]);
    }

    cout << "Норма вектора невязки: " << normVector(residual, n) << '\n';

    delete[] residual;
}


MyType computeNorm1(const MyType *const x, const int n)  // Октаэдрическая норма вектора
{
    MyType sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += fabs(x[i]);
    }
    return sum;
}


MyType computeMatrixNorm1(const MyType *const *const A, const int n)  // Октаэдрическая норма матрицы
{
    MyType max = 0.0;
    for (int j = 0; j < n; ++j) {
        MyType sum = 0.0;
        for (int i = 0; i < n; ++i) {
            sum += fabs(A[i][j]);
        }
        if (sum > max) {
            max = sum;
        }
    }
    return max;
}


MyType computeNormInf(const MyType *const x, const int n)  // Кубическая норма вектора
{
    MyType nax = 0.0;
    for (int i = 0; i < n; ++i) {
        if (fabs(x[i]) > nax) {
            nax = fabs(x[i]);
        }
    }
    return nax;
}


MyType computeMatrixNormInf(const MyType *const *const A, const int n)  // Кубическая норма матрицы
{
    MyType nax = 0.0;
    for (int i = 0; i < n; ++i) {
        MyType sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += fabs(A[i][j]);
        }
        if (sum > nax) {
            nax = sum;
        }
    }
    return nax;
}


int gaussMethod(MyType **const A, MyType *const b, MyType *const x, const int n)  // Метода Гаусса
{
    // Далее: k - диагональ, i - строка, j - столбец
    for (int k = 0; k < n; ++k) {
        // Частичный поиск ведущего элемента по столбцу
        int iMax = k;
        for (int i = k; i < n; ++i) {
            if (fabs(A[i][k]) > fabs(A[iMax][k])) {
                iMax = k;
            }
        }
        // Если на диагонали нулевой элемент
        if (fabs(A[iMax][k]) < eps0) {
            return 1;
        }
        // Меняем строки
        if (iMax != k) {
            MyType *temp = A[k];
            A[k] = A[iMax];
            A[iMax] = temp;

            MyType t = b[k];
            b[k] = b[iMax];
            b[iMax] = t;
        }
        // Прямой ход
        for (int i = k + 1; i < n; ++i) {
            MyType t = A[i][k] / A[k][k];
            for (int j = k; j < n; ++j) {
                A[i][j] -= A[k][j] * t;
            }
            b[i] -= b[k] * t;
        }
    }

    backward(A, b, x, n);

    return 0;
}


int computeInverseMatrix(const MyType *const *const A, MyType *const *const invA,
                         const int n)  // Вычисление обратой матрицы
{
    MyType **A_copy,
            *x_j,
            *e_j;

    A_copy = new MyType *[n];
    x_j = new MyType[n];
    e_j = new MyType[n];

    for (int i = 0; i < n; ++i) {
        A_copy[i] = new MyType[n];
        for (int j = 0; j < n; ++j) {
            A_copy[i][j] = A[i][j];
        }
    }

    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            e_j[i] = (i == j) ? 1 : 0;
        }
        int result = gaussMethod(A_copy, e_j, x_j, n);
        if (result != 0) {
            return result;
        }

        for (int i = 0; i < n; ++i) {
            invA[i][j] = x_j[i];
        }
    }

    for (int i = 0; i < n; ++i) {
        delete[] A_copy[i];
    }
    delete[] A_copy;
    delete[] x_j;
    delete[] e_j;

    return 0;
}


void correction(MyType *const *const A, MyType *const b, const int n)  // Корректировка системы
{
    MyType *sum = new MyType[n];
    for (int i = 0; i < n; ++i) {
        sum[i] = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                sum[i] += fabs(A[i][j]);
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        if (A[i][i] < 0.0 && fabs(A[i][i]) > sum[i]) {
            for (int j = 0; j < n; ++j) {
                A[i][j] = -A[i][j];
            }
            b[i] = -b[i];
        }
    }
    delete[] sum;
}


void LDU(const MyType *const *const A, MyType *const *const L, MyType *const *const D,
         MyType *const *const U, const int n)  // LDU-разложение матрицы
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                D[i][j] = A[i][j];
                L[i][j] = 0;
                U[i][j] = 0;
            } else if (i > j) {
                D[i][j] = 0;
                L[i][j] = A[i][j];
                U[i][j] = 0;
            } else {
                D[i][j] = 0;
                L[i][j] = 0;
                U[i][j] = A[i][j];
            }
        }
    }
}


int fixedPointIterMethod(const MyType *const *const A, const MyType *const b, MyType *const x,
                         const MyType tau, const int n)  // Метод простой итерации
{
    MyType **C,
            *y,
            *x_k,
            *diff;

    C = new MyType *[n];
    for (int i = 0; i < n; ++i) {
        C[i] = new MyType[n];
    }
    y = new MyType[n];
    x_k = new MyType[n];
    diff = new MyType[n];

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            MyType e_ij = (i == j) ? 1 : 0;
            C[i][j] = e_ij - tau * A[i][j];
        }
        y[i] = tau * b[i];
        x_k[i] = y[i];
    }

    cout << "Вектор y: ";
    printVector(y, n);
    cout << "Матрица C:" << '\n';
    printMatrix(C, n);
    MyType cNorm = computeMatrixNormInf(C, n);
    cout << "Норма матрицы C: " << cNorm << '\n';

    int iter = 0;
    do {
        if (iter++ != 0) {
            for (int i = 0; i < n; ++i) {
                x_k[i] = x[i];
            }
        }

        for (int i = 0; i < n; ++i) {
            MyType t = 0.0;
            for (int j = 0; j < n; ++j) {
                t += C[i][j] * x_k[j];
            }
            x[i] = t + y[i];
        }

        for (int i = 0; i < n; ++i) {
            diff[i] = x[i] - x_k[i];
        }
    } while (computeNormInf(diff, n) > (EPS * (1 - cNorm)) / cNorm);

    cout << "Число итераций: " << iter << '\n';
    cout << "Вектор x: ";
    printVector(x, n);
    // todo: норма ошибки
    printResidual(A, b, x, n, computeNormInf);

    for (int i = 0; i < n; ++i) {
        delete[] C[i];
    }
    delete[] C;
    delete[] y;
    delete[] x_k;
    delete[] diff;

    return 0;
}


int jacobiMethod(const MyType *const *const A, const MyType *const b, MyType *const x,
                 const int n)  // Метод Якоби
{
    MyType **C,
            *y,
            *x_k,
            *diff;

    C = new MyType *[n];
    for (int i = 0; i < n; ++i) {
        C[i] = new MyType[n];
    }
    y = new MyType[n];
    x_k = new MyType[n];
    diff = new MyType[n];

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i][j] = (i == j) ? 0 : -A[i][j] / A[i][i];
        }
        y[i] = b[i] / A[i][i];
        x_k[i] = y[i];
    }

    cout << "Вектор y: ";
    printVector(y, n);
    cout << "Матрица C: " << '\n';
    printMatrix(C, n);
    MyType cNorm = computeMatrixNormInf(C, n);
    cout << "Норма матрицы C: " << cNorm << '\n';

    int iter = 0;
    do {
        if (iter++ != 0) {
            for (int i = 0; i < n; ++i) {
                x_k[i] = x[i];
            }
        }

        for (int i = 0; i < n; ++i) {
            MyType t = 0.0;
            for (int j = 0; j < n; ++j) {
                t += C[i][j] * x_k[j];
            }
            x[i] = t + y[i];
        }

        for (int i = 0; i < n; ++i) {
            diff[i] = x[i] - x_k[i];
        }
    } while (computeNormInf(diff, n) > (EPS * (1 - cNorm)) / cNorm);

    cout << "Число итераций: " << iter << '\n';
    cout << "Вектор x: ";
    printVector(x, n);
    // todo: норма ошибки
    printResidual(A, b, x, n, computeNormInf);

    for (int i = 0; i < n; ++i) {
        delete[] C[i];
    }
    delete[] C;
    delete[] y;
    delete[] x_k;
    delete[] diff;

    return 0;
}


int seidelMethod(const MyType *const *const A, MyType *const *const L, MyType *const *const D,
                 MyType *const *const U, MyType *const b, MyType *const x,
                 const MyType omega, const int n)  // Метод релаксации
{
    MyType **B,
            **invB,
            **C,
            **C_L,
            **C_D,
            **C_U,
            *y,
            *x_k,
            *diff;

    B = new MyType *[n];
    invB = new MyType *[n];
    C = new MyType *[n];
    C_L = new MyType *[n];
    C_D = new MyType *[n];
    C_U = new MyType *[n];
    for (int i = 0; i < n; ++i) {
        B[i] = new MyType[n];
        invB[i] = new MyType[n];
        C[i] = new MyType[n];
        C_L[i] = new MyType[n];
        C_D[i] = new MyType[n];
        C_U[i] = new MyType[n];
    }
    y = new MyType[n];
    x_k = new MyType[n];
    diff = new MyType[n];

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            B[i][j] = D[i][j] + omega * L[i][j];
        }
    }

    int result = computeInverseMatrix(B, invB, n);
    if (result != 0) {
        cout << "Матрица B вырождена\n";
        return 1;
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            MyType sum = 0.0;
            for (int k = 0; k < n; ++k) {
                sum += invB[i][k] * A[k][j];
            }
            MyType e_ij = (i == j) ? 1 : 0;
            C[i][j] = e_ij - omega * sum;
        }
    }

    cout << "Матрица C: " << '\n';
    printMatrix(C, n);
    MyType cNorm = computeMatrixNormInf(C, n);
    cout << "Норма матрицы С: " << cNorm << '\n';

    LDU(C, C_L, C_D, C_U, n);
    MyType clNorm = computeMatrixNormInf(C_L, n);
    cout << "||C_L||_inf = " << clNorm << '\n';
    MyType cdNorm = computeMatrixNormInf(C_D, n);
    cout << "||C_D||_inf = " << cdNorm << '\n';
    MyType cuNorm = computeMatrixNormInf(C_U, n);
    cout << "||C_U||_inf = " << cuNorm << '\n';
    cout << "||C_L||_inf + ||C_U||_inf = " << clNorm + cuNorm << '\n';

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i][j] = (i == j) ? 0 : -A[i][j] / A[i][i];
        }
        y[i] = omega * b[i] / A[i][i];
        x_k[i] = y[i];
    }

    int iter = 0;
    do {
        if (iter++ != 0) {
            for (int i = 0; i < n; ++i) {
                x_k[i] = x[i];
            }
        }

        for (int i = 0; i < n; ++i) {
            MyType sum = 0.0;
            for (int j = i + 1; j < n; ++j) {
                sum += C[i][j] * x_k[j];
            }
            x[i] = (1 - omega) * x_k[i] + omega * sum + y[i];
        }

        for (int i = 0; i < n; ++i) {
            MyType sum = 0.0;
            for (int j = 0; j < i; ++j) {
                sum += C[i][j] * x[j];
            }
            x[i] += omega * sum;
        }

        for (int i = 0; i < n; ++i) {
            diff[i] = x[i] - x_k[i];
        }
    } while (computeNormInf(diff, n) > (EPS * (1 - cNorm)) / cuNorm);

    cout << "Число итераций: " << iter << '\n';
    cout << "Вектор x: ";
    printVector(x, n);
    // todo: норма ошибки
    printResidual(A, b, x, n, computeNormInf);

    for (int i = 0; i < n; ++i) {
        delete[] B[i];
        delete[] invB[i];
        delete[] C[i];
        delete[] C_L[i];
        delete[] C_D[i];
        delete[] C_U[i];
    }
    delete[] B;
    delete[] invB;
    delete[] C;
    delete[] C_L;
    delete[] C_D;
    delete[] C_U;
    delete[] y;
    delete[] x_k;
    delete[] diff;

    return 0;
}


int seidelThreeDiag(const MyType *const a, const MyType *const b, const MyType *const c,
                    const MyType *const d, MyType *const x, const MyType omega,
                    const int n)  // Метод релаксации для трехдиагональных матриц
{
    MyType *y,
            *x_k,
            *diff;

    y = new MyType[n];
    x_k = new MyType[n];
    diff = new MyType[n];

    for (int i = 0; i < n; ++i) {
        y[i] = omega * d[i] / b[i];
        x_k[i] = y[i];
    }

    int iter = 0;
    do {
        if (iter++ != 0) {
            for (int i = 0; i < n; ++i) {
                x_k[i] = x[i];
            }
        }

        for (int i = 0; i < n; ++i) {
            MyType sum = 0.0;
            for (int j = i + 1; j < n; ++j) {
                sum += a[j] / b[i] * x_k[j];
            }
            x[i] = (1 - omega) * x_k[i] - omega * sum + y[i];
        }

        for (int i = 0; i < n; ++i) {
            MyType sum = 0.0;
            for (int j = 0; j < i; ++j) {
                sum += a[j] / b[i] * x[j];
            }
            x[i] -= omega * sum;
        }

        for (int i = 0; i < n; ++i) {
            diff[i] = x[i] - x_k[i];
        }
    } while (computeNormInf(diff, n) > EPS);

    cout << "Число итераций: " << iter << '\n';
    cout << "Вектор x: ";
    printVector(x, n);

    delete[] y;
    delete[] x_k;
    delete[] diff;

    return 0;
}