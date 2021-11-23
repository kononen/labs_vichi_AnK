#include "iostream"
#include <fstream>
#include <typeinfo>
#include <cmath>

using namespace std;

const double EPSILON = 1e-8;   //��������� ��������� � 0
const double EPS_eigen = 1e-4;
const double PI = 3.14159265;

typedef double mytipe; //��� ������, �������������� �� ���� ���������

#include "Matr.h"

mytipe* eigenvalues(const int N, mytipe** A);
void Hessenberg_matrix(const unsigned int DIM, mytipe** A);
mytipe** eigenvectors(const int N, mytipe** A, mytipe* lambda);
void Raley(int N, mytipe* lambda, mytipe** e,mytipe** A);
bool QThis_Eigenvalues_A(int N, mytipe* lambda, mytipe** A);
bool QThis_Eigenvectors_A(int N, mytipe* lambda, mytipe** e, mytipe** A);

int main() {
	setlocale(LC_ALL, "Russian");

	ifstream fin;    // ������� ���������� ��� ���������� �� �����
	fin.open("1.txt", ios_base::in | ios_base::app | ios_base::binary);  //������� ����
	unsigned int a;  // ������ ������������
	fin >> a;  //������� �� ����� ������ ������������
	cout << "Dim=" << a << endl;
	const unsigned int DIM = a;  //������ ������������� ������� �����������
	mytipe **A;                      // |������� 
	A = new mytipe *[DIM];          // |������������ 
	mytipe **e;                      // |������� 
	e = new mytipe *[DIM];          // |������������
	for (int i = 0; i < DIM; i++) {   // |������ 
		A[i] = new mytipe[DIM];     // |������� �
		e[i] = new mytipe[DIM];
	}

	
   for (int i = 0; i < DIM; i++)  //��������� �� ����� ������� � �������
	{
		for (int j = 0; j < DIM; j++)
		{
			fin >> A[i][j];  //�������
		}
	}

	fin.close();  //��������� ����
	fin.clear();  //���������� ���� �� �������� ����������

	cout << "\nA=";
	vivod(DIM, A);

	mytipe* lambda = new mytipe[DIM];

	cout << "\n\t\tQR-��������\n\n";
	lambda = eigenvalues(DIM, A);
	cout << endl;
	cout << "����������� �������� ������� A:"; vivod(DIM, lambda);
	//QThis_Eigenvalues_A(DIM,lambda,A);

	cout << "\n\t\t����� �������� ��������\n\n";
	cout << "����������� ������� ������� A:\n\n";
	copy(DIM, e, eigenvectors(DIM, A, lambda));
	cout << endl;
	for (int i = 0; i < DIM; i++) { cout << "e" << i + 1 << " = "; vivod(DIM, e[i]); }
	//QThis_Eigenvectors_A(DIM, lambda, e, A);

	cout << "\n\t\t� ������� ����������� �����\n\n";
	Raley(DIM, lambda, e, A);
	cout << "\n����������� �������� ������� A:";
	vivod(DIM, lambda);
	cout << "\n����������� ������� ������� A:\n";
	for (int i = 0; i < DIM; i++) { cout << "e" << i + 1 << " = "; vivod(DIM, e[i]); }

	//QThis_Eigenvalues_A(DIM, lambda, A);
	//QThis_Eigenvectors_A(DIM, lambda, e, A);

	delete[] lambda;
	for (int i = 0; i < DIM; i++) { delete[] A[i]; delete[] e[i];}
	delete[] A; delete[] e;

	cin.get();
	return 0;
}

mytipe* eigenvalues(const int N, mytipe** A)
{
	mytipe* lambda = new mytipe[N];
	mytipe sigma;
	mytipe**E=new mytipe*[N];
	mytipe**Ak = new mytipe*[N];
	mytipe**Qk = new mytipe*[N];
	mytipe**Rk = new mytipe*[N];
	mytipe* vspom = new mytipe [N-1];
	int iter = 0;
	int I = 0;
	for (int i = 0; i < N; i++) { E[i] = new mytipe[N]; Ak[i] = new mytipe[N]; Qk[i] = new mytipe[N]; Rk[i] = new mytipe[N];}; 
	singlematr(N, E);
	copy(N, Ak, A);
	Hessenberg_matrix(N, Ak);
	cout << "H = "; vivod(N, Ak);
	for (int n = N - 1,m=N; n > 0; n--,m--)
	{ 
		iter = 0;
		do
		{
				sigma = Ak[n][n];
				copy(m, Ak, sum_matr(m, Ak, mult_matrnum(m, E, -sigma)));//A-�����*� (�����)
				QR(m, Qk, Rk, Ak);
				//copy(m, Ak, multimatrix(m, Rk, Qk));
				copy(m, Ak, sum_matr(m, multimatrix(m, Rk, Qk), mult_matrnum(m, E, sigma)));//RkQk+�����*� (�������� �����)
				//if (iter % 5 == 0) { cout << "A="; vivod(N, Ak); }
				if (iter > 200) { cout << "������������!!!\n"; break; }
			    iter++;
			   //for (int k = 0; k < n; k++) { vspom[k] = Ak[k + 1][k]; };
		} while (fabs(Ak[n][n - 1])>EPS_eigen);
		lambda[n] = Ak[n][n];
		I += iter;
		cout << N-n << ") " << iter << " ��������\n";
	}
	cout << N << ") " << 0 << " ��������\n����� ";
	cout << I << " ��������\n";
	lambda[0] = Ak[0][0];

	for (int i = 0; i < N; i++)
	{
	delete[] E[i];
	delete[] Ak[i];
	delete[] Qk[i];
	delete[] Rk[i];}
	delete[] E; delete[] Ak; delete[] Qk; delete[] Rk;
	return lambda;
}

void Hessenberg_matrix(const unsigned int DIM, mytipe** A)
{
	mytipe vspom;   //��������������� ��������� ��� ����� ���������
	mytipe c;  //����������� ����������� �_ij
	mytipe s;  //����������� ����������� �_ij

	for (int k = 0; k < DIM - 2; k++)  //���������� ������ T=Q^(-1) � R
									   //������������ ����� T_ij, ��������� �������
	{
		for (int l = k + 2; l < DIM; l++)  //������������ ����� T_ij, ��������� ������ ��� ������� ����������
		{

			c = A[k + 1][k] / sqrt((A[k + 1][k])*(A[k + 1][k]) + (A[l][k])*(A[l][k]));//���������� c
			s = A[l][k] / sqrt((A[k + 1][k])*(A[k + 1][k]) + (A[l][k])*(A[l][k]));//���������� s

			for (int j = k; j < DIM; j++)
			{
				vspom = A[k + 1][j];
				A[k + 1][j] = c * A[k + 1][j] + s * A[l][j]; //��������� ���� ����� �
				A[l][j] = c * A[l][j] - s * vspom;    //�� ������� �
			}
			for (int j = k; j < DIM; j++)
			{
				vspom = A[j][k + 1];
				A[j][k + 1] = c * A[j][k + 1] + s * A[j][l]; //��������� ���� �������� �
				A[j][l] = c * A[j][l] - s * vspom;    //�� ������ �
			}
		}
	}
}

mytipe** eigenvectors(const int N, mytipe** A,mytipe* lambda)
{
	mytipe** e = new mytipe*[N];//������� ����������� ��������
	mytipe** Ak = new mytipe*[N];//������� �-�������
	mytipe** E = new mytipe*[N];//
	mytipe* e_old= new mytipe [N];
	int iter = 0;
	int I = 0;
	for (int i = 0; i < N; i++) { e[i] = new mytipe[N]; Ak[i] = new mytipe[N]; E[i] = new mytipe[N]; } singlematr(N, E); singlematr(N, e);
	for (int j = 0; j < N; j++)
	{   
		iter = 0;
	copy(N,Ak,sum_matr(N, A, mult_matrnum(N, E, -lambda[j])));
	do
	{
		copy(N, e_old, e[j]);
		findx(N, e[j], Ak);//������� ����
		hod(N, e[j], Ak);
		copy(N, e[j], mult_vectnum(N, e[j], 1 / norm(N, e[j], '2')));//����������	
		iter++;
	} while (fabs(fabs(skalar(N, e[j], e_old)) - 1)>EPS_eigen);//(norm(N, e_old, '2') > EPSILON); //(norm(N, matrvec(N, Ak, e[j]), '2')>EPS);
	cout << j + 1 << ") " << iter << " ��������\n";
	I += iter;
	}
	cout << "����� ��������:" << I<<endl;
	delete[] e_old;
	for (int i = 0; i < N; i++) { delete[] Ak[i];  delete[] E[i]; };
	delete[] E; delete[] Ak;

	return e;
}

void Raley(int N, mytipe* lambda, mytipe** e,mytipe** A)
{
	mytipe** Ak = new mytipe*[N];//������� �-�������
	mytipe** E = new mytipe*[N];//
	mytipe* vspom = new mytipe[N];
	int iter = 0;
	int I = 0;
	for (int i = 0; i < N; i++) { e[i] = new mytipe[N]; Ak[i] = new mytipe[N]; E[i] = new mytipe[N]; } 
	singlematr(N, E); singlematr(N, e);
	for (int j = 0; j < N; j++)
	{
		
		for (int i = 0; i < j; i++) {
			copy(N, vspom, mult_vectnum(N, proj(N, e[i], e[j]), -1));
			copy(N, e[j], sum_vect(N, e[j],vspom));	}//���������������	������-������ 

		copy(N, e[j], mult_vectnum(N, e[j], 1 / norm(N, e[j], '2')));//����������
		iter = 0;
		do
		{
			copy(N, vspom, e[j]);
			lambda[j] = skalar(N, matrvec(N, A, e[j]), e[j]);
			copy(N, Ak, sum_matr(N, A, mult_matrnum(N, E, -lambda[j])));
			findx(N, e[j], Ak);//������� ����
			hod(N, e[j], Ak);
			copy(N, e[j], mult_vectnum(N, e[j], 1 / norm(N, e[j], '2')));//����������	
			iter++;
		} while (fabs(fabs(skalar(N, e[j], vspom)) - 1)>EPS_eigen);//(norm(N, matrvec(N, Ak, e[j]), '2')>EPS);
		cout << j + 1 << ") " << iter << " ��������\n";
		I += iter;
	}

	cout << "����� ��������:" << I << endl;
	delete[] vspom;
	for (int i = 0; i < N; i++) { delete[] Ak[i];  delete[] E[i]; };
	delete[] E; delete[] Ak;
}

bool QThis_Eigenvalues_A(int N, mytipe* lambda, mytipe** A)
{
	mytipe** Ak = new mytipe*[N];//������� �-�������
	mytipe** E = new mytipe*[N];//
	mytipe sym = 0.0;
	for (int i = 0; i < N; i++) { Ak[i] = new mytipe[N]; E[i] = new mytipe[N]; };
	singlematr(N, E);
	for (int k = 0; k < N; k++) {
		for (int j = 0; j < N; j++) {
			sym = 0.0;
			copy(N, Ak, sum_matr(N, A, mult_matrnum(N, E, -lambda[k])));
			for (int i = 0; i < j; i++) { sym += fabs(Ak[j][i]); };
			for (int i = j + 1; i < N; i++) { sym += fabs(Ak[j][i]); };
			if (fabs(Ak[j][j]) - sym > EPS_eigen) {
				cout << endl << j + 1 << " �������� �� �������� �����������.\n" ;
				for (int i = 0; i < N; i++) { delete[] Ak[i];  delete[] E[i]; };
				delete[] E; delete[] Ak;
				return false;
			};
		};
	};
	for (int i = 0; i < N; i++) { delete[] Ak[i];  delete[] E[i]; };
	delete[] E; delete[] Ak;
	cout << "\n��� �������� �������� ������������.\n";
	return true;
}

bool QThis_Eigenvectors_A(int N, mytipe* lambda, mytipe** e, mytipe** A)
{
	mytipe** Ak = new mytipe*[N];//������� �-�������
	mytipe** E = new mytipe*[N];//
	for (int i = 0; i < N; i++) { Ak[i] = new mytipe[N]; E[i] = new mytipe[N]; };
	singlematr(N, E);
	for (int j = 0; j < N; j++) {
		copy(N, Ak, sum_matr(N, A, mult_matrnum(N, E, -lambda[j])));
		if(norm(N, matrvec(N, Ak, e[j]), '2')>10*EPS_eigen){
			cout << endl << j + 1 << " ������ �� �������� �����������.\n";
			for (int i = 0; i < N; i++) { delete[] Ak[i];  delete[] E[i]; };
			delete[] E; delete[] Ak;
			return false;
		};
	};

	for (int i = 0; i < N; i++) { delete[] Ak[i];  delete[] E[i]; };
	delete[] E; delete[] Ak;
	cout << "\n��� ������� �������� ������������.\n";
	return true;
}
