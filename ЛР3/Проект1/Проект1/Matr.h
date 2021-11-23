#pragma once//��������� ��� ����, ����� ���� ����� ������ 1 ���

#include <iostream>

typedef double mytipe;


using namespace std;
mytipe norm(const unsigned int DIM, mytipe** const A, const char flag) //4 ����� ������� � �������� �������� �����
{
	mytipe norm = 0;  //���������� ����� �������
	mytipe vspom = 0;  //��������������� ����������
	if (flag == '1')  //�������������� ����� ������� (max �������� �� ��������)
	{
		for (int i = 0; i < DIM; i++)
		{
			vspom = 0;
			for (int j = 0; j<DIM; j++)
			{
				vspom += abs(A[j][i]); //��������� �������� i-�� ������
			}
			if (norm < vspom) { norm = vspom; };   //��������� �� ��������
		}
		return norm;
	}

	if (flag == 'k')  //���������� ����� ������� (max �������� �� �������)
	{
		for (int j = 0; j < DIM; j++)
		{
			vspom = 0;
			for (int i = 0; i<DIM; i++)
			{
				vspom += abs(A[j][i]);  //��������� �������� j-�� ������
			}
			if (norm < vspom) { norm = vspom; };   //��������� �� ��������
		}

		return norm;
	}

	if (flag == '2')  //������� ����� ������� (������ �� ����� ��������� ���� ���������)
	{
		for (int j = 0; j < DIM; j++)
		{
			for (int i = 0; i<DIM; i++)
			{
				norm += A[j][i] * A[j][i];  //���������� ��������
			}
		}

		return sqrt(norm);  //���������� ������ �� ����� ���������
	}

	if (flag == 'm')  //������������ ����� ������� (n*max)
	{

		for (int j = 0; j < DIM; j++)
		{
			for (int i = 0; i<DIM; i++)
			{
				if (abs(A[j][i]) > norm) { norm = abs(A[j][i]); } //������� ������������ �������
			}
		}
		return DIM*norm;
	}
}

mytipe norm(const unsigned int DIM, const mytipe* const b, const char flag) //3 ����� ������� � �������� �������� �����
{
	mytipe norm = 0;  //���������� ����� �������
	if (flag == 'k')  //���������� ����� (max)
	{

		for (int i = 0; i < DIM; i++)
		{
			if (norm < abs(b[i])) { norm = abs(b[i]); }  //������� ��������
		}
		return norm;
	}

	if (flag == '1') //�������������� ����� (����� ���������)
	{
		for (int j = 0; j < DIM; j++)
		{
			norm += abs(b[j]);  //���������� �������� ���������
		}
		return norm;
	}

	if (flag == '2')  //������� �����  (���������)
	{
		for (int i = 0; i<DIM; i++)
		{
			norm += b[i] * b[i];  //��������� �������� ���������
		}

		return sqrt(norm); //���������� ������ ���������� �� ����� ���������
	}
}

void vivod(const unsigned int DIM, const mytipe* const b, mytipe** const a) //��������� ������ ������� � �������
{
	cout << endl;
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			cout << a[i][j] << " ";
		}
		cout << '|' << b[i] << endl;
	}
	cout << endl;
}

void vivod(const unsigned int DIM, mytipe** const a) //��������� ������ �������
{
	cout << endl;
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			cout << a[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void vivod(const unsigned int DIM, const mytipe* const b) //��������� ������ �������
{
	cout << '(';
	for (int i = 0; i < DIM; i++)
	{
		cout << b[i];
		if (i < DIM - 1) { cout << ','; } //�������� ������
	}
	cout << ")^�\n";
}




mytipe neviaz(const unsigned int DIM, const mytipe* const x, mytipe** const A, const mytipe* const b, const char k) //�������
{
	mytipe* b1;          //��������� ������ ���
	b1 = new mytipe[DIM];  //������, ������� ��������� ��� ����������� ������� � �������
	mytipe nor;  //����� ������� ��������

	for (int i = 0; i < DIM; i++)
	{
		b1[i] = 0;
		for (int j = 0; j < DIM; j++)
		{
			b1[i] += A[i][j] * x[j]; //�������, ������������ �������, ������� ������ �����

		}
		b1[i] -= b[i]; //�������� i-�� ���������� ������������� ������� � ���������
	}

	nor = norm(DIM, b1, k);
	delete[] b1; //������������ ������ ������������� �������
	return nor;
}


void copy(const unsigned int DIM, mytipe** B, mytipe** const A) //����������� ������
{

	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			B[i][j] = A[i][j];
		}
	}
}

void copy(const unsigned int DIM, mytipe* b, const mytipe* a) //����������� ��������
{
	for (int j = 0; j < DIM; j++)
	{
		b[j] = a[j];
	}
}

void trans(const unsigned int DIM, mytipe** a) //���������������� �������
{
	mytipe vspom;
	for (int i = 0; i < DIM; i++)
	{

		for (int j = i + 1; j < DIM; j++)
		{
			vspom = a[i][j];
			a[i][j] = a[j][i];
			a[j][i] = vspom;
		}
	}
}

void multimatrix(const unsigned int DIM, mytipe** const a, mytipe** const b) //������������ ������
{
	double t;
	mytipe** E;
	E = new mytipe *[DIM];              // �������� �����
	for (int i = 0; i < DIM; i++)       // ��� �������
	{                                   //  ���������� �����������
		E[i] = new mytipe[DIM];         //  ������������ ������ a � b
	}

	for (int i = 0; i<DIM; i++) {
		for (int l = 0; l<DIM; l++) {
			t = 0;
			for (int j = 0; j<DIM; j++) {
				t += a[i][j] * b[j][l];  //������������ ������ ������� a �� ������� ������� b
			}
			E[i][l] = t;   //������ ���������� ������������ � �������������� �������
		}
	}

	vivod(DIM, E);  //������� ���������

	for (int i = 0; i < DIM; i++) {   //������� ������������ ������ �������
		delete[] E[i];
	}
	delete[] E;             // ������� ������������ �������
}

void mult_matr(const unsigned int DIM, mytipe** const A1, mytipe** const A2, mytipe** X) //������������ ������
{
	for (int i = 0; i<DIM; i++) {
		for (int l = 0; l<DIM; l++) {
			X[i][l] = 0;
			for (int j = 0; j<DIM; j++) {
				X[i][l] += A1[i][j] * A2[j][l];  //������������ ������ ������� a �� ������� ������� b
			}
		}
	}
}

void mult_matrvect(const unsigned int DIM, mytipe** const A, const mytipe* const b, mytipe* x) //��������� ������� �� ������
{
	for (int i = 0; i<DIM; i++) {
		x[i] = 0;
		for (int j = 0; j < DIM; j++) {
			x[i] += A[i][j] * b[j];
		}
	}
}

void mult_matrnum(const unsigned int DIM, mytipe** const A, const mytipe a, mytipe** X) //�������� ������� �� �����
{
	for (int i = 0; i<DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			X[i][j] = a * A[i][j];
		}
	}
}

void mult_vectnum(const unsigned int DIM, const mytipe* const b, const mytipe a, mytipe* x) //��������� ������� �� �����
{
	for (int i = 0; i<DIM; i++) {
		x[i] = a * b[i];
	}
}

void sum_vect(const unsigned int DIM, const mytipe* const b1, const mytipe* const b2, mytipe* x) //�������� ��������
{
	for (int i = 0; i < DIM; i++) {
		x[i] = b1[i] + b2[i];
	}
}

void sum_matr(const unsigned int DIM, const mytipe** A1, const mytipe** A2, mytipe** X) //�������� ������
{
	for (int i = 0; i<DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			X[i][j] = A1[i][j] + A2[i][j];
		}
	}
}



void singlematr(const unsigned int DIM, mytipe** const E) //���������� ������� ���������
{
	for (int i = 0; i < DIM; i++)      //|������
	{									//|������� �
		for (int j = 0; j < DIM; j++)  //|���������
		{							    //|���
			E[i][j] = 0;				//|����������
		}								//|��������������
		E[i][i] = 1;					//|
	}
}


void matrvec(const unsigned int DIM, mytipe** A, mytipe* b) //��������� ������� �� �������
{
	mytipe* c;
	c = new mytipe[DIM];
	for (int i = 0; i < DIM; i++)
	{
		c[i] = 0;
		for (int j = 0; j < DIM; j++)
		{
			c[i] += A[i][j] * b[j];
		}
	}
	copy(DIM, b, c);
	delete[] c;
}

