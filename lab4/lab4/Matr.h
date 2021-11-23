#pragma once//��������� ��� ����, ����� ���� ����� ������ 1 ���

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


void copy(const unsigned int DIM, mytipe** B, mytipe** A) //����������� ������
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

mytipe** multimatrix(const unsigned int DIM, mytipe** const a, mytipe** const b) //������������ ������
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
			t = 0.0;
			for (int j = 0; j<DIM; j++) {
				t += a[i][j] * b[j][l];  //������������ ������ ������� a �� ������� ������� b
			}
			E[i][l] = t;   //������ ���������� ������������ � �������������� �������
		}
	}

	return E;
}

mytipe** mult_matrnum(const unsigned int DIM, mytipe** const A, const mytipe a) //�������� ������� �� �����
{
	mytipe** E = new mytipe *[DIM];
	for (int i = 0; i < DIM; i++)       // ��� �������
	{                                   //  ���������� �����������
		E[i] = new mytipe[DIM];         //  ������������ ������ a � b
	}

	for (int i = 0; i<DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			E[i][j] = a * A[i][j];
		}
	}
	return E;
}


mytipe* matrvec(const unsigned int DIM, mytipe** A, mytipe* b) //��������� ������� �� �������
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
	return c;
}



mytipe* mult_vectnum(const unsigned int DIM, const mytipe* const b, const mytipe a) //��������� ������� �� �����
{
	mytipe* c;
	c = new mytipe[DIM];
	for (int i = 0; i<DIM; i++) {
		c[i] = a * b[i];
	}
	return c;
}

mytipe* sum_vect(const unsigned int DIM, const mytipe* const b1, const mytipe* const b2) //�������� ��������
{
	mytipe* c;
	c = new mytipe[DIM];
	for (int i = 0; i < DIM; i++) {
		c[i] = b1[i] + b2[i];
	}
	return c;
}

mytipe** sum_matr(const unsigned int DIM, mytipe** A1, mytipe** A2) //�������� ������
{
	mytipe** E;
	E = new mytipe *[DIM];              // �������� �����
	for (int i = 0; i < DIM; i++)       // ��� �������
	{                                   //  ���������� �����������
		E[i] = new mytipe[DIM];         //  ������������ ������ a � b
	}
	for (int i = 0; i<DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			E[i][j] = A1[i][j] + A2[i][j];
		}
	}
	return E;
}



void singlematr(const unsigned int DIM, mytipe** E) //���������� ������� ���������
{
	for (int i = 0; i < DIM; i++)      //|������
	{									//|������� �
		for (int j = 0; j < DIM; j++)  //|���������
		{							    //|���
			E[i][j] = 0.;				//|����������
		}								//|��������������
		E[i][i] = 1.;					//|
	}
}



void hod(const unsigned int DIM, mytipe* b, mytipe** a) //�������� ��� ������
{
	for (int i = (DIM - 1); i >= 0; i--)
	{
		for (int j = i + 1; j < DIM; j++)
		{
			b[i] = b[i] - a[i][j] * b[j]; //��������� �� ������� ������ ����� ���� ��������� ������� ���� �� ������ ����� �������� �� �������� ���������
		}
		b[i] = b[i] / a[i][i]; //������� �������� ������� ������ ����� �� �������� ������� ���������
	}
}

bool findx(const unsigned int DIM, mytipe* b, mytipe** a) //������ ��� ������ � ��������� �� �������������
{
	mytipe eps;     //��������������� ���������� ��� ��������� ��������� (���������� max)
	mytipe* vspom;  //��������������� ������ ��� ������������ ����� �������
	mytipe vspom2;   //��������������� ���������� ��� ������������ ��������� �������
	bool flag = true;  //���������� ����������(�������� ������������� �������)
	mytipe EPS = 1e-20;

	for (int k = 0, h; k < (DIM); k++)
	{
		eps = abs(a[k][k]);  //������ ������������� �������� � ���������� ���������
		vspom = a[k];   //������ ������ �� ������ ������������� ��������(���� ������� � ����� �� ���� �� ����������)
		h = k;   //������ ������ ������ (�� �������, ��������� ����)

		if (abs(a[k][k]) < EPS) { flag = true; } //���� ������ ������� �������, �� ������������ ����
		else { flag = false; } //���� �� �������, �� �������������� ����

		for (int i = (k + 1); i < DIM; i++)
		{
			if (abs(a[i][k]) > eps && abs(a[i][k])>EPS) //���� ������� ������ ����������� max � ������� � �� ����
			{
				eps = a[i][k]; vspom = a[i];  h = i; flag = false;
			}; //��������� ������, �������� max � �������������� ����
		}

		if (!flag) //���� ������� �� �������
		{
			a[h] = a[k];  //����� ��������
			a[k] = vspom;
			vspom2 = b[h]; //����� ���������� �������
			b[h] = b[k];
			b[k] = vspom2;
		}
		else { cout << "������� ���������\n";   return flag; }

		for (int i = k + 1; i < DIM; i++)
		{
			b[i] = b[i] - b[k] / a[k][k] * a[i][k]; //��������� �� ��������� �������
			for (int j = (DIM - 1); j >= 0; j--)
			{
				a[i][j] = a[i][j] - a[k][j] / a[k][k] * a[i][k]; //��������� ���� ��������� ��� ������������
			}
		}

	}
	return flag;
}

void QR(const unsigned int DIM, mytipe** T, mytipe** R, mytipe** A) //QR-��������
{
	mytipe vspom;   //��������������� ��������� ��� ����� ���������
	mytipe c;  //����������� ����������� �_ij
	mytipe s;  //����������� ����������� �_ij
	copy(DIM, R, A);
	singlematr(DIM, T);
	for (int k = 0; k < DIM; k++)  //���������� ������ T=Q^(-1) � R
								   //������������ ����� T_ij, ��������� �� ������� ���������
	{
		for (int i = k + 1; i < DIM; i++)  //������������ ����� T_ij, ������� ������ c ���������� �������
		{
			if (abs(R[i][k]) > EPSILON)
			{
				c = R[k][k] / sqrt((R[k][k])*(R[k][k]) + (R[i][k])*(R[i][k]));//���������� c
				s = R[i][k] / sqrt((R[k][k])*(R[k][k]) + (R[i][k])*(R[i][k]));//���������� s
				for (int j = 0; j < DIM; j++)  //��������� T_ij �� A
				{
					vspom = T[k][j];
					T[k][j] = c * T[k][j] + s * T[i][j]; //��������� ���� ����� �
					T[i][j] = c * T[i][j] - s * vspom;  //�� ������� ���������� �
				}
				for (int j = k; j < DIM; j++)
				{
					vspom = R[k][j];
					R[k][j] = c * R[k][j] + s * R[i][j]; //��������� ���� ����� �
					R[i][j] = c * R[i][j] - s * vspom;    //�� ������� �
				}
			}

		}
	}
	trans(DIM, T);//������� Q
}

mytipe skalar(int N, mytipe* a, mytipe* b)
{
	mytipe r=0.0;
	for (int i = 0; i < N; i++) { r += a[i] * b[i]; };
	return r;
}

mytipe* proj(int N, mytipe* b, mytipe* a)
{
	return mult_vectnum(N, b, skalar(N, a, b) / skalar(N, b, b));
}