#include "/home/yahya/cuda-workspace/Uploaded/C_function.h"

//using namespace std;

void Print_int(int A[], int width) {
	for (int i = 0; i < width; ++i) {
		cout << endl << i << "  : " << A[i];
	}
	cout << endl;
}
void Print(float A[], int width) {
	for (int i = 0; i < width; ++i) {
		cout << endl << i << "  : " << A[i];
	}
	cout << endl;
}

void Print2(float A[], int height, int width) {
	for (int i = 0; i < height; ++i) {
		cout << endl << i << " :  ";
		for (int j = 0; j < width; ++j) {
			cout << A[i * width + j] << "  ";
		}

	}
	cout << endl;
}
float Volume(float x[4], float y[4], float z[4]) {
	float V;
	V = (x[0] * y[2] * z[1] - x[0] * y[1] * z[2] + x[1] * y[0] * z[2]
			- x[1] * y[2] * z[0] - x[2] * y[0] * z[1] + x[2] * y[1] * z[0]
			+ x[0] * y[1] * z[3] - x[0] * y[3] * z[1] - x[1] * y[0] * z[3]
			+ x[1] * y[3] * z[0] + x[3] * y[0] * z[1] - x[3] * y[1] * z[0]
			- x[0] * y[2] * z[3] + x[0] * y[3] * z[2] + x[2] * y[0] * z[3]
			- x[2] * y[3] * z[0] - x[3] * y[0] * z[2] + x[3] * y[2] * z[0]
			+ x[1] * y[2] * z[3] - x[1] * y[3] * z[2] - x[2] * y[1] * z[3]
			+ x[2] * y[3] * z[1] + x[3] * y[1] * z[2] - x[3] * y[2] * z[1]) / 6;
	return V;

}

void Jacobi(float x[4], float y[4], float z[4], float Vol, float J[16]) {
	J[0] = (x[1] * y[2] * z[3] - x[1] * y[3] * z[2] - x[2] * y[1] * z[3]
			+ x[2] * y[3] * z[1] + x[3] * y[1] * z[2] - x[3] * y[2] * z[1])
			/ (6 * Vol);
	J[1] = (y[2] * z[1] - y[1] * z[2] + y[1] * z[3] - y[3] * z[1] - y[2] * z[3]
			+ y[3] * z[2]) / (6 * Vol);
	J[2] = (x[1] * z[2] - x[2] * z[1] - x[1] * z[3] + x[3] * z[1] + x[2] * z[3]
			- x[3] * z[2]) / (6 * Vol);
	J[3] = (x[2] * y[1] - x[1] * y[2] + x[1] * y[3] - x[3] * y[1] - x[2] * y[3]
			+ x[3] * y[2]) / (6 * Vol);
	J[4] = (x[0] * y[3] * z[2] - x[0] * y[2] * z[3] + x[2] * y[0] * z[3]
			- x[2] * y[3] * z[0] - x[3] * y[0] * z[2] + x[3] * y[2] * z[0])
			/ (6 * Vol);
	J[5] = (y[0] * z[2] - y[2] * z[0] - y[0] * z[3] + y[3] * z[0] + y[2] * z[3]
			- y[3] * z[2]) / (6 * Vol);
	J[6] = (x[2] * z[0] - x[0] * z[2] + x[0] * z[3] - x[3] * z[0] - x[2] * z[3]
			+ x[3] * z[2]) / (6 * Vol);
	J[7] = (x[0] * y[2] - x[2] * y[0] - x[0] * y[3] + x[3] * y[0] + x[2] * y[3]
			- x[3] * y[2]) / (6 * Vol);
	J[8] = (x[0] * y[1] * z[3] - x[0] * y[3] * z[1] - x[1] * y[0] * z[3]
			+ x[1] * y[3] * z[0] + x[3] * y[0] * z[1] - x[3] * y[1] * z[0])
			/ (6 * Vol);
	J[9] = (y[1] * z[0] - y[0] * z[1] + y[0] * z[3] - y[3] * z[0] - y[1] * z[3]
			+ y[3] * z[1]) / (6 * Vol);
	J[10] = (x[0] * z[1] - x[1] * z[0] - x[0] * z[3] + x[3] * z[0] + x[1] * z[3]
			- x[3] * z[1]) / (6 * Vol);
	J[11] = (x[1] * y[0] - x[0] * y[1] + x[0] * y[3] - x[3] * y[0] - x[1] * y[3]
			+ x[3] * y[1]) / (6 * Vol);
	J[12] = (x[0] * y[2] * z[1] - x[0] * y[1] * z[2] + x[1] * y[0] * z[2]
			- x[1] * y[2] * z[0] - x[2] * y[0] * z[1] + x[2] * y[1] * z[0])
			/ (6 * Vol);
	J[13] = (y[0] * z[1] - y[1] * z[0] - y[0] * z[2] + y[2] * z[0] + y[1] * z[2]
			- y[2] * z[1]) / (6 * Vol);
	J[14] = (x[1] * z[0] - x[0] * z[1] + x[0] * z[2] - x[2] * z[0] - x[1] * z[2]
			+ x[2] * z[1]) / (6 * Vol);
	J[15] = (x[0] * y[1] - x[1] * y[0] - x[0] * y[2] + x[2] * y[0] + x[1] * y[2]
			- x[2] * y[1]) / (6 * Vol);
}

void shape_function(float B[ele_num * 6 * 12], float J[ele_num][4][4]) {
	int row_B = 6;
	int col_B = 12;
	for (int k = 0; k < ele_num; k++) {
		//B[col_B*row_B*k+col_B*i+j] [k][i][j]
		B[col_B * row_B * k + col_B * 0 + 0] = J[k][0][1];
		B[col_B * row_B * k + col_B * 0 + 3] = J[k][1][1];
		B[col_B * row_B * k + col_B * 0 + 6] = J[k][2][1];
		B[col_B * row_B * k + col_B * 0 + 9] = J[k][3][1];

		B[col_B * row_B * k + col_B * 1 + 1] = J[k][0][2];
		B[col_B * row_B * k + col_B * 1 + 4] = J[k][1][2];
		B[col_B * row_B * k + col_B * 1 + 7] = J[k][2][2];
		B[col_B * row_B * k + col_B * 1 + 10] = J[k][3][2];

		B[col_B * row_B * k + col_B * 2 + 2] = J[k][0][3];
		B[col_B * row_B * k + col_B * 2 + 5] = J[k][1][3];
		B[col_B * row_B * k + col_B * 2 + 8] = J[k][2][3];
		B[col_B * row_B * k + col_B * 2 + 11] = J[k][3][3];

		B[col_B * row_B * k + col_B * 3 + 0] = J[k][0][2];
		B[col_B * row_B * k + col_B * 3 + 3] = J[k][1][2];
		B[col_B * row_B * k + col_B * 3 + 6] = J[k][2][2];
		B[col_B * row_B * k + col_B * 3 + 9] = J[k][3][2];

		B[col_B * row_B * k + col_B * 3 + 1] = J[k][0][1];
		B[col_B * row_B * k + col_B * 3 + 4] = J[k][1][1];
		B[col_B * row_B * k + col_B * 3 + 7] = J[k][2][1];
		B[col_B * row_B * k + col_B * 3 + 10] = J[k][3][1];

		B[col_B * row_B * k + col_B * 4 + 1] = J[k][0][3];
		B[col_B * row_B * k + col_B * 4 + 4] = J[k][1][3];
		B[col_B * row_B * k + col_B * 4 + 7] = J[k][2][3];
		B[col_B * row_B * k + col_B * 4 + 10] = J[k][3][3];

		B[col_B * row_B * k + col_B * 4 + 2] = J[k][0][2];
		B[col_B * row_B * k + col_B * 4 + 5] = J[k][1][2];
		B[col_B * row_B * k + col_B * 4 + 8] = J[k][2][2];
		B[col_B * row_B * k + col_B * 4 + 11] = J[k][3][2];

		B[col_B * row_B * k + col_B * 5 + 0] = J[k][0][3];
		B[col_B * row_B * k + col_B * 5 + 3] = J[k][1][3];
		B[col_B * row_B * k + col_B * 5 + 6] = J[k][2][3];
		B[col_B * row_B * k + col_B * 5 + 9] = J[k][3][3];

		B[col_B * row_B * k + col_B * 5 + 2] = J[k][0][1];
		B[col_B * row_B * k + col_B * 5 + 5] = J[k][1][1];
		B[col_B * row_B * k + col_B * 5 + 8] = J[k][2][1];
		B[col_B * row_B * k + col_B * 5 + 11] = J[k][3][1];
	}
}

void transpose(float A[], int row_A, int col_A, float A_T[]) {
	for (int k = 0; k < ele_num; k++) {
		for (int i = 0; i < col_A; i++) {
			for (int j = 0; j < row_A; j++) {
				A_T[row_A * col_A * k + row_A * i + j] = A[row_A * col_A * k
						+ col_A * j + i];
			}
		}
	}
}

void Mass(vector<float>* M, const float* V, const int* ID_ele) {
	int width = 4 * DOF;
	int pos = 0;
	float value;
	for (int i = 0; i < ele_num; i++) {
		float m2[12 * 12] = { };
		for (int j = 0; j < 12; j++) {
			m2[j * 12 + j] = Ro * V[i] / 4;
		}
		for (int j = 0; j < DOF * 4; j++) {
			for (int k = 0; k < DOF * 4; k++) {
				if (ID_ele[i * width + j] != 0 && ID_ele[i * width + k] != 0) {
					pos = (ID_ele[i * width + j] - 1) * length
							+ (ID_ele[i * width + k] - 1);
					value=(*M).at(pos)+ m2[j * 12 + k];
					(*M).at(pos)=value;
				}
			}
		}
	}
}

void Stiffness(const float* B_ele, const float* B_eleT, const float C[6][6],
		const int* ID_ele, const float* V, float* K) {
	int row_B_ele = 6;
	int col_B_ele = 12;
	int row_B_eleT = 12;
	int col_B_eleT = 6;
	int col_ID_ele = 4 * DOF;
	int pos=0;
		float K1[12][6];
		float K2[12][12];
	for (int i = 0; i < ele_num; i++) {
		float sum = 0;
		for (int i2 = 0; i2 < 12; i2++) {
			for (int j2 = 0; j2 < 6; j2++) {
				for (int k2 = 0; k2 < 6; k2++) {
					pos = col_B_eleT * row_B_eleT * i + col_B_eleT * i2 + k2;
					sum += B_eleT[pos] * C[k2][j2]; // *(*(*(B_minT + i2) + k2))*(*(*(C + k2) + j2));
				}
				K1[i2][j2] = sum;
				sum = 0;
			}
		}
		sum = 0;
		for (int i2 = 0; i2 < 12; i2++) {
			for (int j2 = 0; j2 < 12; j2++) {
				for (int k2 = 0; k2 < 6; k2++) {
					pos = col_B_ele * row_B_ele * i + col_B_ele * k2 + j2;
					sum += K1[i2][k2] * B_ele[pos] * V[i];
				}
				K2[i2][j2] = sum;
				sum = 0;
			}
		}

		for (int n = 0; n < DOF * 4; n++) {
			for (int n1 = 0; n1 < DOF * 4; n1++) {
				if (ID_ele[i * col_ID_ele + n] != 0
						&& ID_ele[i * col_ID_ele + n1] != 0)
					pos = (ID_ele[i * col_ID_ele + n] - 1) * length
							+ (ID_ele[i * col_ID_ele + n1] - 1);
				K[pos] += K2[n][n1];
			}
		}
	}

	int counter2 = 0;
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < length; j++) {
			if (K[i * length + j] != 0) {
				counter2 += 1;
			}
		}
	}
}


void Sparsize(vector<float>* A, const int width_A, float* A_sp, int* RA_sp, int* CK_sp,
		int* non_zero) {
	RA_sp[0] = 0;
	int counter = 0;
	for (int i = 0; i < width_A; i++) {
		for (int j = 0; j < width_A; j++) {
			if (A->at(i * width_A + j) != 0) {
				A_sp[counter] = A->at(i * width_A + j);
				CK_sp[counter] = j;
				counter += 1;
			}
			RA_sp[i + 1] = counter;
		}
	}
	*non_zero = counter;
	cout << counter << "\n";
}

void Sparsize_reduced(vector<float>* B_dyn, float* B_dyn_sp,
		const int* RB_dyn_sp, const int* CB_dyn_sp) {
	for (int i = 0; i < length; i++) {
		for (int k = RB_dyn_sp[i]; k < RB_dyn_sp[i + 1]; k++) {
			B_dyn_sp[k] = B_dyn->at(i * length + CB_dyn_sp[k]);
		}
	}
	//cout<<"Done!"<<endl;
}

void SpVec(const float* A_sp,const int* RA_sp,const int* CA_sp, const float* i_vec,
		int length_vec, float* f_vec) {
	//returns A*i_vec=f_vec

	float sum = 0;
	for (int i = 0; i < length_vec; i++) {
		for (int j = RA_sp[i]; j < RA_sp[i + 1]; j++) {
			sum += i_vec[CA_sp[j]] * A_sp[j];
		}
		//cout<<endl<<"sum= "<<sum;
		f_vec[i] = sum;
		sum = 0;
	}
}


//====================================================================
void Solve(float A_sp[], int Col_A_sp[], int Row_A_sp[], float B[], float U[]) {
	float p[length] = { };
	float r[length] = { };
	float t[length] = { };
	float rho = 0;
	float rhos;
	float alpha;
	float sum = 0;
	bool solved = 0;

	for (int i = 0; i < length; i++) {
		for (int k = Row_A_sp[i]; k < Row_A_sp[i + 1]; k++) {
			sum += (U[Col_A_sp[k]]) * (A_sp[k]);
		}
		r[i] = B[i] - sum;
		p[i] = r[i];
		sum = 0;
	}
	for (int i = 0; i < length; i++) {
		rho += r[i] * r[i];
	}
	for (int j = 0; j < length; j++) {
		if (solved == 0) {
			sum = 0;
			for (int i = 0; i < length; i++) {
				for (int k = Row_A_sp[i]; k < Row_A_sp[i + 1]; k++) {
					sum += (p[Col_A_sp[k]]) * (A_sp[k]);
				}
				t[i] = sum;
				sum = 0;
			}

			float PT = 0;
			for (int i = 0; i < length; i++) {
				PT += p[i] * t[i];
			}
			alpha = rho / PT;
			for (int i = 0; i < length; i++) {
				U[i] += alpha * p[i];
				r[i] -= alpha * t[i];
			}

			rhos = rho;
			rho = 0;
			for (int i = 0; i < length; i++) {
				rho += r[i] * r[i];
			}
			if ((rho) < error) {
				solved = 1;
				//cout << endl << "Solved in " << j << " Steps!" << "\n";
			}
			for (int i = 0; i < length; i++) {
				p[i] = r[i] + (rho / rhos) * p[i];
			}

		}
	}
	//cout << "rho    " << rho << "\n";
	//cout << "U(56)= " << *(U + 56) << "\n" << "\n";
}

void Solve_diag(float A_sp[], float B[], float U[]) {
	for (int i = 0; i < length; ++i) {
		U[i] = B[i] / A_sp[i];
	}

	//cout << "rho    " << rho << "\n";
	//cout << "U(56)= " << *(U + 56) << "\n" << "\n";
}

//void Multi(float x[], int Ix, int Jx, float y[], int Iy, int Jy, float z[]) {
//	//Rerurns x*y & Jx=Iy
//	float sum = 0;
//	for (int i = 0; i < Ix; i++) {
//		for (int j = 0; j < Jy; j++) {
//			for (int k = 0; k < Iy; k++) {
//				sum += x[i * Jx + k] * y[k * Jy + j];
//			}
//			z[i * Iy + j] = sum;
//			sum = 0;
//		}
//	}
//}

//void Change_B(float B_dyn[length * length], float B_dyn_sp[sparse_data],
//		float M[length * length], float K_mu[length * length],
//		float K_la[length * length], float la, float mu) {
//	for (int i = 0; i < length; i++) {
//		for (int j = 0; j < length; j++) {
//			B_dyn[i * length + j] = 2 * M[i * length + j]
//					- (la * K_la[i * length + j] + mu * K_mu[i * length + j])
//							* dt * dt;
//		}
//	}
//	int Column[length];
//	int Row[length + 1];
//	int sss = 0;
//	Sparsize(B_dyn, length, B_dyn_sp, Row, Column, &sss);
//}
