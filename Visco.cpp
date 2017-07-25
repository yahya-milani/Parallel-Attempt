//============================================================================
// Name        : Code.cpp
// Author      : Yahya Milani
// Version     :
// Copyright   : Use as much as you like with author's name for noncommercial cases ONLY
// Description : C++, Ansi-style
//============================================================================

#include "/home/yahya/cuda-workspace/Uploaded/C_function.h"

int main() {
	//Read Nodes=================================================/
	ifstream NODE;
	cout << "Nodes initiated" << endl;
	NODE.open("nodes2.txt");
	float node[node_num * DOF] = { };
	//int row_node = node_num;
	int col_node = DOF;
	int i = 0;
	float data1, data2, data3, data0;
	while (NODE >> data0 >> data1 >> data2 >> data3) {
		node[i * col_node + 0] = data1;
		node[i * col_node + 1] = data2;
		node[i * col_node + 2] = data3;
		i += 1;
	}
	NODE.close();
	cout << i << " Nodes Created successfully and the last node is ==> "
			<< node[(i - 1) * col_node + 0] << " "
			<< node[(i - 1) * col_node + 1] << " "
			<< node[(i - 1) * col_node + 2] << endl;
	//Create Nodes DOF =================================================/
	int ID_node[node_num * DOF] = { };
	int row_ID_node = node_num;
	int col_ID_node = DOF;
	int counter = 1;
	bool bound;
	for (int i = 0; i < row_ID_node; i++) {
		bound = 1;
		for (int j = 0; j < bound_num; j++) {
			if (i == boundary[j]) {
				bound = 0;
			}
		}
		if (bound == 1) {
			for (int k = 0; k < DOF; k++)
				ID_node[i * col_ID_node + k] = counter + k;
			counter += DOF;
		} else {
			for (int k = 0; k < DOF; k++)
				ID_node[i * col_ID_node + k] = 0;
		}
	}
	int LastNodeID = (row_ID_node - 1) * col_ID_node + 0;
	cout << "Nodes ID Created successfully and the last node's IDs are ==> "
			<< *(ID_node + LastNodeID) << " " << *(ID_node + LastNodeID + 1)
			<< " " << *(ID_node + LastNodeID + 2) << endl;
	//Read Elements=================================================/
	ifstream element;
	element.open("elements2.txt");
	i = 0;
	int ele[ele_num * 4] = { };
	cout << endl << "Elements initiated" << endl;
	//int row_ele = ele_num;
	int col_ele = 4;
	int data4, data5, data6, data7, data00;
	while (element >> data00 >> data4 >> data5 >> data6 >> data7) {
		ele[i * col_ele + 0] = data4 - 1;
		ele[i * col_ele + 1] = data5 - 1;
		ele[i * col_ele + 2] = data6 - 1;
		ele[i * col_ele + 3] = data7 - 1;
		i += 1;
	}
	element.close();
	cout << i << " Elements Created successfully and the last elements are ==> "
			<< ele[(i - 1) * col_ele + 0] << " " << ele[(i - 1) * col_ele + 1]
			<< " " << ele[(i - 1) * col_ele + 2] << " "
			<< ele[(i - 1) * col_ele + 3] << endl;
	//Create Elemental DOF =================================================/
	int ID_ele[ele_num * 4 * DOF] = { };
	int row_ID_ele = ele_num;
	int col_ID_ele = 4 * DOF;
	for (int i = 0; i < row_ID_ele; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < DOF; k++) {
				ID_ele[i * col_ID_ele + j * DOF + k] = ID_node[(ele[i * col_ele
						+ j]) * col_ID_node + k];
			}
		}

	}
	cout << "Elements ID Created successfully!" << endl;
	//=================================================/
	float x[4], y[4], z[4];
	float V[ele_num] = { };
	float B[ele_num][4][4] = { };
	int NegVol = 0;
	float TotVol = 0;
	//set B matrix and volume for every element
	for (int i = 0; i < ele_num; i++) {
		for (int j = 0; j < 4; j++) {
			x[j] = node[(ele[i * col_ele + j]) * col_node + 0];
			y[j] = node[(ele[i * col_ele + j]) * col_node + 1];
			z[j] = node[(ele[i * col_ele + j]) * col_node + 2];
		}

		// Compute Volume of each element so that none has negative volume and change indexes if there is one
		V[i] = Volume(x, y, z);
		if (V[i] < 0) {
			NegVol += 1;
			float temp = x[2];
			x[2] = x[3];
			x[3] = temp;
			temp = y[2];
			y[2] = y[3];
			y[3] = temp;
			temp = z[2];
			z[2] = z[3];
			z[3] = temp;
			V[i] = Volume(x, y, z);

			// Change ID after changing element indexes
			for (int i2 = 2 * DOF; i2 < 3 * DOF; i2++) {
				int temp2 = ID_ele[i * col_ID_ele + i2];
				ID_ele[i * col_ID_ele + i2] = ID_ele[i * col_ID_ele + i2 + DOF];
				ID_ele[i * col_ID_ele + i2 + DOF] = temp2;
			}
		}
		TotVol += V[i];
		float J[16];
		Jacobi(x, y, z, V[i], J);
		for (int var = 0; var < 4; ++var) {
			for (int var2 = 0; var2 < 4; ++var2) {
				B[i][var][var2] = J[var * 4 + var2];
			}

		}

	}
	cout << endl << "Negative volumes for " << NegVol
			<< " elements have been taken care of." << endl;
	cout << "Volume, Jacobian and B Computed for " << ele_num
			<< " elements and the total Volume is ==> " << TotVol << endl;

	//For end of Geometrical issues========================================================

	//=====================================================================================
	//===========================ALL Shit Happens Here=====================================
	//=====================================================================================
	float B_ele[ele_num * 6 * 12] = { };
	int row_B_ele = 6;
	int col_B_ele = 12;
	shape_function(B_ele, B);

	float B_eleT[ele_num * 12 * 6] = { };
	transpose(B_ele, row_B_ele, col_B_ele, B_eleT);

	size_t length2;
	length2 = length * length * sizeof(float);
	//length2=length*length;
	//float M[length*length] = { };
	//float *M = (float *) malloc(length2);
	vector<float> M(length * length);
	//float *M=new float[length*length];
	Mass(&M, V, ID_ele);
	cout << endl << "Mass created successfully last entry is ==> "
			<< M[length * length - 1] << endl;
	//free(M);
	float *K_la = (float *) malloc(length2);
	float *K_mu = (float *) malloc(length2);
	float C_la[6][6] = { { 1, 1, 1, 0, 0, 0 }, { 1, 1, 1, 0, 0, 0 }, { 1, 1, 1,
			0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0,
			0, 0 } };
	float C_mu[6][6] = { { 2, 0, 0, 0, 0, 0 }, { 0, 2, 0, 0, 0, 0 }, { 0, 0, 2,
			0, 0, 0 }, { 0, 0, 0, 1, 0, 0 }, { 0, 0, 0, 0, 1, 0 }, { 0, 0, 0, 0,
			0, 1 } };
	Stiffness(B_ele, B_eleT, C_la, ID_ele, V, K_la);
	cout << endl << "K lambda created successfully last entry is ==> "
			<< K_la[length * length - 1] << endl;

	Stiffness(B_ele, B_eleT, C_mu, ID_ele, V, K_mu);
	cout << "K mu created successfully last entry is ======> "
			<< K_mu[length * length - 1] << endl;

	//Dynamic Matrices ============================================================/
	vector<float> A_dyn(length * length);
	vector<float> B_dyn(length * length);
	vector<float> C_dyn(length * length);
	for (int i = 0; i < length; ++i) {
		for (int j = 0; j < length; ++j) {

			A_dyn[i * length + j] = M[i * length + j]
					+ 0.5 * dt * damp_coeff * M[i * length + j];
			B_dyn[i * length + j] = 2 * M[i * length + j]
					- (1 * K_la[i * length + j] + 1 * K_mu[i * length + j]) * dt
							* dt;
			C_dyn[i * length + j] = 0.5 * dt * damp_coeff * M[i * length + j]
					- M[i * length + j];
		}
	}

	//Sparsize K ============================================================/
//	float K_sp[sparse_data] = { }; //[51484];
//	int RK_sp[length + 1] = { };
//	int CK_sp[sparse_data] = { }; // [51484];
	int xxx = 0;
//	Sparsize(K, length, K_sp, RK_sp, CK_sp, &xxx);
//	cout << xxx << "\n";
//
	float A_dyn_sp[length] = { };
	int RA_dyn_sp[length + 1] = { };
	int CA_dyn_sp[length] = { };
	Sparsize(&A_dyn, length, A_dyn_sp, RA_dyn_sp, CA_dyn_sp, &xxx);
	cout << xxx << "\n";
//
	float B_dyn_sp[sparse_data] = { };
	int RB_dyn_sp[length + 1] = { };
	int CB_dyn_sp[sparse_data] = { };
	Sparsize(&B_dyn, length, B_dyn_sp, RB_dyn_sp, CB_dyn_sp, &xxx);
	cout << xxx << "\n";

	float C_dyn_sp[length] = { };
	int RC_dyn_sp[length + 1] = { };
	int CC_dyn_sp[length] = { };
	Sparsize(&C_dyn, length, C_dyn_sp, RC_dyn_sp, CC_dyn_sp, &xxx);
//    Print_int(RB_dyn_sp,length+1);
	cout << xxx << "\n";
//	//Force =================================
	float F[length] = { };
	F[1291] = -100; //	F[686] = -1000;

	C_dyn.clear();
	A_dyn.clear();

	clock_t start;
	clock_t start2;
	clock_t start3;
	float duration=0;
	float duration2=0;
	start = clock();
	float U1[length] = { };
	float U2[length] = { };
	float U_temp[length] = { };
	float R_side1[length] = { };
	float R_side[length] = { };
	int Step_change = dt_change / dt;
	int Step_save = dt_save / dt;
	for (int i = 0; i < tF / dt; i++) {
		if ((i % Step_save) == 0) {
			//cout << U2[1291] << endl;
		}
//		cout << "step:  " << i + 2 << endl;
//====================================================================
//		if ((i % Step_change) == 0) {
//			float mu1 = 1;
//			float la1 = 1;
//			for (int i1 = 0; i1 < length; ++i1) {
//				for (int j1 = 0; j1 < length; ++j1) {
//					B_dyn.at(i1 * length + j1) = 2 * M[i1 * length + j1]
//							- (la1 * K_la[i1 * length + j1]
//									+ mu1 * K_mu[i1 * length + j1]) * dt * dt;
//				}
//			}
//			Sparsize_reduced(&B_dyn, B_dyn_sp, RB_dyn_sp, CB_dyn_sp);
//
//		}
//====================================================================
		SpVec(B_dyn_sp, RB_dyn_sp, CB_dyn_sp, U2, length, R_side1);
		SpVec(C_dyn_sp, RC_dyn_sp, CC_dyn_sp, U1, length, R_side);

		for (int i1 = 0; i1 < length; i1++) {
			R_side[i1] += R_side1[i1] + F[i1] * dt * dt;
		}

		//Solve_diag(A_dyn_sp, R_side, U_temp);

		for (int i1 = 0; i1 < length; i1++) {
			U1[i1] = U2[i1];
			U2[i1] = U_temp[i1];
			U_temp[i1] = 0;
		}
	}
	duration = (clock() - start) / (float) CLOCKS_PER_SEC;
	cout << "All Done!" << "\n";
	cout << "time: " << duration << endl;

	return 0;
}
