/*
 * C_function.h
 *
 *  Created on: Jul 25, 2017
 *      Author: yahya
 */

#ifndef C_FUNCTION_H_
#define C_FUNCTION_H_

//#include <stdafx.h>
#include <stdio.h>
//#include <tchar.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
#include <stdlib.h>
#include <vector>
using namespace std;

//Global Constants=================================
const int ele_num = 1615; //96; //1615;
const int node_num = 496; //35; //496;
const int bound_num = 3;
const int boundary[3] = { 17, 83, 92 }; //{ 0, 3, 6, 7 }; //{ 2,3,5,6 };
const int DOF = 3;
const int length = DOF * (node_num - bound_num);
const float Ro = 1000;
const float damp_coeff = 0;
const float mu = 782.7586;
const float E = 2000;
const float v = 0.45;
const float la = 7866.7;
const float tF = 2;
const float dt = 0.001;
const float dt_change = 0.5;
const float dt_save = 0.1;
const int sparse_data = 50625;
const float error = 0.01;

void Print_int(int A[], int width);
void Print(float A[], int width);
void Print2(float A[], int height, int width);
float Volume(float x[4], float y[4], float z[4]);
void Jacobi(float x[4], float y[4], float z[4], float Vol, float J[16]);
void shape_function(float B[], float J[][4][4]);
void transpose(float A[], int row_A, int col_A, float A_T[]);
void Mass(vector<float>* M, const float* V, const int* ID_ele);
void Stiffness(const float* B_ele, const float* B_eleT,const float C[6][6],const int* ID_ele,
		const float* V, float* K);
void Solve(float A_sp[], int Col_A_sp[], int Row_A_sp[], float B[],
		float U[]);
void Sparsize(vector<float>* A, const int width_A, float* A_sp, int* RA_sp, int* CK_sp,
		int* non_zero);
void Sparsize_reduced(vector<float>* B_dyn, float* B_dyn_sp,
		const int* RB_dyn_sp, const int* CB_dyn_sp);
void SpVec(const float* A_sp,const int* RA_sp,const int* CA_sp, const float* i_vec,
		int length_vec, float* f_vec);
void Solve(float A_sp[], int Col_A_sp[], int Row_A_sp[], float B[],
		float U[]);
void Solve_diag(float A_sp[], float B[], float U[]);




#endif /* C_FUNCTION_H_ */
