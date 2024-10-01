// Код для решения статической задачи и динамической задачи с учетом сил трения 
// для 8-узлового параллелепипеда в методе конечных элементов

#include <iostream>
#include <fstream>

#include "Class_Matrix.h"

#include "Operator_overloading.h"

#include "Other_Function.h"

#include "Function&Solver_MFE.h"

#include "matrix_MFE.h"

int main()
{
	// Свойства материала
	const double modulusElastic = 7e10;
	const double coefficientPuasson = 0.33;
	const double dencity = 2700.0;
	const double time = 2.0;
	const double deltaT = 0.000001;
	// Размеры конечного элемента
	const double length = 0.5, width = 0.06, heigth = 0.05;

	// Матрица жесткости
	Matrix<double> matrixStiff = matrixStiffness(length, width, heigth, modulusElastic, coefficientPuasson);
	const size_t rowsMatrixStiffness = matrixStiff.sizeRows();

	//Статическая задача
	// Вектор узловых усилий
	std::vector<double> forceStatic(rowsMatrixStiffness);
	forceStatic[0] = 2.5e3;
	forceStatic[3] = 2.5e3;
	forceStatic[12] = 2.5e3;
	forceStatic[15] = 2.5e3;

	Matrix<double> copyMatrixStiffness = matrixStiff;

	boundConditionStatic(copyMatrixStiffness);
	//Решение для статической задачи
	std::vector<double> displacement = displacementStatic(copyMatrixStiffness, forceStatic);
	//std::cout << displacement;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Динамическая задача

	Matrix<double> matrixMass = massDiagonal(rowsMatrixStiffness, dencity, length, width, heigth);
	//std::cout << matrixMass;

	//std::vector<double> initialDisplacement(rowsMatrixStiffness);
	std::vector<double> initialDisplacement = std::move(displacement);

	/*initialDisplacement[0] = 0.0;
	initialDisplacement[3] = 0.0;
	initialDisplacement[12] = 0.0;
	initialDisplacement[15] = 0.0;*/

	/*initialDisplacement[0] = 0;
	initialDisplacement[3] = 0;
	initialDisplacement[12] = 0;
	initialDisplacement[15] = 0;*/

	std::vector<double> initialSpeed(rowsMatrixStiffness);
	initialSpeed[0] = 0.0;
	initialSpeed[3] = 0.0;
	initialSpeed[12] = 0.0;
	initialSpeed[15] = 0.0;

	std::vector<double> initialAcceliration(rowsMatrixStiffness);
	initialAcceliration[0] = 0.0;
	initialAcceliration[3] = 0.0;
	initialAcceliration[12] = 0.0;
	initialAcceliration[15] = 0.0;

	Matrix<double> displacements = displacementDinamic (matrixStiff, matrixMass, initialDisplacement,
													   initialSpeed, initialAcceliration, time, deltaT);

	const size_t rowsD = displacements.sizeRows();
	std::ofstream fout("result.txt");
	for (size_t i = 0; i < rowsD; ++i)
	{
		fout << "1: " << displacements[i][0] << "  " << "2: " << "  " <<
			displacements[i][1] << "  " << "5: " <<  displacements[i][4] <<
			"  " << "6: " << displacements[i][6] << "\n";
	}
	fout.close();

	//size_t steps = static_cast<size_t> (time / deltaT);
	//std::ofstream fout("stepstime.txt");
	//for (size_t i = 0; i < steps; ++i)
	//{
	//	fout << deltaT * i << "\n";
	//}
	//fout.close();
}
