#pragma once

// Матрица упругих постоянных
Matrix<double> matrixElasticConstants(const size_t& rows, const size_t& columns,
	const double& modulusElastic, const double& coefficientPuasson)
{
	Matrix<double> matrix(rows, columns);
	for (size_t i = 0; i < rows; ++i)
	{
		for (size_t j = 0; j < columns; ++j)
		{
			if (i < 3)
			{
				if (i == j)
					matrix[i][j] = modulusElastic * (1.0 - coefficientPuasson) /
					((1.0 + coefficientPuasson) * (1.0 - 2.0 * coefficientPuasson));
				else if (j < 3)
					matrix[i][j] = modulusElastic * coefficientPuasson /
					((1.0 + coefficientPuasson) * (1.0 - 2.0 * coefficientPuasson));
			}
			else
			{
				if (i == j)
					matrix[i][j] = modulusElastic / (2.0 * (1.0 + coefficientPuasson));
			}
		}
	}
	return matrix;
}

Matrix<double> quadraticPoints(const size_t& rows, const size_t& columns, const Matrix<double>& localCoordinate)
{
	Matrix<double> matrix(rows, columns);
	for (size_t i = 0; i < rows; ++i)
	{
		for (size_t j = 0; j < columns; ++j)
		{
			matrix[i][j] = sqrt(3.0) / 3.0 * localCoordinate[i][j];
		}
	}
	return matrix;
}

// Транспонированная матрица дифференцирования на матрицу упругих постоянных
Matrix<double> BtDMatrix(const Matrix<double>& Btransp, const Matrix<double>& D)
{
	const size_t size_c_Btransp = Btransp.sizeColumns();
	const size_t size_c_D = D.sizeColumns();
	Matrix<double> result(size_c_Btransp, size_c_D);

	for (size_t i = 0; i < size_c_Btransp; ++i)
	{
		for (size_t j = 0; j < size_c_D; ++j)
		{
			for (size_t k = 0; k < size_c_D; ++k)
			{
				result[i][j] += Btransp[k][i] * D[k][j];
			}
		}
	}
	return result;
}

// Диагональная матрица масс
Matrix<double> massDiagonal(const size_t& sizeMatrixMass, const double& dencity,
	const double& length, const double& width, const double& heigth)
{
	Matrix<double> matrixMass(sizeMatrixMass, sizeMatrixMass);
	size_t quantitiyNode = sizeMatrixMass / 3;
	double mass = dencity * length * width * heigth;
	for (size_t i = 0; i < sizeMatrixMass; ++i)
	{
		matrixMass[i][i] = mass / quantitiyNode;
	}
	return matrixMass;
}

//// Совместная матрица масс
Matrix<double> massJoint(const size_t& sizeMatrixMass, const double& dencity, const double& length,
						 const double& width, const double& heigth)
{
	const double determinantMatrixJacobian = length * width * heigth / 8.0;

	//локальные координаты
	Matrix<double> localCoordinates = { {1.0 ,1.0 ,-1.0, -1.0, 1.0, 1.0, -1.0, -1.0 },
		{ -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0 }, { -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0 } };

	const size_t m = localCoordinates.sizeRows();

	// Матрица квадратичных точек
	Matrix<double> matrixQuadraticPoints = quadraticPoints(m, 8, localCoordinates);

	const size_t n = matrixQuadraticPoints.sizeColumns();
	Matrix<double> matrixMass(sizeMatrixMass, sizeMatrixMass);


	for (size_t i = 0; i < n; ++i)
	{
		for (size_t j = 0; j < n; ++j)
		{
			for (size_t h = 0; h < n; ++h)
			{
				matrixMass[3 * i][3 * j] += functionForm(localCoordinates[0][i], localCoordinates[1][i],
					localCoordinates[2][i], matrixQuadraticPoints[0][h],
					matrixQuadraticPoints[1][h], matrixQuadraticPoints[2][h]) *
					functionForm(localCoordinates[0][j], localCoordinates[1][j],
						localCoordinates[2][j], matrixQuadraticPoints[0][h],
						matrixQuadraticPoints[1][h], matrixQuadraticPoints[2][h]);

				matrixMass[3 * i + 1][3 * j + 1] += functionForm(localCoordinates[0][i], localCoordinates[1][i],
					localCoordinates[2][i], matrixQuadraticPoints[0][h],
					matrixQuadraticPoints[1][h], matrixQuadraticPoints[2][h]) *
					functionForm(localCoordinates[0][j], localCoordinates[1][j],
						localCoordinates[2][j], matrixQuadraticPoints[0][h],
						matrixQuadraticPoints[1][h], matrixQuadraticPoints[2][h]);

				matrixMass[3 * i + 2][3 * j + 2] += functionForm(localCoordinates[0][i], localCoordinates[1][i],
					localCoordinates[2][i], matrixQuadraticPoints[0][h],
					matrixQuadraticPoints[1][h], matrixQuadraticPoints[2][h]) *
					functionForm(localCoordinates[0][j], localCoordinates[1][j],
						localCoordinates[2][j], matrixQuadraticPoints[0][h],
						matrixQuadraticPoints[1][h], matrixQuadraticPoints[2][h]);
			}

			matrixMass[3 * i][3 * j] *= determinantMatrixJacobian * dencity;

			matrixMass[3 * i + 1][3 * j + 1] *= determinantMatrixJacobian * dencity;

			matrixMass[3 * i + 2][3 * j + 2] *= determinantMatrixJacobian * dencity;
		}
	}
	return matrixMass;
}


Matrix<double> matrixStiffness(const double& length, const double& width, const double& heigth,
	const double& modulusEl, const double& coeffPuas)
{
	const double determinantMatrixJacobian = length * width * heigth / 8.0;

	Matrix<double> localCoordinates = { {1.0 ,1.0 ,-1.0, -1.0, 1.0, 1.0, -1.0, -1.0 },
		{ -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0 }, { -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0 } };

	const size_t n = localCoordinates.sizeRows();

	// Матрица упругих постоянных
	Matrix<double> matrixElasticConst = matrixElasticConstants(6, 6, modulusEl, coeffPuas);
	const size_t m = matrixElasticConst.sizeRows();

	// Матрица квадратичных точек
	Matrix<double> matrixQuadraticPoints = quadraticPoints(n, 8, localCoordinates);
	const size_t h = matrixQuadraticPoints.sizeColumns();

	// Матрица жесткости
	const size_t rowsMatrixStiffnes = n * h;
	Matrix<double> matrixStiff(rowsMatrixStiffnes, rowsMatrixStiffnes);

	for (size_t l = 0; l < h; l++)
	{
		Matrix<double> matrixDifferentiation(m, rowsMatrixStiffnes);
		size_t g = 0;
		for (size_t k = 0; k < rowsMatrixStiffnes - 2; k += 3)
		{
			matrixDifferentiation[0][k] = functionFormDerivativeKsi(localCoordinates[0][g], localCoordinates[1][g],
				localCoordinates[2][g], matrixQuadraticPoints[1][l], matrixQuadraticPoints[2][l], length);

			matrixDifferentiation[1][k + 1] = functionFormDerivativeEtta(localCoordinates[0][g], localCoordinates[1][g],
				localCoordinates[2][g], matrixQuadraticPoints[0][l], matrixQuadraticPoints[2][l], width);

			matrixDifferentiation[2][k + 2] = functionFormDerivativePsi(localCoordinates[0][g], localCoordinates[1][g],
				localCoordinates[2][g], matrixQuadraticPoints[0][l], matrixQuadraticPoints[1][l], heigth);
			g++;
		}

		g = 0;
		for (size_t k = 0; k < rowsMatrixStiffnes - 2; k += 3)
		{
			matrixDifferentiation[3][k] = functionFormDerivativeEtta(localCoordinates[0][g], localCoordinates[1][g],
				localCoordinates[2][g], matrixQuadraticPoints[0][l], matrixQuadraticPoints[2][l], width);

			matrixDifferentiation[3][k + 1] = functionFormDerivativeKsi(localCoordinates[0][g], localCoordinates[1][g],
				localCoordinates[2][g], matrixQuadraticPoints[1][l], matrixQuadraticPoints[2][l], length);
			g++;
		}

		g = 0;
		for (size_t k = 0; k < rowsMatrixStiffnes - 2; k += 3)
		{
			matrixDifferentiation[4][k + 1] = functionFormDerivativePsi(localCoordinates[0][g], localCoordinates[1][g],
				localCoordinates[2][g], matrixQuadraticPoints[0][l], matrixQuadraticPoints[1][l], heigth);

			matrixDifferentiation[4][k + 2] = functionFormDerivativeEtta(localCoordinates[0][g], localCoordinates[1][g],
				localCoordinates[2][g], matrixQuadraticPoints[0][l], matrixQuadraticPoints[2][l], width);
			g++;
		}

		g = 0;
		for (size_t k = 0; k < rowsMatrixStiffnes - 2; k += 3)
		{
			matrixDifferentiation[5][k] = functionFormDerivativePsi(localCoordinates[0][g], localCoordinates[1][g],
				localCoordinates[2][g], matrixQuadraticPoints[0][l], matrixQuadraticPoints[1][l], heigth);
			matrixDifferentiation[5][k + 2] = functionFormDerivativeKsi(localCoordinates[0][g], localCoordinates[1][g],
				localCoordinates[2][g], matrixQuadraticPoints[1][l], matrixQuadraticPoints[2][l], length);
			g++;
		}

		Matrix<double> transposeMatrixDiffbyMatrixEC = BtDMatrix(matrixDifferentiation, matrixElasticConst);

		for (size_t i = 0; i < rowsMatrixStiffnes; ++i)
		{
			for (size_t j = 0; j < rowsMatrixStiffnes; ++j)
			{
				double product = 0;
				for (size_t k = 0; k < m; ++k)
				{
					product += transposeMatrixDiffbyMatrixEC[i][k] * matrixDifferentiation[k][j];
				}
				matrixStiff[i][j] += product * determinantMatrixJacobian;
			}
		}
	}

	return matrixStiff;
}