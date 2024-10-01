#pragma once
#include <algorithm>

// функция формы
double functionForm (const double& localCoordinateKsi, const double& localCoordinateEtta,
					const double& localCoordinatePsi, const double& quadraticPointKsi,
					const double& quadraticPointEtta, const double& quadraticPointPsi)
{
	double N = 1.0 / 8.0 * (1.0 + localCoordinateKsi * quadraticPointKsi) * 
			   (1.0 + localCoordinateEtta * quadraticPointEtta) *
			   (1.0 + localCoordinatePsi * quadraticPointPsi);
	return N;
}

// производаня функции формы по ksi
double functionFormDerivativeKsi (const double& localCoordinateKsi, const double& localCoordinateEtta,
								 const double& localCoordinatePsi, const double& quadraticPointEtta,
								 const double& quadraticPointPsi, const double& length)
{
	double derivativeFormKsi = 1.0 / 4.0 * localCoordinateKsi / length * 
							   (1.0 + localCoordinateEtta * quadraticPointEtta) * 
							   (1.0 + localCoordinatePsi * quadraticPointPsi);
	return derivativeFormKsi;
}

// производаня функции формы по etta
double functionFormDerivativeEtta (const double& localCoordinateKsi, const double& localCoordinateEtta,
								  const double& localCoordinatePsi, const double& quadraticPointKsi,
								  const double& quadraticPointPsi, const double& width)
{
	double derivativeFormEtta = 1.0 / 4.0 * (1.0 + localCoordinateKsi * quadraticPointKsi) *
								localCoordinateEtta / width * (1.0 + localCoordinatePsi * quadraticPointPsi);
	return derivativeFormEtta;
}

// производаня функции формы по psi
double functionFormDerivativePsi (const double& localCoordinateKsi, const double& localCoordinateEtta,
								 const double& localCoordinatePsi, const double& quadraticPointKsi,
								 const double& quadraticPointEtta, const double& h)
{
	double derivativeFormPsi = 1.0 / 4.0 * (1.0 + localCoordinateKsi * quadraticPointKsi) *
							   (1.0 + localCoordinateEtta * quadraticPointEtta) * localCoordinatePsi / h;
	return derivativeFormPsi;
}

// Наложений условий симметрии для статической задачи
void boundConditionStatic (Matrix<double>& stiffness)
{
	stiffness[2][2] *= 10.0e10;
	stiffness[5][5] *= 10.0e10;
	stiffness[8][8] *= 10.0e10;
	stiffness[11][11] *= 10.0e10;

	stiffness[6][6] *= 10.0e10;
	stiffness[9][9] *= 10.0e10;
	stiffness[18][18] *= 10.0e10;
	stiffness[21][21] *= 10.0e10;

	stiffness[1][1] *= 10.0e10;
	stiffness[10][10] *= 10.0e10;
	stiffness[13][13] *= 10.0e10;
	stiffness[22][22] *= 10.0e10;
}

// Наложений условий симметрии для динамической заачи
void boundConditionsDinamic (Matrix<double>& matrixStiffness,
							Matrix<double>& matrixMass, Matrix<double>& displacement,
							std::vector<double>& speed,
							std::vector<double>& acceleration, std::vector<double>& force)
{
	std::vector<size_t> indexsBoundConditions = { 1, 2, 5, 6, 8, 9, 10 ,11, 13, 18, 21, 22 };
	size_t rowsStiffness = matrixStiffness.sizeRows();
	size_t displacementRows = displacement.sizeRows();

	size_t j = 0;
	for (size_t i = 0; i < rowsStiffness; ++i)
	{
		if (std::find(indexsBoundConditions.begin(), indexsBoundConditions.end(), i) != indexsBoundConditions.end())
		{
			matrixStiffness.eraseMatrix(j);
			matrixMass.eraseMatrix(j);
			speed.erase(speed.begin() + j);
			acceleration.erase(acceleration.begin() + j);
			force.erase(force.begin() + j);
		}
		else
			j++;
	}

	size_t StiffnessRowsNew = matrixStiffness.sizeRows();
	for (size_t i = 0; i < StiffnessRowsNew; ++i)
	{
		j = 0;
		size_t k = 0;
		while (j != StiffnessRowsNew)
		{
			if (std::find(indexsBoundConditions.begin(), indexsBoundConditions.end(), k) != indexsBoundConditions.end())
			{
				matrixStiffness[i].erase(matrixStiffness[i].begin() + j);
				matrixMass[i].erase(matrixMass[i].begin() + j);

				if (i < displacementRows)
					displacement[i].erase(displacement[i].begin() + j);
			}
			else
				j++;
			k++;
		}
	}
}

// Решатель для статической задачи
std::vector<double> displacementStatic(const Matrix<double>& stiffness,
	const std::vector<double>& force)
{
	const size_t size_MR = stiffness.sizeRows();
	std::vector<double> x(size_MR);
	Matrix<double> L = matrixL(stiffness);

	for (size_t i = 0; i < size_MR; ++i)
	{
		double sum = 0;
		for (size_t j = 0; j < i; ++j)
		{
			sum += L[i][j] * x[j];
		}
		x[i] = (force[i] - sum) / L[i][i];
	}

	std::vector<double> q(size_MR);

	for (size_t i = size_MR; i-- > 0;)
	{
		double sum = 0;
		for (size_t j = i + 1; j < size_MR; ++j)
		{
			sum += L[j][i] * q[j];
		}
		q[i] = (x[i] - sum) / L[i][i];
	}
	return q;
}

// Решатель для динамической задачи
Matrix<double> displacementDinamic (Matrix<double>& matrixStiffness, Matrix<double>& matrixMass,
								   const std::vector<double>& initialDisplacements,
								   const std::vector<double>& initialSpeed,
								   const std::vector<double>& initialAcceleration,
								   const double& Time,const double& deltaT)
{
	const size_t rowsStiffness = matrixStiffness.sizeRows();
	size_t displacementsRows = static_cast<size_t> (Time / deltaT);

	double alpha = 0.25;
	double delta = 0.5;

	Matrix<double> displacements(displacementsRows, rowsStiffness);
	displacements[0] = initialDisplacements;

	std::vector<double> speedOld = initialSpeed;
	std::vector<double> speedNew = {};

	std::vector<double> accelerationOld = initialAcceleration;
	std::vector<double> accelerationNew = {};

	std::vector<double> force(rowsStiffness);

	size_t n = 0;
	double t = 0.0;
	double sign = 0;

	std::cout << "Input 1 or 2, or 3, where 1 is oscillations without a driving force, 2 is dry friction and 3 is viscous friction";
	while (n != 1 and n != 2 and n != 3)
	{
		std::cin >> n;
	}

	if (n == 1)
	{
		std::vector<size_t> indexsNormalReaction = { 4, 7, 14, 16, 17, 19, 20, 23 };
		double coeffDryFrictionRest = 0.0;
		double coeffDryFrictionSliding = 0.0;
		double normalReac = 400.0;
		double averagePointsSpeedOld = 0.0;

		std::cout << "Input coefficient of dry friction at rest" << "\n";
		std::cin >> coeffDryFrictionRest;
		std::cout << "Input coefficient of dry friction at sliding" << "\n";
		std::cin >> coeffDryFrictionSliding;

		// Задание всестороннего обжатия
		for (auto iter = indexsNormalReaction.begin(); iter != indexsNormalReaction.end(); ++iter)
		{
			force[*iter] = -1.0 * normalReac;
		}

		boundConditionsDinamic (matrixStiffness, matrixMass, displacements, speedOld, accelerationOld, force);

		for (size_t i = 0; i < displacementsRows - 1; i++)
		{
			double elasticForce = -1.0 * matrixStiffness[0] * displacements[i] - matrixStiffness[1] * displacements[i] -
								  matrixStiffness[4] * displacements[i] - matrixStiffness[6] * displacements[i];
			double eps = 1.0e-4;
			double averagePointsSpeed = (speedOld[0] + speedOld[1] + speedOld[4] + speedOld[6]) / 4.0;
			double averagePointsAcceleration = (accelerationOld[0] + accelerationOld[1] + accelerationOld[4] +
												accelerationOld[6]) / 4.0;
			//std::cout << averagePointsAcceleration << "   " << averagePointsSpeed << "\n";

			if ((averagePointsSpeed > 0.0 and abs(averagePointsSpeed) > eps) or
				(abs(averagePointsSpeed) < eps and elasticForce > 0.0))
				sign = -1.0;
			else if ((averagePointsSpeed < 0.0 and abs(averagePointsSpeed) > eps) or
				(abs(averagePointsSpeed) < eps and elasticForce < 0.0))
				sign = 1.0;

			double coeffDryFriction = 0.0;
			
			if (abs(averagePointsSpeed) < eps)
				coeffDryFriction = coeffDryFrictionRest;
			else
				coeffDryFriction = coeffDryFrictionSliding;
			

			force[1] = sign * coeffDryFriction * normalReac;
			force[4] = sign * coeffDryFriction * normalReac;
			force[6] = 2.0 * sign * coeffDryFriction * normalReac;

			// ПРОВЕРКА НА ЗАЛИПАНИЯ - ПРОЕКЦИЯ ВСЕХ СИЛ НА ОСЬ KSI
			double frictionForce = 4.0 * coeffDryFriction * normalReac;

			if (abs(elasticForce) <= abs(frictionForce) and abs(averagePointsSpeed) < eps and
				abs(averagePointsSpeedOld) < eps)
			{
				displacements[i + 1] = displacements[i];
			}
			else
			{
				std::vector<double> b = alpha * pow(deltaT, 2) * force + matrixMass * displacements[i] + 
										deltaT * matrixMass * speedOld - (alpha - 0.5) * pow(deltaT, 2) * matrixMass * accelerationOld;

				Matrix<double> a = matrixMass + alpha * pow(deltaT, 2) * matrixStiffness;
				displacements[i + 1] = solverGauss(a, b);


				accelerationNew = (1.0 / (alpha * pow(deltaT, 2))) * (displacements[i + 1] - displacements[i]) -
								   (1.0 / (alpha * deltaT)) * speedOld + (1.0 - 1.0 / (2.0 * alpha)) * accelerationOld;

				speedNew = (delta / (alpha * deltaT)) * (displacements[i + 1] - displacements[i]) + 
						   (1.0 - delta / alpha) * speedOld + (1.0 - delta / (2.0 * alpha)) * deltaT * accelerationOld;

				speedOld = speedNew;
				accelerationOld = accelerationNew;

			}
			averagePointsSpeedOld = averagePointsSpeed;
		}
	}
	else if (n == 2)
	{
		std::vector<size_t> indexsNormalReaction = { 4, 7, 14, 16, 17, 19, 20, 23 };
		double coeffDryFrictionRest = 0.0;
		double coeffDryFrictionSliding = 0.0;
		double omega = acos(-1.0);
		double normalReac = 100.0;
		double averagePointsSpeedOld = 0.0;
		double averagePointsSpeedReal = 0.0;

		std::cout << "Input coefficient of dry friction at rest" << "\n";
		std::cin >> coeffDryFrictionRest;
		std::cout << "Input coefficient of dry friction at sliding" << "\n";
		std::cin >> coeffDryFrictionSliding;

		for (auto iter = indexsNormalReaction.begin(); iter != indexsNormalReaction.end(); ++iter)
		{
			force[*iter] = -normalReac;
		}

		boundConditionsDinamic (matrixStiffness, matrixMass, displacements, speedOld, accelerationOld, force);

		for (size_t i = 0; i < displacementsRows - 1; ++i)
		{
			double driveForce = 150.0 * cos(omega * t);
			double elasticForce = -1.0 * matrixStiffness[0] * displacements[i] - 
				matrixStiffness[1] * displacements[i] - matrixStiffness[4] * displacements[i] -
				matrixStiffness[6] * displacements[i];

			double eps = 1.0e-4;
			double averagePointsSpeed = (speedOld[0] + speedOld[1] + speedOld[4] + speedOld[6]) / 4;

			if ((averagePointsSpeed > 0.0 and abs(averagePointsSpeed) > eps) or
				(abs(averagePointsSpeed) < eps and (4 * driveForce + elasticForce > 0.0)))
				sign = -1.0;
			else if ((averagePointsSpeed < 0.0 and abs(averagePointsSpeed) > eps) or
					(abs(averagePointsSpeed) < eps and (4 * driveForce + elasticForce < 0.0)))
				sign = 1.0;

			double coeffDryFriction = 0.0;

			if (abs(averagePointsSpeed) < eps)
					coeffDryFriction = coeffDryFrictionRest;
				else
					coeffDryFriction = coeffDryFrictionSliding;


			// ПРОВЕРКА НА ЗАЛИПАНИЯ
			double frictionForce = 4.0 * sign * coeffDryFriction * normalReac;

			bool elasticCondition = ((abs(elasticForce) - abs(4 * driveForce + frictionForce) < -eps) and
				driveForce * sign >= 0.0) and abs(averagePointsSpeed) < eps and
				abs(averagePointsSpeedOld) < eps and abs(averagePointsSpeedReal) < eps;

			bool driveCondition = ((abs(4 * driveForce) - abs(elasticForce + frictionForce)) < -eps and
				elasticForce * sign >= 0.0) and abs(averagePointsSpeed) < eps and
				abs(averagePointsSpeedOld) < eps and abs(averagePointsSpeedReal) < eps;
				
			bool DrivePlusElasticCondition = ((abs(4 * driveForce + elasticForce) - abs(frictionForce) < -eps) and
				driveForce * sign <= 0.0 and sign * elasticForce <= 0.0) and abs(averagePointsSpeed) < eps and
				abs(averagePointsSpeedOld) < eps and abs(averagePointsSpeedReal) < eps;

			if (elasticCondition or driveCondition or DrivePlusElasticCondition)
			{
				displacements[i + 1] = displacements[i];
			}
			else 
			{
				force[0] = driveForce;
				force[1] = sign * coeffDryFriction * normalReac + driveForce;
				force[4] = sign * coeffDryFriction * normalReac + driveForce;
				force[6] = 2.0 * sign * coeffDryFriction * normalReac + driveForce;

				std::vector<double> b = alpha * pow(deltaT, 2) * force + matrixMass * displacements[i] +
										deltaT * matrixMass * speedOld -
										(alpha - 0.5) * pow(deltaT, 2) * matrixMass * accelerationOld;

				Matrix<double> a = matrixMass + alpha * pow(deltaT, 2) * matrixStiffness;
				displacements[i + 1] = solverGauss(a, b);


				accelerationNew = (1.0 / (alpha * pow(deltaT, 2))) * (displacements[i + 1] - displacements[i]) - 
								   (1.0 / (alpha * deltaT)) * speedOld + (1.0 - 1.0 / (2.0 * alpha)) * accelerationOld;

				speedNew = (delta / (alpha * deltaT)) * (displacements[i + 1] - displacements[i]) +
							(1.0 - delta / alpha) * speedOld + (1.0 - delta / (2.0 * alpha)) * deltaT * accelerationOld;

				speedOld = speedNew;
				accelerationOld = accelerationNew;
			}

			averagePointsSpeedOld = averagePointsSpeedReal;
			averagePointsSpeedReal = averagePointsSpeed;
			t += deltaT;
		}
	}
	else
	{
		double omega = 2.0 * acos(-1.0);
		double betta = 0.0;
		std::cout << "Input coefficient betta" << "\n";
		std::cin >> betta;

		boundConditionsDinamic (matrixStiffness, matrixMass, displacements, speedOld, accelerationOld, force);

		for (size_t i = 0; i < displacementsRows - 1; ++i)
		{
			force[0] = 55.0 * cos(omega * t);
			force[1] = 55.0 * cos(omega * t);
			force[4] = 55.0 * cos(omega * t);
			force[6] = 55.0 * cos(omega * t);

			std::vector<double> b = alpha * pow(deltaT, 2) * force + matrixMass * displacements[i] + 
									delta * deltaT * (betta * matrixMass) * displacements[i] + deltaT * matrixMass * speedOld + 
									(delta - alpha) * pow(deltaT, 2) * (betta * matrixMass) * speedOld - 
									(alpha - 0.5) * pow(deltaT, 2) * matrixMass * accelerationOld -
									pow(deltaT, 3) * (alpha - delta * 0.5) * (betta * matrixMass) * accelerationOld;

			Matrix<double> a = matrixMass + alpha * pow(deltaT, 2) * matrixStiffness + delta * deltaT * (betta * matrixMass);
			displacements[i + 1] = solverGauss(a, b);

			accelerationNew = (1.0 / (alpha * pow(deltaT, 2))) * (displacements[i + 1] - displacements[i]) -
							(1.0 / (alpha * deltaT)) * speedOld + (1.0 - 1.0 / (2.0 * alpha)) * accelerationOld;

			speedNew = (delta / (alpha * deltaT)) * (displacements[i + 1] - displacements[i]) +
					   (1.0 - delta / alpha) * speedOld + (1.0 - delta / (2.0 * alpha)) * deltaT * accelerationOld;

			speedOld = speedNew;
			accelerationOld = accelerationNew;
			t += deltaT;
		}
	}

	return displacements;
}
