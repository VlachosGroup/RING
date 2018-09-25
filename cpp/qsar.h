#ifndef QSAR_H
#define QSAR_H

#include "molecule.h"

class Surfqsar
{
	protected:
		double KHZeroCalc();
		Molecule* mol;
		int NumEOUnits;
		double AICSecond();
		double RelativeOxygensCount();
		double KierShapeThird();
		double ABInfoContentZero();
		double StructInfoContentOne();
		
	public:
		Surfqsar(Molecule*, int);
		double calculateCMC();
		double CMCerror();
		double calculateCP();
		double CPerror();
		double calculateHLB();
		double calculateSurfTens();
};

#endif
