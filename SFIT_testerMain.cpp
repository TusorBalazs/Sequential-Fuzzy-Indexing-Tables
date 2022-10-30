#include <stdlib.h>
#include <stdio.h>
#include "windows.h"
#include "SFIT.h"

using namespace std;

int main()
{
	SFIT classifier;
	
	classifier.train("data/wisctrain.dat", 0, 0.9);
	classifier.runDiagnostics();
	classifier.evaluateNoDisplay("data/wisctrain2.dat");
	classifier.evaluateNoDisplay("data/wisctest2.dat");

	cout << "\n Job's done!";
	char c;
	cin >> c;
	return 0;
}