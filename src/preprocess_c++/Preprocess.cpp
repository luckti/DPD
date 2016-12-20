/* -------------------------------------------------------------------
Preprocess.cpp
generate the configuration of wall, chain and free dpd parti-
cels for simulating polymer solutions. Configuration data can
be input from the terminal and wtitten into "fort.10" at the 
first time and then read from it next time.
-------------------------------------------------------------------*/

#include <direct.h>
#include "stdafx.h"

extern int inptconf(string);
extern int setParticles();
extern int optconf(string);

int main()
{
	
	//_mkdir("./data");
	inptconf("./data/lint.dat");	
	setParticles();
	optconf("./data/initconf");

	//system("pause");
	return 0;
}