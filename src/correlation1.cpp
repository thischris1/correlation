//============================================================================
// Name        : correlation1.cpp
// Author      : Christian Müller
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
// System includes
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdio.h>
#include <unistd.h>
// manyparticle includes

#include <LLLlib/LLLlib.h>
#include <LLLlib/yosEigenState.hpp>
#include <LLLlib/yosBasis.hpp>
#include <LLLlib/Types.hpp>
#include <geometry/CxPeriodicPosition.hpp>
#include <LLLlib/LLLlib.h>

#include <geometry/CxVortex.hpp>
#include <utils/CDistribution.hpp>
#include <utils/CDistanceDistribution.h>
#include <geometry/CxPosition.hpp>
#include <utils/CRandomizerGsl.h>

#include <utils/CFileParser.h>
#include <utils/CLine.h>
#include <utils/C2dDistribution.hpp>

#include <LLLlib/valForColumn.hpp>
#include "utilities.h"

using namespace std;

/*int main() {
 cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
 glLogger.error("Hello World 2");
 yosBasis * aBasis = new yosBasis(5,15,0);
 yosEigenState aState(aBasis, 0);
 return 0;
 }*/

/*
 void findZeros (yosBasis &myBasis,
 persist_eigSt eig,
 posVector elecPos,
 std::string zeroFilename="./zeros.dat")
 {

 try
 {
 yosEigenState searchState (&myBasis,eig,elecPos);
 searchState.findZeros();

 }
 catch( CxErrors &e)
 {
 e.print();
 exit (-1);
 }

 int zeroCount = 0;
 posVector zeroValues = searchState.getShallowZeros();
 if (zeroValues.size() <3) {
 glLogger.error("Did not find enough zeros");
 return;
 }
 std::vector <CxPeriodicPosition> zeros;
 vector <CxPeriodicPosition> electrons = searchState.getElectronPositions();
 glLogger.info("Number of zeros is (%i)",  zeroCount);

 std::ofstream of(zeroFilename.c_str());
 if (zeroValues.size() !=0)
 {
 of << "# Xposition        YPosition \n";
 for (unsigned int index = 0; index < electrons.size(); index++)
 {
 of << " #Electron "<< index << ": ";
 of <<  electrons[index].getXPosition() << "  " << electrons[index].getYPosition()<<endl;

 }

 for (unsigned int index = 0; index < zeroValues.size(); index++)
 {
 of << zeroValues[index].getXPosition() << "  " << zeroValues[index].getYPosition()<<" ";
 int nearPos =  zeroValues[index].getNearestPosition(electrons);
 of << electrons[nearPos].getXPosition();
 of << electrons[nearPos].getYPosition() << "\t";
 of << (electrons[nearPos]-zeroValues[index]).vabs()<<endl;

 }
 }

 for (unsigned int elecIndex = 0; elecIndex < electrons.size(); elecIndex ++)
 {
 double minimum = 1.0;
 int minIndex = 0;
 for ( int zeroIndex = 0; zeroIndex < zeroCount; zeroIndex++)
 {
 double tempDiff = (zeroValues[zeroIndex]-electrons[elecIndex]).vabs();

 if ( tempDiff< minimum )
 {
 if ( tempDiff  > 1e-06  )
 {
 minimum = tempDiff;;
 minIndex = zeroIndex;
 }
 }
 }

 (*outStream) << " "<< (electrons[elecIndex]-zeroValues[minIndex]).vabs();
 }

 }

 */
/*!


 \brief Helper function shows the usage of the program

 */
static void usage() {
	cerr << "Bad parameter usage, exiting\n";
	cerr << "Possible (and mandatory) options are \n \t -s [stateFileName]\n";
	cerr << "\t -b [baseFileName] \n";
	cerr
			<< "optional parameters \n\t-d [debugLevel] \n\t -f logFileName -n number of steps\n";
	exit(-1);
}
/*
 !
 \brief parses the command line and calls evaluate
 */
int main(int argc, char ** argv) {
	std::string positionFile, statefile, baseFile;

	const char* myOpt = "v:f:s:p:b:e:n:";
	int retVal = '?', monteCarloSteps = 500;
	bool hasStateFile = false;
	bool hasBasisFile = false;
	bool eeCorrelationOnly = false;
	bool isBinaryFile = true;
	ERROR("");
	glLogger.setLogLevel(ERROR);
	opterr = 0;
	glLogger.setFileName("correlationLog.txt");
	while ((retVal = getopt(argc, argv, myOpt)) != EOF) {
		switch (retVal) {
		case 's':
			statefile = optarg;
			hasStateFile = true;
			break;
		case 'b':
			baseFile = optarg;
			hasBasisFile = true;
			break;
		case 'v':
			cerr << "Loglevel = " << optarg << std::endl;
			glLogger.setLogLevel(std::string(optarg));
			glLogger.disableConsole();
			break;
		case 'f':
			glLogger.setFileName(optarg);
			break;
		case 'e':
			eeCorrelationOnly = true;
			break;
		case 'n':
			monteCarloSteps = atoi(optarg);
			break;
		case 'a':
			isBinaryFile = false;
			break;
		case '?':
			usage();
			break;
		}
	}

	if (!(hasBasisFile && hasStateFile)) {
		usage();
	}

	cerr << "Statefile = " << statefile << " positionFile =" << positionFile;
	cerr << "Basefile =" << baseFile << endl;
	cerr << "Entering zero search program \n";
	glLogger.error("Number of MC steps = %d", monteCarloSteps);
	yosBasis *basis = 0;
	yosEigenState *yosEigSt = 0;

	try {

		readBasisAndState(statefile, baseFile, basis, yosEigSt,
				monteCarloSteps, eeCorrelationOnly, isBinaryFile);
	} catch (CxErrors &e) {
		cerr << "An error was thrown \n";
		e.print();
		glLogger.error("Error %s was thrown", e.getMessage());
		exit(-1);
	}

	return (0);

}

