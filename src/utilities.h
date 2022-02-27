#ifndef UTILITIES_H
#define UTILITIES_H
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


//void findZeros (yosBasis &myBasis, persist_eigSt eig,  std::vector<CxPeriodicPosition> elecPos, std::string zeroFilename="./zeros.dat");
// std::vector<CxPeriodicPosition> readElectronVector(std::string fileName);


bool circleDensity(float radius, float radEpsilon, float angleEpsilon, int nr_columns, FILE * outf_dens, valForColumn * pValForColumn[]);


bool ringDensity(float radiusMax, float radiusMin, float radEpsilon, float angleEpsilon, int nr_columns, FILE * outf_dens, valForColumn * pValForColumn[]);

bool readBasisAndState( 
			const std::string & stateFileName,
		    const std::string & basisFileName,
	    	yosBasis *myBasis,
			yosEigenState *yosEig,
			int monteCarloSteps,
			bool eeCorrelationOnly = false, bool binaryFile = true);
	
bool monteCarlo(yosEigenState *m_state);



bool readSpectrumFromState( const std::string & stateFileName,
			    const std::string & basisFileName,
			    bool isBinary);


yosBasis readBasis(const std::string & basisFileName,   bool isBinary);


float getMinimalDistance (  std::vector<CxPeriodicPosition> & electrons, 
			    std::vector<CxPeriodicPosition> & vortices);


int removeAllPositionsWithinX(  std::vector<CxPeriodicPosition> & electrons, 
				std::vector<CxPeriodicPosition> & vortices,
				const double & distance);
  
int removeMostDistantVortices ( std::vector<CxPeriodicPosition> & electrons, 
				std::vector<CxPeriodicPosition> & vortices,
				const int & numOfVor);


float getMinimalDistance (  std::vector<CxPeriodicPosition> & electrons, 
			    std::vector<CxVortex> & vortices);


int removeAllPositionsWithinX(  std::vector<CxPeriodicPosition> & electrons, 
				std::vector<CxVortex> & vortices,
				const double & distance);
  
int removeMostDistantVortices ( std::vector<CxPeriodicPosition> & electrons, 
				std::vector<CxVortex> & vortices,
				const int & numOfVor);


void debugVortices(std::vector<CxVortex> & vortices, const std::string & message = "");

void debugPositions(std::vector<CxPeriodicPosition> & vortices, const std::string & message = "") ;

void writeSpectrum (std::vector<double> spectrum, std::string fileName);
#endif
