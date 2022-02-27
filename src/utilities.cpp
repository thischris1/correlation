#include "utilities.h"
#include <LLLlib/LLLlib.h>
#include <LLLlib/yosEigenState.hpp>
#include <LLLlib/yosBasis.hpp>
#include <LLLlib/Types.hpp>
#include <LLLlib/persist_eigSt.hpp>
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

// System includes
// boost
#include <boost/filesystem.hpp>
using namespace boost::filesystem;
using namespace std;
void debugVortices(std::vector<CxVortex> & vortices, const std::string & message )
{
	if (glLogger.getLogLevel() != DEBUG)
	{
		return;
	}
	glLogger.debug(message.c_str());
	std::vector<CxVortex>::iterator debugIt = vortices.begin();

	while (debugIt != vortices.end())
	{
		glLogger.debug("Vortex position (%f), (%f), cell X (%f), Y = (%f)",
				debugIt->getXPosition(), debugIt->getYPosition(),
				debugIt->get_xCellSize(), debugIt->get_yCellSize());
		debugIt++;
	}
	
}

void debugPositions(std::vector<CxPeriodicPosition> & positions, const std::string & message )
{
	
	if (glLogger.getLogLevel() != DEBUG)
	{
		return;
	}
	glLogger.debug(message.c_str());
	std::vector<CxPeriodicPosition>::iterator debugIt = positions.begin();
	
	while (debugIt != positions.end())
		{
			glLogger.debug("Vortex position (%f), (%f)",
					debugIt->getXPosition(), debugIt->getYPosition());
			debugIt++;
		}
	
}


yosBasis readBasis(const std::string & basisFileName,   bool isBinary)
{

	if (basisFileName.size() == 0)
	{

		throw (new CxFileNotFoundError("basis-state file not found", __FILE__,__LINE__));
	}



	/*
    Start constructing basis, reading state and positions
	 */

	glLogger.info("reading basis");
	FILE *basisFile = fopen(basisFileName.c_str(),"rb");

	if (!basisFile)
	{
		glLogger.error("Could not open basis file");
		throw (new CxFileNotFoundError("Could not open basis file", __FILE__,__LINE__));
	}
	yosBasis myBasis(basisFile, 0);
	fclose (basisFile);
	return  (myBasis);

}

/*!
  \brief starts the evaluation part by reading basis and states from files
  \param stateFileName name of the stateFile
  \param basisFileName name f the basis file 
  \param positionFileName name of position file
 */

bool circleDensity(float radius, float radEpsilon, float angleEpsilon, 
		int nr_columns, FILE * outf_dens, valForColumn * pValForColumn[])
{ 
	glLogger.debug("Start circleDensity with radisu (%f), radiusStep (%f), angleStep(%f)", radius, radEpsilon, angleEpsilon);
	double charge = 0.0;
	bool firstStep = true;
	float angleStep = 2*M_PI/angleEpsilon;
	for ( float cRadius = 0.0f; cRadius <= radius; cRadius  = cRadius +radEpsilon)
	{

		double areaElement = 2*M_PI*radEpsilon*(radEpsilon+cRadius)/ angleEpsilon;
		glLogger.debug("Radius now (%f)", cRadius);
		if (!firstStep)
		{
			for ( float angle = 0; angle <=2*M_PI; angle = angle+angleStep)
			{

				float x = (0.5 + (cRadius*cos(angle)));
				float y = (0.5 + (cRadius*sin(angle)));
				fprintf (outf_dens, "%f %f %f ", x, y, areaElement);
				for (int k=0; k<nr_columns; k++)
				{
					double xtmp = (pValForColumn[k])->getValForColumn(x,y);
					fprintf (outf_dens, " %f \n", xtmp);
					// calculating "center-of-mass" for this operator

					if (k == 0) 
					{

						// sum up density
						charge = charge+ areaElement*xtmp;
					}
				}
			}
		}
		else 
		{ // first step, area is a circle

			float area = M_PI * radEpsilon*radEpsilon;
			double xtmp = (pValForColumn[0])->getValForColumn(0.5f,0.5f);
			fprintf (outf_dens, " %f \t %f \t %f \t %f\n ", 0.5f, 0.5f, xtmp,area);
			charge = charge+ area*xtmp;
			firstStep = false;
		}
	}
	fprintf (outf_dens, "# total charge %f \n", charge);   
	return (true);
}


bool ringDensity(float radiusMax, float radiusMin, float radEpsilon, 
		float angleEpsilon, int nr_columns, FILE * outf_dens, valForColumn * pValForColumn[])

{
	double charge = 0.0;
	for ( float cRadius = radiusMin; cRadius <= radiusMax; cRadius  = cRadius +radEpsilon)
	{
		double areaElement = cRadius*radEpsilon*angleEpsilon*0.5;
		for ( float angle = 0; angle <=2*M_PI; angle = angle+angleEpsilon)
		{

			float x = (0.5 + (cRadius*cos(angle)));
			float y = (0.5 + (cRadius*sin(angle)));
			fprintf (outf_dens, "%f %f", cRadius, angle);
			for (int k=0; k<nr_columns; k++)
			{
				double xtmp = (pValForColumn[k])->getValForColumn(x,y);
				fprintf (outf_dens, " %f", xtmp);
				if (k == 0) {

					// sum up density
					charge = charge+ areaElement*xtmp;
				}


				// calculating "center-of-mass" for this operator
			}
		}
		fprintf (outf_dens, "\n");
		fprintf (outf_dens, "# Total charge in sampled area (%f)", charge);

	}
	return (true);
}
/*!
  \brief function performs MC runs for the given eiegenstate in order to calculate distribution functions.
  \param m_state ptr to mp-state to be studied
  \return true if everything went fine
 */


/*!
  \brief starts the evaulation part by reading basis and states from files
  \param stateFileName name of the stateFile
  \param basisFileName name f the basis file 
  \param positionFileName name of position file
 */

bool readBasisAndState( const std::string  &stateFileName,
		const std::string & basisFileName,
		yosBasis *myBasis,
		yosEigenState *yosEig,
		int monteCarloSteps,
		bool eeCorrelationOnly,
		bool binaryFile)
{
	glLogger.info("Entering evaluate");
	boost::filesystem::path full_path( boost::filesystem::current_path() );
	std::cout << "Current path is : " << full_path << std::endl;

	std::string sStatefile = stateFileName;
	std::cerr << "HERE 0 \n";
	ERROR("HERE0");
	if (yosEig == 0)
	{
		std::cerr << ("HEREyosEig Null");

	}
	if (stateFileName.empty())
	{
		throw (new CxFileNotFoundError("mp-state filename empty", __FILE__,__LINE__));


	}
	if (basisFileName.empty())
	{

		throw (new CxFileNotFoundError("basis-state file not found", __FILE__,__LINE__));
	}



	/*
    Start constructing basis, reading state and positions
	 */

	glLogger.info("reading basis");
	path basisLocalPath (basisFileName);
	path basisPath = full_path /= basisFileName;
	std::cerr <<("HERE1");
	if (!exists(basisPath))
	{
		glLogger.error("Basis File does not exist");
		std::cout << "Current path is : " << basisPath << std::endl;
		glLogger.error(basisFileName.c_str());
		return false;
	}
	FILE *basisFile = fopen(basisFileName.c_str(),"rb");

	if (!basisFile)
	{
		glLogger.error("Could not open basis file");
		throw (new CxFileNotFoundError("Could not open basis file", __FILE__,__LINE__));
	}
	myBasis = new yosBasis (basisFile, 0);
	fclose (basisFile);

	persist_eigSt eig (myBasis, true);
	persist_eigSt tempEig (myBasis,true);
	std::ifstream fsVec;
	std::vector <double> spectrum;
	//! \todo change that!
	std::cerr << "Before Reading state \n";

	int iSt = 1;

	DEBUG ("Reading binary vector-file...\n");
	eig.set_binary(binaryFile);

	if (binaryFile)
	{
		fsVec.open (stateFileName.c_str(), std::ios::binary | std::ios::in);
	}
	else
	{
		fsVec.open(stateFileName.c_str(),std::ios::in);

	}
	if ( !fsVec )
	{
		fsVec.close();
		throw (new CxFileNotFoundError(stateFileName.c_str(),__FILE__,__LINE__ ));

	}

	int i = 0;
	while (!fsVec.eof())
	{
		if (i == 0)
		{
			fsVec>> (eig);
			spectrum.push_back(eig.getEn());
			break;
		}
		else {
			fsVec >> (tempEig);
			spectrum.push_back(tempEig.getEn());
		}

		i++;


	}
	std::string spectrumFileName("spectrum_");
//	spectrumFileName = spectrumFileName+stateFileName;
//	writeSpectrum(spectrum,spectrumFileName);
/*
	for(int i=0; i<iSt && !fsVec.eof(); i++) 
	{
		fsVec >> (eig);
	}
*/
	if (fsVec.eof())
	{
		fsVec.close();
		throw (new CxErrors("End of file reached, state not found",
				__FILE__,__LINE__));

	}
	fsVec.close();


	glLogger.info ("Imported state with energy %f (norm^2=%f)", eig.getEn(), eig.normsq()); 

	/*

  Search for fixed electron positions in a file "sampleElectrons.dat" in current working directory
	 */
	CFileParser * myParser = CFileParser::getInstance();
	std::vector<CxPosition> fixedPos;
	bool hasFixedPositions = false; //! True if fixed positions exist
	int numFixedPositions = 0;
	if (myParser->readFile("./sampleElectrons.dat"))
	{
		glLogger.debug("Found a file with electron fixedPosition");
		while (myParser->hasNextLine())
		{
			CLine tempLine = myParser->getNextLine();
			std::cerr << "Line =[" << tempLine.getLine(true)<<"]\n";

			CxPosition tempPos(tempLine.getFloats().at(0), 
					tempLine.getFloats().at(1));
			fixedPos.push_back(tempPos);
			hasFixedPositions = true;
			numFixedPositions++;
			glLogger.debug("added a position");
		}

	}
	glLogger.info("Found %d fixedPosition", numFixedPositions);

	/*

  Setup distributions, initialize 
	 */
	std::cerr << "Done with reading \n";

	CDistribution elecVortexD(1000);
	CDistribution vortexVortexD(1000);
	CDistribution elecElecD(1000);
	C2dDistribution elecVortex2d(200,200);
	C2dDistribution vortex2d(100,100);
	CDistribution * fixedD = new CDistribution[numFixedPositions];
	C2dDistribution *fixed2d = new C2dDistribution[numFixedPositions];
	CDistribution * fixedRaw = new   CDistribution[numFixedPositions];
	CDistribution * fixedwoPauli = new  CDistribution[numFixedPositions];
	for (int index = 0; index < numFixedPositions; index++)
	{
		fixedD[index] =  CDistribution(1000);
		fixed2d[index] = C2dDistribution(200,200);
		fixedRaw[index] = CDistribution(1000);
		fixedwoPauli[index] = CDistribution(1000);
	}

	std::vector <CxPosition> electrons(myBasis->getNe()-1- numFixedPositions);
	CRandomizer *myRand = new CRandomizerGsl();
	/*
    MC Loop 
	 */
	 bool loopControl = true;
	 int mcCounter = 0;
	 int badStepsCounter =0;
	 float oldWeight = 0.0f;
	 CxPeriodicPosition cellCenter(0.5, 0.5, 1.0);
	 int rejectCounter = 0;
	 while (loopControl)
	 {
		 std::cerr << "Starting loop \n";
		 glLogger.error("start While looP with index (%d)", mcCounter);
		 if ( myRand->randomPositions(electrons) != (myBasis->getNe()-1 - numFixedPositions))
		 {
			 return (false);
		 }

		 std::vector<CxPeriodicPosition> fullElectrons;
		 std::vector<CxPosition>::iterator elecIt = electrons.begin();
		 while (elecIt != electrons.end())
		 {
			 fullElectrons.push_back(CxPeriodicPosition(*elecIt, 1.0));
			 elecIt++;
		 }

		 if (hasFixedPositions)
		 {
			 std::vector <CxPosition>::iterator  vecIter = fixedPos.begin();

			 while (vecIter != fixedPos.end())
			 {
				 fullElectrons.push_back(*vecIter);
				 vecIter++;
			 }
		 }
		 if (glLogger.getLogLevel()== DEBUG)
		 {
			 debugPositions(fullElectrons);
			 
		 }
		 for (unsigned int index = 0; index < fullElectrons.size(); index ++)
		 {
			 glLogger.info("Utilies receives pos (%f), (%f) index (%d)", 
					 fullElectrons[index].getXPosition(),
					 fullElectrons[index].getYPosition(),
					 index );


		 }
		 std::cerr << "create eigenstate \n";
		 yosEigenState m_state(myBasis, eig, fullElectrons);
		 std::cerr << "have eigenstate \n";
		 m_state.initializeSearch();
		 /*
	Check wether step is acceptable
		  */
		 if (mcCounter > 0)  // Not at first step
		 {
			 float newWeight = m_state.getWeight();
			 float random = myRand->fRand();

			 glLogger.info("Old weight = (%f), newWeight = (%f), random =(%f)", oldWeight, newWeight,random);
			 if (newWeight/oldWeight < random )
			 {	  
				 // Not accepted step
				 rejectCounter ++;
				 glLogger.info("Rejected step");
				 continue;

			 }
			 oldWeight = newWeight;	  

		 }

		 glLogger.info("Accepted step, search zeros");
		 std::vector<CxPeriodicPosition> zeroPos;
		 std::vector<CxVortex> vortexPos;
		 // Valid step, now find zeros
		 if (eeCorrelationOnly == false)
		 {
			 m_state.findZeros();
			 zeroPos =  m_state.getShallowZeros();
			 vortexPos = m_state.getShallowVortices();
			 glLogger.error("Back from search zeros with %d zeros", vortexPos.size()); 
			 debugVortices(vortexPos);
			
			 debugPositions(fullElectrons,"All Electrons");
			 
			 double pauliDistance = getMinimalDistance(fullElectrons, vortexPos);
			 if (fabs(pauliDistance) < 5e-4) 
			 {
				 pauliDistance = 5e-4;
			 }

			 if (hasFixedPositions)
			 {
				 glLogger.info("Adding values for fixed electron");
				 for (unsigned int fixIndex = 0; fixIndex < fixedPos.size(); fixIndex++)
				 {
					 glLogger.info("Fixed Electron with index (%d), at (%f), (%f)",
							 fixIndex, fixedPos.at(fixIndex).getXPosition(),
							 fixedPos.at(fixIndex).getYPosition());
					 CxPeriodicPosition elecPos (fixedPos.at(fixIndex), 1.0f);
					 // Loop over all vortices
					 std::vector<CxVortex>::iterator vortexIt = vortexPos.begin();
					 while (vortexIt !=vortexPos.end())
					 {
						 // Calculate distance add to distribution TODO!
						 glLogger.info("Adding vortex Position x +(%f), y=(%f)",
								 vortexIt->getXPosition(), vortexIt->getYPosition());
						 fixedRaw[fixIndex].addValue((elecPos - (*vortexIt)).vabs());
						 
						 vortexIt++;
					 }

				 }
				 glLogger.info("Finished fixed electron");

			 }


			 removeAllPositionsWithinX(fullElectrons, vortexPos, pauliDistance);
			 glLogger.error("After removing Pauli zeros we are left with  %d", 
					 vortexPos.size()); 
			 
			 if (hasFixedPositions)
			 {
				 glLogger.info("Adding values for fixed electron");
				 for (unsigned int fixIndex = 0; fixIndex < fixedPos.size(); fixIndex++)
				 {
					 glLogger.info("Fixed Electron with index (%d), at (%f), (%f)",
							 fixIndex, fixedPos.at(fixIndex).getXPosition(),
							 fixedPos.at(fixIndex).getYPosition());
					 CxPeriodicPosition elecPos (fixedPos.at(fixIndex), 1.0f);
					 // Loop over all vortices
					 std::vector<CxVortex>::iterator vortexIt = vortexPos.begin();
					 while (vortexIt !=vortexPos.end())
					 {
						 // Calculate distance add to distribution TODO!
						 glLogger.info("Adding vortex Position x +(%f), y=(%f)",
								 vortexIt->getXPosition(), vortexIt->getYPosition());
						 fixedwoPauli[fixIndex].addValue((elecPos - (*vortexIt)).vabs());

						 vortexIt++;
					 }

				 }
				 glLogger.info("Finished fixed electron");

			 }

			 /* 
	     Now remove the center of mass vortices

			  */
			 /*
			  * Remove it for now
			  */
			 /*
	    removeMostDistantVortices(fullElectrons, vortexPos, m_state.getNe()-1);
	    debugIt = vortexPos.begin();
	    while (debugIt != vortexPos.end())
	    {
	    glLogger.info("remaining Vortex position (%f), (%f), cell X (%f), Y = (%f)",
	    debugIt->getXPosition(), debugIt->getYPosition(),
	    debugIt->get_xCellSize(), debugIt->get_yCellSize());
	    debugIt++;
	    }

	    glLogger.error("After removing center of mass vortices (%d) vortices are left", vortexPos.size());
			  */
			 if ( vortexPos.size() < m_state.getElectronPositions().size())
			 {
				 /*
		Too few zeros found, discarding step
		Count illegal steps
				  */
				 glLogger.info("Too few zeros found, discarding step");
				 badStepsCounter++;
				 continue;
			 }
		 }
		 mcCounter++;    
		 /*

      Calculate distributions now
		  */
		 if (eeCorrelationOnly == false)
		 {

			 for (unsigned int vortex1Index = 0; vortex1Index < vortexPos.size(); vortex1Index++)
			 {

				 CxPeriodicPosition vortex1( vortexPos[vortex1Index]);
				 vortex2d.addValue(vortex1);
				 glLogger.info("Outer vortex loop with index (%d) and positions (%f) (%f)", 
						 vortex1Index,
						 vortex1.getXPosition(),
						 vortex1.getYPosition());
				 for (unsigned int elecIndex = 0; elecIndex < electrons.size(); elecIndex++)
				 {
					 glLogger.info("Electron pos loop, electron No. (%d) at (%f), (%f)",
							 elecIndex,
							 electrons[elecIndex].getXPosition(),
							 electrons[elecIndex].getYPosition());
					 CxPeriodicPosition elecPos (electrons[elecIndex], 1.0f);
					 CxPeriodicPosition differenz = elecPos - vortex1;

					 differenz = differenz + CxPeriodicPosition(0.5,0.5);
					 elecVortex2d.addValue(differenz);
					 glLogger.info(" Add Position (%f), (%f) to 2 d distro", differenz.getXPosition(), differenz.getYPosition()); 
					 elecVortexD.addValue((elecPos -vortex1).vabs());
					 for (unsigned int elecIndex2 = elecIndex+1; elecIndex2 < electrons.size(); elecIndex2++)
					 {
						 CxPeriodicPosition elecPos2(electrons[elecIndex2]);

						 elecElecD.addValue((elecPos2-elecPos).vabs());
					 }
				 }
				 for (unsigned int vortex2Index = vortex1Index+1; 
				 vortex2Index  < vortexPos.size(); 
				 vortex2Index++)
				 {
					 glLogger.info("Inner vortex loop with index (%d) and positions (%f) (%f)", 
							 vortex2Index,
							 vortexPos[vortex2Index].getXPosition(),
							 vortexPos[vortex2Index].getYPosition());
					 CxPeriodicPosition vortex2( vortexPos[vortex2Index]);    
					 //	      vortexVortexD.addValue( (vortex1- vortexPos[vortex2Index]).vabs());
					 vortexVortexD.addValue((vortex1-vortex2).vabs());
				 }
				 /* 
		 if we have fixed electrons, calculate extra distribution now
				  */
				 if (hasFixedPositions)
				 {
					 glLogger.info("Adding values for fixed electron");
					 for (unsigned int fixIndex = 0; fixIndex < fixedPos.size(); fixIndex++)
					 {
						 glLogger.info("Fixed Electron with index (%d), at (%f), (%f)",
								 fixIndex, fixedPos.at(fixIndex).getXPosition(),
								 fixedPos.at(fixIndex).getYPosition());
						 CxPeriodicPosition elecPos (fixedPos.at(fixIndex), 1.0f);
						 // Loop over all vortices
						 std::vector<CxVortex>::iterator vortexIt = vortexPos.begin();
						 while (vortexIt !=vortexPos.end())
						 {
							 // Calculate distance add to distribution TODO!
							 glLogger.info("Adding vortex Position x +(%f), y=(%f)",
									 vortexIt->getXPosition(), vortexIt->getYPosition());
							 fixedD[fixIndex].addValue((elecPos - (*vortexIt)).vabs());
							 fixed2d[fixIndex].addValue((*vortexIt));
							 vortexIt++;
						 }

					 }
					 glLogger.info("Finished fixed electron");

				 }
			 }
		 }
		 else
		 {
			 for (unsigned int elecIndex = 0; elecIndex < electrons.size(); elecIndex++)
			 {
				 glLogger.debug("Electron pos loop, electron No. (%d) at (%f), (%f)",
						 elecIndex,
						 electrons[elecIndex].getXPosition(),
						 electrons[elecIndex].getYPosition());
				 CxPeriodicPosition elecPos (electrons[elecIndex]);

				 for (unsigned int elecIndex2 = elecIndex+1; elecIndex2 < electrons.size(); elecIndex2++)
				 {
					 CxPeriodicPosition elecPos2(electrons[elecIndex2]);

					 elecElecD.addValue((elecPos2-elecPos).vabs());
				 }
			 }
		 }

		 glLogger.error("Monte carlo step (%d)", mcCounter);
		 if (mcCounter >  monteCarloSteps)
		 {
			 loopControl = false;
		 }
		 else 
		 {
			 loopControl = true;
		 }
		 if (rint(mcCounter/1000)*1000 == mcCounter) // log temporary fileMo
		 {
			 // Write temporary files
			 if (eeCorrelationOnly == false)
			 {
				 std::ofstream eVStream("./evCorrelation_temp.dat");
				 eVStream << "# Electron vortex correlation \n";
				 eVStream << "# Loop : "<<mcCounter << "\n";
				 eVStream << "# analyzing state " << stateFileName <<std::endl;
				 elecVortexD.writeToStream(eVStream);
				 eVStream.close();
				 std::ofstream vorvorStream("./vortexVortexCorrelation_temp.dat");
				 vorvorStream << "# Vortex vortex correlation \n";
				 vorvorStream << "# Loop : "<<mcCounter << "\n";
				 vorvorStream << "# analyzing state " << stateFileName <<std::endl;
				 vortexVortexD.writeToStream(vorvorStream);
				 vorvorStream.close();
				 std::ofstream eV2dStream("./evCorrelation2d.dat");
				 eV2dStream << "# Electron vortex correlation \n";	
				 eV2dStream << "# Loop : "<<mcCounter << "\n";
				 eV2dStream << "# analyzing state " << stateFileName <<std::endl;
				 elecVortex2d.writeToStream(eV2dStream);
				 if (hasFixedPositions)
				 {
					 std::string fileNameBas("fixedElecVortexDistr");
					 for (unsigned int index = 0; index < fixedPos.size(); index++)
					 {
						 glLogger.error("Save fixed Elec values");
						 fileNameBas = fileNameBas.append("a");
						 std::string fileName = fileNameBas.append(".dat");
						 std::ofstream fixedDistStream(fileName.c_str());
						 fixedDistStream << "# electron vortex distribution for electron [" << index;
						 fixedDistStream << "] at position "<< fixedPos[index].getXPosition()<<" " << fixedPos[index].getYPosition() << std::endl;
						 fixedDistStream << "# analyzing state " << stateFileName <<std::endl;
						 fixedDistStream << "# Loop : "<<mcCounter << "\n";
						 fixedD[index].writeToStream(fixedDistStream);
						 fixedDistStream.close();
						 // raw vortices

						 std::string fileName2 = fileNameBas.append("_raw.dat");
						 std::ofstream fixedRawDistStream(fileName2.c_str());
						 fixedRawDistStream << "# electron vortex distribution for electron [" << index;
						 fixedRawDistStream << "] at position "<< fixedPos[index].getXPosition()<<" " << fixedPos[index].getYPosition() << std::endl;
						 fixedRawDistStream << "# analyzing state " << stateFileName << " all vortices" <<std::endl;
						 fixedRawDistStream << "# Loop : "<<mcCounter << "\n";

						 fixedRaw[index].writeToStream(fixedRawDistStream);
						 fixedRawDistStream.close();


						 // raw wo Puali
						 std::string fileName3 = fileNameBas.append("_raw_woPauli.dat");
						 std::ofstream fixedRawWoPauliDistStream(fileName3.c_str());
						 fixedRawWoPauliDistStream << "# electron vortex distribution for electron [" << index;
						 fixedRawWoPauliDistStream << "] at postion "<< fixedPos[index].getXPosition()<<" " << fixedPos[index].getYPosition() << std::endl;
						 fixedRawWoPauliDistStream << "# analyzing state " << stateFileName << " removed Pauli vortices" <<std::endl;
						 fixedRawWoPauliDistStream << "# Loop : "<<mcCounter << "\n";

						 fixedwoPauli[index].writeToStream(fixedRawWoPauliDistStream);
						 fixedRawDistStream.close();


						 // 2 d the same
						 fileName = fileNameBas.append("_2d");
						 fileName =  fileName.append("a.dat");
						 std::ofstream fix2dStream(fileName.c_str());
						 fix2dStream << "# electron vortex distribution for electron [" << index;
						 fix2dStream << "] at postion "<< fixedPos[index].getXPosition()<<" " << fixedPos[index].getYPosition() << std::endl;
						 fix2dStream << "# analyzing state " << stateFileName <<std::endl;
						 fix2dStream << "# Loop : "<<mcCounter << "\n";
						 fixed2d[index].writeToStream(fix2dStream);
						 fix2dStream.close();
					 }



				 } // end hasFixedPositions
			 }
			 std::ofstream eEStream("./elecElecCorrelation_temp.dat");
			 eEStream << "# Electron electron correlation \n";
			 eEStream << "# Loop : "<<mcCounter << "\n";
			 eEStream << "# analyzing state " << stateFileName <<std::endl;
			 elecElecD.writeToStream(eEStream);
			 eEStream.close();
		 }
	 }
	 // Done, now writing distributions
	 std::ofstream eVStream("./evCorrelation.dat");
	 eVStream << "# Electron vortex correlation \n";
	 eVStream << "# Total number of accepetd MC steps : "<<mcCounter << "\n";
	 eVStream << "# analyzing state " << stateFileName <<std::endl;
	 if (!elecVortexD.writeToStream(eVStream)) return (false);
	 eVStream.close();
	 std::ofstream vorvorStream("./vortexVortexCorrelation.dat");
	 vorvorStream << "# Total number of accepetd MC steps : "<<mcCounter << "\n";
	 vorvorStream << "# analyzing state " << stateFileName <<std::endl;
	 if (!vortexVortexD.writeToStream(vorvorStream)) return (false);
	 vorvorStream.close();
	 std::ofstream eEStream("./elecElecCorrelation.dat");
	 eEStream << "# Electron electron correlation \n";
	 eEStream << "# Loop : "<<mcCounter << "\n";
	 eEStream << "# analyzing state " << stateFileName <<std::endl;
	 elecElecD.writeToStream(eEStream);
	 eEStream.close();
	 std::ofstream eV2dStream("./evCorrelation2d.dat");
	 eV2dStream << "# Electron vortex correlation \n";	
	 eV2dStream << "# Loop : "<<mcCounter << "\n";
	 eV2dStream << "# analyzing state " << stateFileName <<std::endl;
	 elecVortex2d.writeToStream(eV2dStream);
 	std::ofstream vortex2dStream("./evCorrelation2d.dat");
	 vortex2dStream << "#  vortex distribution \n";	
	 vortex2dStream << "# Loop : "<<mcCounter << "\n";
	 vortex2dStream << "# analyzing state " << stateFileName <<std::endl;
	 vortex2d.writeToStream(vortex2dStream);
	vortex2dStream.close();
	

	 if (eeCorrelationOnly == false)
	 {
		 std::ofstream eVStream("./evCorrelation_temp.dat");
		 eVStream << "# Electron vortex correlation \n";
		 eVStream << "# Loop : "<<mcCounter << "\n";
		 eVStream << "# analyzing state " << stateFileName <<std::endl;
		 elecVortexD.writeToStream(eVStream);
		 eVStream.close();
		 std::ofstream vorvorStream("./vortexVortexCorrelation_temp.dat");
		 vorvorStream << "# Vortex vortex correlation \n";
		 vorvorStream << "# Loop : "<<mcCounter << "\n";
		 vorvorStream << "# analyzing state " << stateFileName <<std::endl;
		 vortexVortexD.writeToStream(vorvorStream);
		 vorvorStream.close();
	 }
	 if (hasFixedPositions)
	 {
		 std::string fileNameBas("fixedElecVortexDistr");
		 for (unsigned int index = 0; index < fixedPos.size(); index++)
		 {
			 glLogger.error("Save fixed Elec values");
			 std::string tempName =  fileNameBas.append("a");
			 std::string fileName = tempName.append(".dat");
			 std::ofstream fixedDistStream(fileName.c_str());
			 fixedDistStream << "# electron vortex distribution for electron [" << index;
			 fixedDistStream << "] at postion "<< fixedPos[index].getXPosition()<<" " << fixedPos[index].getYPosition() << std::endl;
			 fixedDistStream << "# analyzing state " << stateFileName <<std::endl;
			 fixedDistStream << "# Loop : "<<mcCounter << "\n";
			 fixedD[index].writeToStream(fixedDistStream);
			 fixedDistStream.close();
			 // raw vortices

			 std::string fileName2 = fileNameBas.append("_raw.dat");
			 std::ofstream fixedRawDistStream(fileName2.c_str());
			 fixedRawDistStream << "# electron vortex distribution for electron [" << index;
			 fixedRawDistStream << "] at postion "<< fixedPos[index].getXPosition()<<" " << fixedPos[index].getYPosition() << std::endl;
			 fixedRawDistStream << "# analyzing state " << stateFileName << " all vortices" <<std::endl;
			 fixedRawDistStream << "# Loop : "<<mcCounter << "\n";

			 fixedRaw[index].writeToStream(fixedRawDistStream);
			 fixedRawDistStream.close();


			 // raw wo Puali
			 std::string fileName3 = fileNameBas.append("_raw_woPauli.dat");
			 std::ofstream fixedRawWoPauliDistStream(fileName3.c_str());
			 fixedRawWoPauliDistStream << "# electron vortex distribution for electron [" << index;
			 fixedRawWoPauliDistStream << "] at postion "<< fixedPos[index].getXPosition()<<" " << fixedPos[index].getYPosition() << std::endl;
			 fixedRawWoPauliDistStream << "# analyzing state " << stateFileName << " removed Pauli vortices" <<std::endl;
			 fixedRawWoPauliDistStream << "# Loop : "<<mcCounter << "\n";

			 fixedwoPauli[index].writeToStream(fixedRawWoPauliDistStream);
			 fixedRawDistStream.close();


			 // 2 d the same
			 fileNameBas = fileNameBas.append("_2d");
			 fileName =  fileNameBas.append("a");
			 std::ofstream fix2dStream(fileName.c_str());
			 fix2dStream << "# electron vortex distribution for electron [" << index;
			 fix2dStream << "] at postion "<< fixedPos[index].getXPosition()<<" " << fixedPos[index].getYPosition() << std::endl;
			 fix2dStream << "# analyzing state " << stateFileName <<std::endl;
			 fix2dStream << "# Loop : "<<mcCounter << "\n";
			 fixed2d[index].writeToStream(fix2dStream);
			 fix2dStream.close();
		 }
	 } // end hasFixedPositions
	 std::cerr << "Leaving monteCarlo \n";

	 return (true);




}


bool readSpectrumFromState( const std::string & stateFileName,
		const std::string & basisFileName,
		bool isBinary)

{
	glLogger.info("Entering readSpectrumFromState");
	std::string sStatefile(stateFileName);



	/*
    Start constructing basis, reading state and positions
	 */

	glLogger.info("reading basis");
	FILE *basisFile = fopen(basisFileName.c_str(),"rb");

	if (!basisFile)
	{
		glLogger.error("Could not open basis file");

		throw (new CxFileNotFoundError("Could not open basis file", __FILE__,__LINE__));
	}
	yosBasis *myBasis = new yosBasis (basisFile, 0);
	fclose (basisFile);

	persist_eigSt eig (myBasis, isBinary);

	std::ifstream fsVec;
	//! \todo change that!


	if (isBinary)
	{
		glLogger.info ("Reading binary vector-file...\n");

		fsVec.open (stateFileName.c_str(), 
				std::ios::binary | std::ios::in);
		if ( !fsVec )
		{
			fsVec.close();
			throw (new CxFileNotFoundError(stateFileName.c_str(),__FILE__,__LINE__ ));

		}

		eig.set_binary(isBinary);
		int stateCount = 0;
		while( !fsVec.eof()) 
		{
			fsVec >> (eig);
			std::cerr <<"  "<< eig.getEn()<<" "; 
			stateCount++;
			// Print state number and energy here!
		}
		std::cerr <<std::endl;

	}

	fsVec.close();
	return (true);
} 

/*!
  \brief calculates the minimal distance of Positions in electrons from Positions in vortices
  can be used to find the displacement of Pauli vortices from electrons

 */
float getMinimalDistance ( std::vector<CxPeriodicPosition> & electrons, 
		std::vector<CxPeriodicPosition> & vortices)
{
	double retVal = 1.0;
	glLogger.info("Entering getMinimalDistance");
	for (std::vector<CxPeriodicPosition>::iterator elecIt = electrons.begin(); 
	elecIt != electrons.end();
	elecIt++)
	{
		glLogger.info("Elec it x(%f), y(%f), xCell, (%f), yCell (%f)", 
				elecIt->getXPosition(), elecIt->getYPosition(),
				elecIt->get_xCellSize(), elecIt->get_yCellSize());
		for (std::vector<CxPeriodicPosition>::iterator vortexIt = vortices.begin();
		vortexIt != vortices.end();
		vortexIt++)
		{
			glLogger.info("Vortex it x(%f), y(%f), xCell, (%f), yCell (%f)", 
					vortexIt->getXPosition(), vortexIt->getYPosition(),
					vortexIt->get_xCellSize(), vortexIt->get_yCellSize());
			double differenz = ((*vortexIt) - (*elecIt)).vabs();
			if (differenz < retVal)
			{
				retVal = differenz;
			}
		}
	}

	glLogger.error("getMinimalDistance returns (%f)", retVal);
	return (retVal);
}

int removeAllPositionsWithinX(  std::vector<CxPeriodicPosition> & electrons, 
		std::vector<CxPeriodicPosition> & vortices,
		const double & distance)
{
	int removeCount = 0;
	glLogger.error("entering removeAllPositionswithin (%f)", distance);
	for (std::vector<CxPeriodicPosition>::iterator elecIt = electrons.begin(); 
	elecIt != electrons.end();
	elecIt++)
	{
		glLogger.error("Electron position (%f), (%f), xCell (%f), yCell (%f), distance = (%f)", 
				elecIt->getXPosition(),
				elecIt->getYPosition(), 
				elecIt->get_xCellSize(),
				elecIt->get_yCellSize(),
				distance);


		std::vector<CxPeriodicPosition>::iterator vortexIt = vortices.begin();


		while (vortexIt != vortices.end())
		{
			glLogger.error("Vortex (%f), (%f), cell (%f) x (%f)",
					vortexIt->getXPosition(),
					vortexIt->getYPosition(),
					vortexIt->get_xCellSize(),
					vortexIt->get_yCellSize());
			double differenz = ((*vortexIt) - (*elecIt)).vabs();
			glLogger.error("Differenz is %f", differenz);
			if (fabs(differenz) < fabs( distance))
			{
				glLogger.error("Remove an vortex from list (%f), (%f)", 
						vortexIt->getXPosition(),
						vortexIt->getYPosition());

				vortices.erase(vortexIt);
				removeCount++;
				vortexIt = vortices.begin();
			}
			else 
			{

				vortexIt++;
			}

		}
	}
	return(removeCount);
}

int removeMostDistantVortices ( std::vector<CxPeriodicPosition> & electrons, 
		std::vector<CxPeriodicPosition> & vortices,
		const int & numOfVor)
{
	glLogger.debug("Entering removeMostDistantVortices with (%d) vortices", numOfVor);
	double maxDistance = 0.0;
	if (numOfVor > (int)vortices.size())
	{
		glLogger.error("Too many vortices requestesd to be deleted (%d), have (%d)",
				numOfVor, vortices.size());
		return 0;
	}
	// find maximal Distance, store position in an iterator
	for (int count = 0; count < numOfVor ; count ++)
	{

		std::vector<CxPeriodicPosition>::iterator maximalVortex = vortices.end();
		std::vector<CxPeriodicPosition>::iterator vortIndex = vortices.begin();
		while (vortIndex != vortices.end())
		{
			int eIndex = vortIndex->getNearestPosition(electrons);
			double diff = ((*vortIndex) - electrons.at(eIndex)).vabs();
			glLogger.info("Found an vortex with distance (%f) at pos (%d)", diff, eIndex);
			if (diff > maxDistance)
			{
				glLogger.info("Found an vortex with distance (%f)", diff);
				maximalVortex = vortIndex;
			}
			vortIndex ++;
		}

		glLogger.info("Remove vortex at (%f), (%f)", 
				maximalVortex->getXPosition(),
				maximalVortex->getYPosition());
		vortices.erase(maximalVortex);

	}

	//  remove 

	return (numOfVor);
}


//! Now for vortices

/*!
  \brief calculates the minimal distance of Positions in electrons from Positions in vortices
  can be used to find the displacement of Pauli vortices from electrons

 */
float getMinimalDistance ( std::vector<CxPeriodicPosition> & electrons, 
		std::vector<CxVortex> & vortices)
{
	double retVal = 1.0;
	glLogger.info("Entering getMinimalDistance");
	for (std::vector<CxPeriodicPosition>::iterator elecIt = electrons.begin(); 
	elecIt != electrons.end();
	elecIt++)
	{
		glLogger.info("Elec it x(%f), y(%f), xCell, (%f), yCell (%f)", 
				elecIt->getXPosition(), elecIt->getYPosition(),
				elecIt->get_xCellSize(), elecIt->get_yCellSize());
		for (std::vector<CxVortex>::iterator vortexIt = vortices.begin();
		vortexIt != vortices.end();
		vortexIt++)
		{
			glLogger.info("Vortex it x(%f), y(%f), xCell, (%f), yCell (%f)", 
					vortexIt->getXPosition(), vortexIt->getYPosition(),
					vortexIt->get_xCellSize(), vortexIt->get_yCellSize());
			double differenz = (static_cast<CxPeriodicPosition>(*vortexIt) - (*elecIt)).vabs();
			if (differenz < retVal)
			{
				retVal = differenz;
			}
		}
	}

	glLogger.error("getMinimalDistance returns (%f)", retVal);
	return (retVal);
}

int removeAllPositionsWithinX(  std::vector<CxPeriodicPosition> & electrons, 
		std::vector<CxVortex> & vortices,
		const double & distance)
{
	int removeCount = 0;
	glLogger.error("entering removeAllPositionswithin (%f) w voritces", distance);
	for (std::vector<CxPeriodicPosition>::iterator elecIt = electrons.begin(); 
	elecIt != electrons.end();
	elecIt++)
	{
		glLogger.info("Electron position (%f), (%f), xCell (%f), yCell (%f), distance = (%f)", 
				elecIt->getXPosition(),
				elecIt->getYPosition(), 
				elecIt->get_xCellSize(),
				elecIt->get_yCellSize(),
				distance);


		std::vector<CxVortex>::iterator vortexIt = vortices.begin();


		while (vortexIt != vortices.end())
		{
			glLogger.debug("Vortex (%f), (%f), cell (%f) x (%f)",
					vortexIt->getXPosition(),
					vortexIt->getYPosition(),
					vortexIt->get_xCellSize(),
					vortexIt->get_yCellSize());
			double differenz = (static_cast<CxPeriodicPosition>(*vortexIt) - (*elecIt)).vabs();
			if (fabs(differenz) < fabs( distance))
			{
				glLogger.debug("Remove an vortex from list (%f), (%f)", 
						vortexIt->getXPosition(),
						vortexIt->getYPosition());

				vortices.erase(vortexIt);
				removeCount++;
				vortexIt = vortices.begin();
			}
			else 
			{
				vortexIt++;
			}

		} // end while vortexit
	} // end fo 

	return(removeCount);
}

int removeMostDistantVortices ( std::vector<CxPeriodicPosition> & electrons, 
		std::vector<CxVortex> & vortices,
		const int & numOfVor)
{
	glLogger.debug("Entering removeMostDistantVortices with (%d) vortices", numOfVor);
	double maxDistance = 0.0;
	if (numOfVor > (int)vortices.size())
	{
		glLogger.error("Too many vortices requestesd to be deleted (%d), have (%d)",
				numOfVor, vortices.size());
		return 0;
	}
	// find maximal Distance, store position in an iterator
	for (int count = 0; count < numOfVor ; count ++)
	{

		std::vector<CxVortex>::iterator maximalVortex = vortices.end();
		std::vector<CxVortex>::iterator vortIndex = vortices.begin();
		while (vortIndex != vortices.end())
		{
			int eIndex = vortIndex->getNearestPosition(electrons);
			double diff = (static_cast<CxPeriodicPosition>(*vortIndex) - electrons.at(eIndex)).vabs();
			glLogger.debug("Found an vortex with distance (%f) at pos (%d)", diff, eIndex);
			if (diff > maxDistance)
			{
				glLogger.info("Found an vortex with distance (%f)", diff);
				maximalVortex = vortIndex;
			}
			vortIndex ++;
		} // end of while

		glLogger.debug("Remove vortex at (%f), (%f)", maximalVortex->getXPosition(), maximalVortex->getYPosition());
		vortices.erase(maximalVortex);

	} // end of for loop 

	//  remove 

	return (numOfVor);
}
void writeSpectrum (std::vector<double> spectrum, std::string fileName)
{
// file open
	std::ofstream specFile (fileName.c_str());

	for(std::vector<double>::iterator it = spectrum.begin(); it != spectrum.end(); ++it) {


		specFile << (double)*it;
		it++;
	}
	specFile.close();


}
