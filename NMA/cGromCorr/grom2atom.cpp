//! String and io management modules
#include <iostream>
#include <fstream>
#include <regex>

//! Exception modules
#include <stdexcept>
#include <assert.h>

//! Math modules
#include <cmath>
#include <algorithm>
#include <vector>

#include "grom2atom.h"



void ReadOneTimeGrom2Atom(std::streampos const * const  currentPos, std::streampos & outPos, std::vector<float > & Xi, std::vector<std::string > & atomTypei, float & systemSize_nm, std::string const& fileName){

	std::ifstream gromFile;
	gromFile.open(fileName);

	// Check for nullpointer
	if (currentPos != 0){
		gromFile.seekg(*currentPos);
	}
	
	std::string line;

	int atomTotal = nan("");
	int atomCount = 0;

	bool encounteredHeader = false;

	if (gromFile.is_open()){
		while (getline(gromFile,line)){
			
			std::regex time_regex(".*t=.*");
			std::regex atom_regex(".*SOL.*");

	// 1
			if (std::regex_match(line, time_regex)){
				if (encounteredHeader){
					assert (std::none_of(Xi.begin(), Xi.end(), [](float x){return std::isnan(x);}));
					assert (atomCount == atomTotal);
					break;
				}

				// When you encounter the time, take the following line as the number of atoms in the timeslice
				encounteredHeader = true;
				std::cout << line << std::endl;
				if (getline(gromFile,line)){
					try{
						atomTotal = Convert<int>(line);
					}
					catch(std::invalid_argument & e){
						throw InputFileMalformed("Encountered non-integer entry when number-of-atoms was excpected after time-header line.");
					}
				}
				else{
					throw InputFileMalformed("Encountered EOF when number-of-atoms was expected after time-header line.");
				}

				Xi.assign(atomTotal*3,nan(""));
				atomTypei.assign(atomTotal,"");
			}

	// 2
	//
	// WARNING: This section works by chance, because the input format is very spectific. Search here if bugs arise.
			else if (!std::regex_match(line, atom_regex) && encounteredHeader){
				std::cout << line << std::endl;

				// Extract system size
				std::istringstream iss(line);
				std::string sub;
				iss >> sub;
				systemSize_nm = std::atof(sub.c_str());

				// Set stream pos after this line
				std::streampos outPos_try;
				outPos_try = gromFile.tellg(); //Can throw error
				outPos = outPos_try;
			}

	// 3 
			else if (std::regex_match(line, atom_regex) && encounteredHeader){
				std::istringstream iss(line);

				int loop_count = 0;
				do{
					std::string sub;
					iss >> sub;
					if (loop_count == 1) atomTypei[atomCount] = sub;
					//std::cout << "Substring: " << sub << std::endl;
					if (loop_count == 3)  Xi[atomCount*3 + 0] = std::atof(sub.c_str());
					if (loop_count == 4)  Xi[atomCount*3 + 1] = std::atof(sub.c_str());
					if (loop_count == 5)  Xi[atomCount*3 + 2] = std::atof(sub.c_str());
					//if (loop_count == 6)  Vi[atomCount][0] = sub;
					//if (loop_count == 7)  Vi[atomCount][0] = sub;
					//if (loop_count == 8)  Vi[atomCount][0] = sub;
					loop_count++;
				} while(iss);

				//std::cout << Xi[atomCount][0] << std::endl;
				atomCount++;
			}
	// End switch 
		}
		gromFile.close();
	}

	//If file opening failed:
	else{
		throw ReadError("Open failed on file passed in.");
	}
}




void GenerateMassAndSize(std::vector<std::string > const& atomTypei, std::vector<float > & atomMassi_amu, std::vector<float > & atomSizei_nm){

	atomMassi_amu.assign(atomTypei.size(), 0);
	atomSizei_nm.assign(atomTypei.size(), 0);

	for (int i = 0 ; i < atomTypei.size() ; i++){
		if (atomTypei[i].compare("OW") == 0){
			atomMassi_amu[i] = 15.9994;
			atomSizei_nm[i] = .3;
		}
		else if (atomTypei[i].compare("HW1") == 0){
			atomMassi_amu[i] = 1.008;
			atomSizei_nm[i] = .2;
		}
		else if (atomTypei[i].compare("HW2") == 0){
			atomMassi_amu[i] = 1.008;
			atomSizei_nm[i] = .2;
		}
		else{
			throw InputFileMalformed("Atom name " + atomTypei[i] + " from input file not found in lookup table.");
		}
	}
}

float DensityContributionErf(float dist_nm[3], float systemSize_nm, float halfBin_nm, float atomSize_nm, float CUTOFFsq_nm2){
	//Compute nearest periodic neighbors
	float halfBox_nm = systemSize_nm / 2;
	bool isCutoff = false;
	float dsq_nm2 = 0;

	for (int hat = 0 ; hat < 3 ; hat++){
		dist_nm[hat] = std::fmod(dist_nm[hat], systemSize_nm);

		if (dist_nm[hat] < -halfBox_nm){
			dist_nm[hat] += systemSize_nm;
		}
		if (dist_nm[hat] > halfBox_nm){
			dist_nm[hat] -= systemSize_nm; 
		}

		dsq_nm2 = dist_nm[hat] * dist_nm[hat];
		if (dsq_nm2 >= CUTOFFsq_nm2){
			return 0;
		}
	}

	// Compute the contribution to the density field
	float densityContribution = 1.;
	float halfBin_red = halfBin_nm / atomSize_nm;

	for (int hat = 0 ; hat < 3 ; hat++){

		float absDist_red = std::abs(dist_nm[hat])/atomSize_nm;

		densityContribution *= - std::erf(-halfBin_red - absDist_red ) 
												 	 + std::erf(+halfBin_red - absDist_red );
		
		//if (atomCount < 10 && boxCount < 10){
		//	std::cout << "absDist_red["	<< hat << "]: "<<  absDist_red << ", " << densityContribution << "\t";
		//}
		
	}
	return densityContribution / 8.;

}

void ComputeErfDensity_v1(std::vector<float > const & Xi_nm,
													std::vector<float > const & atomMassi_amu,
													std::vector<float > const & atomSizei_nm,
													float systemSize_nm,
													int Nb,
													std::vector<float > & density_amu){

	density_amu.assign(Nb*Nb*Nb, 0);
	double CUTOFF_nm = 1;
	double CUTOFFsq_nm2= CUTOFF_nm * CUTOFF_nm;

	float halfBox_nm = systemSize_nm / 2.;
	float halfBin_nm = halfBox_nm / Nb;


	std::cout << "SystemSize (nm), HalfSystem (nm):"
						<< systemSize_nm << "," << halfBox_nm << std::endl ;


	for (int atomCount = 0 ; atomCount < Xi_nm.size()/3 ; atomCount++){
		for (int boxCount = 0 ; boxCount < density_amu.size() ; boxCount++){
			float dist_nm[3] = {};

			dist_nm[0] = (float(boxCount % Nb)      * (systemSize_nm/Nb)) - Xi_nm[atomCount*3 + 0];
			dist_nm[1] = (float(boxCount / Nb % Nb) * (systemSize_nm/Nb)) - Xi_nm[atomCount*3 + 1];
			dist_nm[2] = (float(boxCount / Nb / Nb) * (systemSize_nm/Nb)) - Xi_nm[atomCount*3 + 2];
			
/* Refactored into DensityContributionErf():
			//Compute nearest periodic neighbors
			bool isCutoff = false;
			float dsq_nm2 = 0;

			for (int hat = 0 ; hat < 3 ; hat++){
				dist_nm[hat] = std::fmod(dist_nm[hat], systemSize_nm);

				if (dist_nm[hat] < -halfBox_nm){
					dist_nm[hat] += systemSize_nm;
				}
				if (dist_nm[hat] > halfBox_nm){
					dist_nm[hat] -= systemSize_nm; 
				}

				dsq_nm2 = dist_nm[hat] * dist_nm[hat];
				if (dsq_nm2 >= CUTOFFsq_nm2){
					isCutoff = true;
				}
			}

			// Compute the contribution to the density field
			if (!isCutoff){
				float densityContribution = 1.;
				float halfBin_red = halfBin_nm / atomSizei_nm[atomCount];

				for (int hat = 0 ; hat < 3 ; hat++){

					float absDist_red = std::abs(dist_nm[hat])/atomSizei_nm[atomCount];

					densityContribution *= - std::erf(-halfBin_red - absDist_red ) 
															 	 + std::erf(+halfBin_red - absDist_red );
					
					//if (atomCount < 10 && boxCount < 10){
					//	std::cout << "absDist_red["	<< hat << "]: "<<  absDist_red << ", " << densityContribution << "\t";
					//}
					
				}
				densityContribution /= 8.;
				densityContribution *= atomMassi_amu[atomCount];
	
				density_amu[boxCount] += densityContribution;

				//if (atomCount < 10 && boxCount < 10){
				//	std::cout << "atomCount, boxCount, <dist>:"	<< atomCount << ", " << boxCount  << ", "
				//						<< "<" << dist_nm[0] << ", " << dist_nm[1] << ", " << dist_nm[2] << ">: "
				//						<< densityContribution << std::endl;
				//}
			}
*/

			float densityContribution = DensityContributionErf(dist_nm, systemSize_nm, halfBin_nm, atomSizei_nm[atomCount],CUTOFFsq_nm2);
			density_amu[boxCount] +=  densityContribution * atomMassi_amu[atomCount];
		}
	}
}


void ComputeErfMomentum_v1(std::vector<float > const & Xi_nm,
													 std::vector<float > const & Vi_nm_ps,
													 std::vector<float > const & atomMassi_amu,
													 std::vector<float > const & atomSizei_nm,
													 float systemSize_nm,
													 int Nb,
													 std::vector<float > & momentum_amunm_ps){

	momentum_amunm_ps.assign(Nb*Nb*Nb*3, 0);
	double CUTOFF_nm = 1.;
	double CUTOFFsq_nm2= CUTOFF_nm * CUTOFF_nm;

	float halfBox_nm = systemSize_nm / 2.;
	float halfBin_nm = halfBox_nm / Nb;


	std::cout << "SystemSize (nm), HalfSystem (nm):"
						<< systemSize_nm << "," << halfBox_nm << std::endl ;


	for (int atomCount = 0 ; atomCount < Xi_nm.size()/3 ; atomCount++){
		for (int boxCount = 0 ; boxCount < Nb*Nb*Nb ; boxCount++){
			float dist_nm[3] = {};

			dist_nm[0] = (float(boxCount % Nb)      * (systemSize_nm/Nb)) - Xi_nm[atomCount*3 + 0];
			dist_nm[1] = (float(boxCount / Nb % Nb) * (systemSize_nm/Nb)) - Xi_nm[atomCount*3 + 1];
			dist_nm[2] = (float(boxCount / Nb / Nb) * (systemSize_nm/Nb)) - Xi_nm[atomCount*3 + 2];
			
			float densityContribution = DensityContributionErf(dist_nm, systemSize_nm, halfBin_nm, atomSizei_nm[atomCount],CUTOFFsq_nm2);
			densityContribution *= atomMassi_amu[atomCount];

			for (int hat = 0 ; hat < 3 ; hat++){
				momentum_amunm_ps[boxCount*3 + hat] += densityContribution * Vi_nm_ps[atomCount*3+ hat];
			}

		}
	}
}
























/*
void AtomLookup::AddAtom(std::string atomName, float atomMass_amu, float atomSize_nm){
	m_atomNames.push_back(atomName);
	m_atomMass_amu.push_back(atomMass_amu);
	m_atomSize_nm.push_back(atomSize_nm);
	assert(m_atomNames.size() == m_atomMass_amu.size() && m_atomNames.size() == m_atomSize_nm.size());
}

void AtomLookup::GetData(std::string const& atomNameIn, float & massOut_amu, float & sizeOut_nm){
	int i = 0;
	while (i < m_atomNames.size() && ! atomNames[i].compare(atomNameIn)){
		i++;
	}

	if (i == m_atomNames.size()){
		throw InputFlieMalformed("Atom name " << atomNameIn << " from input file not found in lookup table.")
	}
	massOut_amu = m_atomMass_amu[i];
	sizeOut_nm  = m_atomSize_nm[i];
}

std::vector<float > AtomLookup::ConvertNameToMass(std::vector<std::string > & atomNames){
	std::vector<float > masses(atomNames.size(), 0);
	return masses;
}

std::vector<float > AtomLookup::ConvertNameToSize(std::vector<std::string > & atomNames){
	std::vector<float > sizes(atomNames.size(), 0);
	return sizes;
}
*/

