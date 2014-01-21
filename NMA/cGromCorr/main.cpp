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
#include <numeric>

//! cGromCorr modules
#include "grom2atom.h"


// --------------------------------- SYMM CORRELATION -------------------------------------------
// --------------------------------- SYMM CORRELATION -------------------------------------------
// --------------------------------- SYMM CORRELATION -------------------------------------------
// --------------------------------- SYMM CORRELATION -------------------------------------------
// --------------------------------- SYMM CORRELATION -------------------------------------------
// --------------------------------- SYMM CORRELATION -------------------------------------------

int main_gr(int n_samples){

	std::cout << "Computing <p(0)p(dr)>" << std::endl;
	bool posIsDefined = false;
	std::streampos currentPos;
	bool resign = false;

	int Nb = 10;
	int numBoxes = Nb * Nb * Nb;
	std::vector<float > density_av(Nb*Nb*Nb, 0);
	std::vector<float > prod_sum(Nb*Nb*Nb * Nb*Nb*Nb, 0);
	std::vector<float > corr_av(Nb*Nb*Nb * Nb*Nb*Nb, 0);


	for (int timeCount = 0 ; timeCount < n_samples ; timeCount++){

		
// --------------------------------- READ FILE -------------------------------------------
		
		std::vector<float > Xi;
		std::vector<std::string > atomTypei;
		float systemSize_nm = 0;
	


		try{
			std::streampos outputPos;
			std::streampos* ptrPos = posIsDefined? &currentPos : 0;
			ReadOneTimeGrom2Atom(ptrPos,
													 outputPos,
													 Xi, 
													 atomTypei,
													 systemSize_nm,
													 "water_traj.gro");
			currentPos = outputPos;
			posIsDefined = true;
		}
		catch (ReadError& e){
			assert(!resign);
			std::cout << "Caught ReadError: " << e.what() << "\n 	Trying to read the line again." <<std::endl;
			resign = true;
		}
		catch (InputFileMalformed& e){
			std::cout << "Caught InputFileMalformed: " << e.what() << "\n  Closing program. I should really tell you where, but I can't do that super well. The best you get is timeCount = "<<timeCount <<std::endl;
			break;
		}


// --------------------------------- CONVERT ATOMTYPE -------------------------------------------
		
		
		std::vector<float > atomMassi(atomTypei.size(), 0);
		std::vector<float > atomSizei(atomTypei.size(), 0);

		GenerateMassAndSize(atomTypei, atomMassi, atomSizei);


		
// --------------------------------- COMPUTE DENSITY -------------------------------------------

		std::vector<float > density_amu;

		ComputeErfDensity_v1(Xi, atomMassi, atomSizei, systemSize_nm, Nb, density_amu);

		float binsize_nm = (systemSize_nm/float(Nb));
		float scale_amu_TO_kg_m3 = 1.66 /binsize_nm /binsize_nm /binsize_nm;

		for (int boxCount = 0 ; boxCount < density_av.size() ; boxCount++){
			density_amu[boxCount] *= scale_amu_TO_kg_m3;
			density_av[boxCount]  += density_amu[boxCount];
		}

// --------------------------------- COMPUTE CORRELATION --------------------------------------------------

		float meanFrameDensity = std::accumulate(density_amu.begin(), density_amu.end(), 0.)/numBoxes;
		std::cout << "Mean frame density: " << meanFrameDensity << std::endl;

		for (int i = 0 ; i < numBoxes ; i++){
			for (int j = 0 ; j < numBoxes ; j++){
				corr_av[abs(i - j)] += (density_amu[i]-meanFrameDensity) * (density_amu[j]-meanFrameDensity);
			}
		}

	}

// --------------------------------- RESCALE AND CORRECT AVERAGES --------------------------------------
	
	for (int boxCount = 0 ; boxCount < numBoxes ; boxCount++){
		density_av[boxCount] /= n_samples;
	}
	
	float meanRMS = 0;
	for (int i = 0 ; i < numBoxes ; i++){
		corr_av[i] /= n_samples; 
	}
	meanRMS = corr_av[0];

	for (int i = 0 ; i < numBoxes ; i++){
		corr_av[i]  /= meanRMS;
	}


// --------------------------------- OUTPUT FILE -------------------------------------------

	std::ofstream outputStream;
	std::string outputName;
	
	outputName = "CsvOut/GrCorrelation.csv";
	outputStream.open(outputName, std::ofstream::out);
	for (int i = 0 ; i < density_av.size() ; i++){
		if (i!=0){
			outputStream << ",";
		}
		outputStream << density_av[i];
	}
	outputStream << std::endl;
	for (int i = 0 ; i < numBoxes;  i++){
		if (i!=0){
			outputStream << ",";
		}
		outputStream << corr_av[i];
	}
	outputStream << std::endl;
	//std::cout << std::endl;
	outputStream.close();
	std::cout << "Wrote the 3D density and [translationally symmetric] correlation in 1D" << std::endl;

	return 0;
}

// --------------------------------- ASYMM CORRELATION -------------------------------------------
// --------------------------------- ASYMM CORRELATION -------------------------------------------
// --------------------------------- ASYMM CORRELATION -------------------------------------------
// --------------------------------- ASYMM CORRELATION -------------------------------------------
// --------------------------------- ASYMM CORRELATION -------------------------------------------
// --------------------------------- ASYMM CORRELATION -------------------------------------------

int main_new(int n_samples){
	std::cout << "Computing <p(r1)p(r2)>" << std::endl;
	bool posIsDefined = false;
	std::streampos currentPos;
	bool resign = false;

	int Nb = 10;
	int numBoxes = Nb * Nb * Nb;
	std::vector<float > density_av(Nb*Nb*Nb, 0);
	std::vector<float > prod_sum(Nb*Nb*Nb * Nb*Nb*Nb, 0);
	std::vector<float > corr_av(Nb*Nb*Nb * Nb*Nb*Nb, 0);


	for (int timeCount = 0 ; timeCount < n_samples ; timeCount++){

		
// --------------------------------- READ FILE -------------------------------------------
		
		std::vector<float > Xi;
		std::vector<std::string > atomTypei;
		float systemSize_nm = 0;
	


		try{
			std::streampos outputPos;
			std::streampos* ptrPos = posIsDefined? &currentPos : 0;
			ReadOneTimeGrom2Atom(ptrPos,
													 outputPos,
													 Xi, 
													 atomTypei,
													 systemSize_nm,
													 "water_traj.gro");
			currentPos = outputPos;
			posIsDefined = true;
		}
		catch (ReadError& e){
			assert(!resign);
			std::cout << "Caught ReadError: " << e.what() << "\n 	Trying to read the line again." <<std::endl;
			resign = true;
		}
		catch (InputFileMalformed& e){
			std::cout << "Caught InputFileMalformed: " << e.what() << "\n  Closing program. I should really tell you where, but I can't do that super well. The best you get is timeCount = "<<timeCount <<std::endl;
			break;
		}


// --------------------------------- CONVERT ATOMTYPE -------------------------------------------
		
		
		std::vector<float > atomMassi(atomTypei.size(), 0);
		std::vector<float > atomSizei(atomTypei.size(), 0);

		GenerateMassAndSize(atomTypei, atomMassi, atomSizei);


		
// --------------------------------- COMPUTE DENSITY -------------------------------------------

		std::vector<float > density_amu;

		ComputeErfDensity_v1(Xi, atomMassi, atomSizei, systemSize_nm, Nb, density_amu);

		float binsize_nm = (systemSize_nm/float(Nb));
		float scale_amu_TO_kg_m3 = 1.66 /binsize_nm /binsize_nm /binsize_nm;

		for (int boxCount = 0 ; boxCount < density_av.size() ; boxCount++){
			density_amu[boxCount] *= scale_amu_TO_kg_m3;
			density_av[boxCount]  += density_amu[boxCount];
		}

// --------------------------------- COMPUTE CORRELATION --------------------------------------------------

		float meanFrameDensity = std::accumulate(density_amu.begin(), density_amu.end(), 0.)/numBoxes;
		std::cout << "Mean frame density: " << meanFrameDensity << std::endl;

		for (int i = 0 ; i < numBoxes ; i++){
			for (int j = 0 ; j < numBoxes ; j++){
				corr_av[(i * numBoxes) + j] += (density_amu[i]-meanFrameDensity) * (density_amu[j]-meanFrameDensity);
				prod_sum[(i * numBoxes) + j] += density_amu[i] * density_amu[j];
			}
		}

	}

// --------------------------------- RESCALE AND CORRECT AVERAGES --------------------------------------
	
	for (int boxCount = 0 ; boxCount < numBoxes ; boxCount++){
		density_av[boxCount] /= n_samples;
	}
	
	float meanRMS = 0;
	float meanRMS_v2 = 0;
	for (int i = 0 ; i < numBoxes ; i++){
		for (int j = 0 ; j < numBoxes ; j++){
			corr_av[(i * numBoxes) + j] /= n_samples; 
			prod_sum[(i * numBoxes) + j] /= n_samples;
			prod_sum[(i * numBoxes) + j] -= (density_av[i] * density_av[j]);
		}
		meanRMS 	 += corr_av [(i* numBoxes) + i];
		meanRMS_v2 += prod_sum[(i* numBoxes) + i];
	}
	meanRMS /= numBoxes;


	for (int i = 0 ; i < numBoxes ; i++){
		for (int j = 0 ; j < numBoxes ; j++){
			corr_av[(i * numBoxes) + j]  /= meanRMS;
			prod_sum[(i * numBoxes) + j] /= meanRMS_v2;
		}
	}


// --------------------------------- OUTPUT FILE -------------------------------------------

	std::ofstream outputStream;
	std::string outputName;
	
	outputName = "CsvOut/DensityCorrAvLin.csv";
	outputStream.open(outputName, std::ofstream::out);
	for (int i = 0 ; i < density_av.size() ; i++){
		if (i!=0){
			outputStream << ",";
		}
		outputStream << density_av[i];
	}
	outputStream << std::endl;
	for (int i = 0 ; i < numBoxes;  i++){
		for (int j = 0 ; j < numBoxes ; j++){
			if (j!=0){
				outputStream << ",";
			}
			outputStream << corr_av[(i * numBoxes) + j];
		}
		outputStream << std::endl;
	}
	outputStream << std::endl;
	for (int i = 0 ; i < numBoxes;  i++){
		for (int j = 0 ; j < numBoxes ; j++){
			if (j!=0){
				outputStream << ",";
			}
			outputStream << prod_sum[(i * numBoxes) + j];
		}
		outputStream << std::endl;
	}
	outputStream << std::endl;
	//std::cout << std::endl;
	outputStream.close();
	std::cout << "Wrote the linear stream, density and correlation" << std::endl;


	outputName = "CsvOut/densityAv2D.csv";
	outputStream.open( outputName, std::ofstream::out);
	for (int i = 0 ; i < Nb ; i++){
		for (int j = 0 ; j < Nb ; j++){
			for (int k = 0 ; k < Nb ; k++){
				outputStream << i << "\t" <<  j << "\t" << k << "\t" << density_av[i + j * Nb + k * Nb * Nb] << std::endl;
			}
		}
		outputStream<< "\n";
	}
	outputStream.close();
	std::cout << "Wrote the 3d Stream." << std::endl;


	return 0;
}


// --------------------------------- MAIN -------------------------------------------
// --------------------------------- MAIN -------------------------------------------
// --------------------------------- MAIN -------------------------------------------
// --------------------------------- MAIN -------------------------------------------
// --------------------------------- MAIN -------------------------------------------
// --------------------------------- MAIN -------------------------------------------


int main(int argc, char** argv){
	if (argc != 2){
		throw std::invalid_argument("You must supply an argument to ./parse specifying the number of lines to read.");
	}
	int n_entries = std::atoi(argv[1]);
	return main_gr(n_entries);
}
