#ifndef __GROM2ATOM_H_DEF__
#define __GROM2ATOM_H_DEF__


class InputFileMalformed : public std::logic_error
{ public: explicit InputFileMalformed(const std::string& what_arg) : std::logic_error(what_arg){} };

class ReadError : public std::runtime_error
{ public: explicit ReadError(const std::string& what_arg) : std::runtime_error(what_arg){} };


//! Convert: Converts a string into another data format, and throws std::invalid_argument if the string was malformed.
template< typename T > inline T Convert(const std::string& str)
{
	std::istringstream iss(str);
	T obj;
			
	iss >> std::ws >> obj >> std::ws;

	if(!iss.eof()) throw std::invalid_argument("Convert was not given an appropriate string");

	return obj; 
}



/*! ReadOneTimeGrom2Atom: reads a time slice from file into Xi, Vi, type_i, and 
 *	
 *	\input 	std::stremapos * currentPos 	Position in file to start from. 0 if beginning of file
 	*					std::stremapos & outPos 			Reference to output variable, outputs the non-atom line before the next time-header
 *					std::vector<float> Xi					Position vector to be populated, all data is overwritten
 *					std::vector<string> atomTypei	Atom type vector to be populated, all data is overwritten
 *					float				systemSize_nm			The output system size for the current frame
 *					std::string filename 					Name of the file to read from
 *
 *	\throws ReadError
 *					InputFileMalformed
 *			basic exception guarantee
 *
 *
 *	\return void
 */
void ReadOneTimeGrom2Atom(std::streampos const* const currentPos, 
													std::streampos & outPos,
												 	std::vector<float > & Xi,
												 	std::vector<std::string > & atomTypei,
													float &systemSize_nm,
												 	std::string const& fileName = "md_traj.gro");


/*! GenerateMassAndSize: Converts the array of atom types into arrays of their physical masses and "sizes"
 *
 * \input 	atomTypei 
 * 					atomMassi_amu : REQ size must be equal to atomTypei size
 * 					atomSize_nm : REQ size must be equal to atomTypei size
 *
 * \throws  InputFileMalformed
 *		basic exception guarantee
 *
 * */
void GenerateMassAndSize(std::vector<std::string > const& atomTypei,
												 std::vector<float > & atomMassi_amu,
												 std::vector<float > & atomSizei_nm);


/*! ComputeErfDensity_v1: first implementation of a coarsegraining procedure 
 *
 * \input 	std::vector<float > const & Xi_nm						The vector of atom positions, in nm
						std::vector<float > const & atomMassi_amu 	The vector of atom masses, in amu
						std::vector<float > const & atomSizei_nm		The vector of atom sizes, in nm
						float systemSize_nm													The length of the system (assumes cubic)
						int Nb																			The number of boxes to use for a single length
						std::vector<float > & density_amu						The output density field for this configuration, all data is overwritten

 * \throws  Nothing new; std::vectors may throw.
 *
 * */
void ComputeErfDensity_v1(std::vector<float > const & Xi_nm,
													std::vector<float > const & atomMassi_amu,
													std::vector<float > const & atomSizei_nm,
													float systemSize_nm,
													int Nb,
													std::vector<float > & density_amu);


/*! ComputeErfDensity_v1: first implementation of a coarsegraining procedure 
 *
 * \input 	std::vector<float > const & Xi_nm						The vector of atom positions, in nm
          	std::vector<float > const & Vi_nm_ps			  The vector of atom positions, in nm per ps
						std::vector<float > const & atomMassi_amu 	The vector of atom masses, in amu
						std::vector<float > const & atomSizei_nm		The vector of atom sizes, in nm
						float systemSize_nm													The length of the system (assumes cubic)
						int Nb																			The number of boxes to use for a single length
						std::vector<float > & momentum_amunm_ps			The output density field for this configuration, all data is overwritten

 * \throws  Nothing new; std::vectors may throw.
 *
 * */
void ComputeErfMomentum_v1(std::vector<float > const & Xi_nm,
													 std::vector<float > const & Vi_nm_ps,
													 std::vector<float > const & atomMassi_amu,
													 std::vector<float > const & atomSizei_nm,
													 float systemSize_nm,
													 int Nb,
													 std::vector<float > & momentum_amunm_ps);

													

























/*
struct AtomParams{
	std::vector<float > Xi_nm;
	std::vector<float > Vi_nm_ps;
	std::vector<std::string > atomTypei;
	float time_ps;
	float boxLength_nm;
	int atomTotal;
};

class AtomLookup{
	public:
		AtomLookup(): m_atomNames(), m_atomMass_amu(), m_atomSize_nm(){}
		void AddAtom(std::string const& atomName, float atomMass_amu, float atomSize_nm);
		 *! GetData 	inputs the mass and the size into the last two arguments by dictionary-style lookup.
		 *	
		 *	\throw 		InputFileMalformed
		 *		strong throw guarantee
		 *
		void GetData(std::string const& atomName, float & massOut_amu, float & sizeOut_nm);
		std::vector<float > ConvertNameToMass(std::vector<std::string > & atomNames);
		std::vector<float > ConvertNameToSize(std::vector<std::string > & atomNames);

	private:
		std::vector<std::string > m_atomNames;
		std::vector<float > m_atomMass_amu;
		std::vector<float > m_atomSize_nm;

}
*/

#endif //__GROM2ATOM_H_DEF__
