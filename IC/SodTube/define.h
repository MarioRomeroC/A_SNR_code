#ifndef DEFINE_H_INCLUDE
#define DEFINE_H_INCLUDE

///DEBUGGING
//Written here only for the sake of less writting time

#define DEBUG 0
#define STEP_BY_STEP 0

///FILE DIRECTORIES

#define OUTPUT_PATH "Outputs/"  //Path to folder in which code outputs will be saved. put the '/' behind or weird things may happen.
#define KEY_ON 0		//If on, shows a key at first line in files
#define FORCED_PRECISION 4 //Number of digit that OUTPUT file will have. If set to 0 or lower, output precision is 8 and 16 digits for float and double, respectively.
               //This parameter does NOT affect code calculations, only outputs.

///INITIAL CONDITIONS

#define IC_FROM_FILE 1      // Self-explanatory. If set to 0, you still need a file to define an ambient (which is one row)

#define CENTRAL_EXPLOSION 1 //Defines a central explosion as a source term only at first iteration
#define MASS_ON 0 	    //If 0, no mass is added to density. Put 1 to add mass, but note that this feature (adding mass) is not fully tested.

#define SHOCK_RADIUS .1 //Shock radius, in pc, at starting time

#define DELAY_AMR 0. //Time, in kyr, in which some refinement criteria before that time will be ignored
                      //^^ If you start with a poor resolved central explosion (i.e.: ALWAYS) and you interpolate, the interpolation error will lead to a wrong result!

///AMR PARAMETERS
//To deactivate AMR, put 0 in MAX_REF
#define NCELLS  512   //Initial number of cells (at INI_REF, important!)
#define MAX_REF 0       //Maximum level of refinement (If 0, this is equivalent to turn AMR off)
#define INI_REF MAX_REF //Initial level of refinement of ALL cells. (I suggest starting with INI_REF = MAX_REF for maximum accuracy at initial time)

#define INTERPOL_RULE 0     //0 uses natural spline interpolation, 1 uses linear piecewise, and 2 both (spline for interior, linear for shock and ambient)
                              //^^ This only affects variables of NUM_EQUATIONS, other ones gets linear piecewise (they are results of the code, not variables to compute)

///ANALYSIS

#define ENTROPY_FOR_SHOCKS 0 //Use entropy instead of density to find shocks

///SIMULATION PARAMETERS
#define FLOAT_PRECISION 1   //If set to 0, grackle gr_type is 'float' and 'double' otherwise. I think is more useful for whole simulation.

#define SYMMETRY 0        //(0) Cartesian (1) Cylindrical and (2) Spherical coordinates

#define RBOUNDARY 0       //Boundary condition for last cell: (0) Reflection, (1) Copycat
#define LBOUNDARY 0       //Boundary condition for first cell. Options are the same as RBOUNDARY. If SYMMETRY>0, LBOUNDARY should be a reflection (0) to avoid funky issues

#define FAILURE -1        //If some subroutine fails (specifically functions from 'Physics/gasPhysics.cpp'), it will return this number
#define SOFT_ZERO 0.0     //In case you have troubles with hard zeros, you may change this number to a very small one (e.g.: 1e-200).
#define EPSILON_FLOAT 1E-4//To avoid numerical truncation errors if we compare two floating point numbers

#define CFL_NUMBER  0.3   //Should be always less than 0.5. Guarantees stability for hydrodinamics
#define CFL_NUMBER2 0.1   //Same as CFL_NUMBER, but it's used for cooling instead of hydrodynamics
#define ETA_E 0.0         //Used for dual energy formalism, if thermal energy / energy is less than this number, entropy system is used (have a number of 0.0 or higher to avoid negative pressures)

///PHYSICAL CONSTANTS
#define ADIAB_CONST 7./5.            //(initial) Adiabatic constant, gamma for friends
#define INIT_METALLICITY 0.02        //(DEPRECIATED) Initial metallicity of cells

#define MEANWEIGHT 0.85              //(DEPRECIATED) Molecular weight

#define PI 3.141592653589793

#define SOLAR_MASS 1.98847e33        // in units of g
#define HYDROGEN_MASS 1.67262171e-24 // in units of g. Mass of hydrogen
#define BOLTZMANN_CTE 1.3806504e-16  // in units of erg/K. Boltzmann constant
#define PARSEC 3.0856775814914E+18   // in units of cm

#define VACUUM_PERM 1.0              //Vacuum Permeability In CGS units (currently unused)

#define SOLAR_METALLICITY INIT_METALLICITY       //Self-explanatory

///EQUATION OF STATE
#define CREATE_MAP 0 //If 1, it creates a table of pressures, temperatures, densities and several quantities in a single file

#define EOS_FILE "IC/Example/EoS.txt" //"CoolingFunctions/EoS_V4p2.txt" //EoS_19.txt"  //Path to your equation of state file. First column must be density (g/cm3), second is pressure (erg/cm3) and third temperature (K)
					      //^If some cell is out of the range of this table, it will extrapolate
#define EOS_INTERPOLATION 0		//If 0, the EoS table is interpolated using the logarithm. Usual bilinear interpolation is 1		

///COOLING PARAMETERS
#define ACTIVATE_COOLING 0        //If 1, it reads the cooling file and uses it as a table
#define COOLING_FILE "IC/Example/Cooling.txt" //"CoolingFunctions/CoolingTable_V4p2.txt" //Cool_19.txt" //Cool_21.txt" //Path to your cooling rate file in erg/cm3/s. For temperatures out of the range given, it will extrapolate
							  //^The table you have to provide has density (g/cm3) as first column, pressure (erg/cm3) as second and the cooling rate (erg/cm3/s) as the third
							  //^^Also make sure that the cooling rate HAS POSITIVE VALUES.
								
#define COOLING_INTERPOLATION 0	  //If 0, the cooling table is interpolated using the logarithm. Usual bilinear interpolation is 1
#define MAX_COOLING_TEMPERATURE 1e10  //For temperatures above this threshold, the code will ignore cooling.
					//^For a Sedov explosion, you can obtain, at the first timesteps when the shock structure is forming, T>1e10. In absence of cooling, T decrease to reasonable values
					//^^However, if cooling is present, you can get a Sedov explosion with a lower initial energy due to the cooling of these T>1e10 K cells.
					//^^^If you have a cooling file that covers a different range of temperatures, please check this parameter!

///HEATING PARAMETERS
#define ACTIVATE_HEATING 0        //If 1, it reads the heating file and uses it as a table
#define HEATING_FILE "IC/Example/Heating.txt" //Path to your heating rate file in erg/s. For temperatures out of the range given, it will extrapolate
						    //^The table you have to provide has density (g/cm3) as first column, pressure (erg/cm3) as second and the heating rate (erg/s) as the third

#define HEATING_INTERPOLATION 1 	//If 0, the heating table is interpolated using the logarithm. Usual bilinear interpolation is 1
#define MAX_HEATING_TEMPERATURE MAX_COOLING_TEMPERATURE //For temperature above this threshold, the code will ignore heating
                                    //If you use the usual Koyama&Inutsuka02 heating, this is not needed. However, for complex heatings, this is to avoid overheating due to extrapolation

#define TEMPERATURE_CUTOFF 10.0        //For cells with lower temperatures, it will remove cooling/heating

#endif
