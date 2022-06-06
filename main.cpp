/**
    Original version: March 2018 
    This version: June 2022
    Mario Romero
**/

#include <cmath>
#include <fstream>
#include <string>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "time.h"
#include "define.h"

//In this code, I define my own type, called my_float, to include both float and double precision.
//User can define FLOAT_PRECISION in define.h
#if FLOAT_PRECISION == 0
    #define my_float float
    #if FORCED_PRECISION <= 0
        #define NUMBER_PRECISION 8
    #else
        #define NUMBER_PRECISION FORCED_PRECISION
    #endif
//Next lines define the precision of number given in the output, not for the calculations
#else
    #define my_float double
    #if FORCED_PRECISION <= 0
        #define NUMBER_PRECISION 16
    #else
        #define NUMBER_PRECISION FORCED_PRECISION
    #endif
#endif

#define INTERPOLATOR Spline

#define NUM_EQUATIONS 4 //Number of differential equations to solve (ideally 4: mass, energy and momentum density plus one for modified entropy needed for dual energy formalism)

using namespace std;

//Define classes and functions
class Cell;
class Flux;
class Source;
class Slopes; //## CHECK IF YOU TRULY NEED SLOPES CLASS
class Spline;
class Bilinear;

//Functions of 'Physics/gasPhysics.cpp'
my_float pressure(my_float,my_float,my_float,my_float,my_float=0.0);
my_float pressure(my_float,my_float,my_float);
my_float temperature(my_float,my_float,my_float);
my_float number_density(my_float,my_float);
my_float internal_energy(my_float, my_float);
my_float specific_internal_energy(my_float, my_float, my_float);
my_float energy_density(my_float,my_float,my_float,my_float);
my_float entropy(my_float,my_float,my_float);
my_float magnetosonic_speed(my_float,my_float,my_float,my_float=0.0);

//Functions of 'Physics/thermodynamics.cpp'
void update_thermodynamics(Cell*, const int = NCELLS);

//Functions of 'Physics/hydrodynamics.cpp'
my_float max_value(my_float,my_float,my_float);
my_float min_value(my_float,my_float,my_float);
my_float slope_limiter(my_float,my_float,my_float);

//Functions of 'Schemes/timeEvolution.cpp'
my_float set_timestep(Cell,my_float);

//Functions of 'Schemes/boundaries.cpp'
void create_boundary(Cell&, Cell&, Cell&, Cell&, Cell*, const int = NCELLS);

//Functions of 'IO/read_input.cpp'
Bilinear load_table(string);

//Functions of 'Schemes/diagnostics.cpp'
void grid_diagnostics(Cell*,const int);

my_float slope_angle(my_float,my_float,my_float,my_float);
my_float relative_error(my_float,my_float);
bool in_front(Cell*,Cell*);

//Now add own files
#include "Physics/gasPhysics.cpp"
#include "Classes/cellclass.h"
#include "Classes/bilinear.h"
#include "Classes/fluxclass.h"
#include "Classes/internomials.h"
#include "Classes/Splineclass.h"
#include "Schemes/AMR_aux.cpp"
#include "Physics/thermodynamics.cpp"
#include "Physics/hydrodynamics.cpp"
#include "Schemes/interpolations.cpp"
#include "Schemes/AMR_main.cpp"
#include "Schemes/refinement_criteria.cpp"
#include "Schemes/timeEvolution.cpp"
#include "Schemes/boundaries.cpp"

#include "IO/read_input.cpp"
#include "IO/ConversorString.cpp"
#include "IO/set_ambient.cpp"
#include "IO/writeFile.cpp"
#include "Schemes/diagnostics.cpp"

int main(int argc, char** argv){

    time_t real_tstart = time(NULL);
    //Create number of available cells
    #if AMR_ON == 1
        const int N = NCELLS * pow(2,MAX_REF - INI_REF); //MAXIMUM number of cells (in case that all your cells have max refinement, which is unlikely.
    #else
        const int N = NCELLS;
    #endif
    const int Nneighbors = 4; //Number of neighbours (and ghost cells) to take into account

    //Check if input was given correctly
    if(argc < 4){
        cout<<"syntaxis: ./main [your file] [initial time in kyr] [end time in kyr] [other times (in kyr) in which you want an output]"<<endl;
        cout<<"Last array of times must be ordered and is optional to add"<<endl;
        throw runtime_error("Bad syntax");
    }
    const my_float year = 365.25*24.*3600.*1000.; //kyr -> s conversion
    my_float t0 = todouble(argv[2])*year; //Initial time
    my_float tf = todouble(argv[3])*year; //Last time
    
    my_float cstm_dt = (tf-t0)/1e3; //Custom timestep. It is an initial guess, and it can be less if stability of the code is compromised
    my_float t  = t0;
    my_float dt;
    my_float next_dt;
    my_float tmp_dt;
    const my_float AMR_delay_t = DELAY_AMR*year;

    int output_count = argc-4; //If you want more outputs than initial and end time, this number is positive and accounts how many files are left
    int curr_arg = 4;
    my_float next_output_t = ( curr_arg < argc ) ? todouble(argv[curr_arg])*year : t0-1.; //Second result is to make sure that you will never get an extra output

    //Declare all my code needed variables
    Cell* grid = new Cell[N];
    Cell* prd_grid = new Cell[N]; //Predictor values
    Cell* cor_grid = new Cell[N]; //Corrector values

    Flux* fluxes  = new Flux[Nneighbors]; //=[prev_E_flux, curr_W_flux, curr_E_flux, next_W_flux]

    //Ghost cells
    Cell first;
    Cell second;
    Cell penult;
    Cell last;
    
    //Initialize Equation of State (=Temperature table)
    init_EoS_table();
    #if ACTIVATE_COOLING == 1
        //Initialize cooling table
        init_cooling_table();
    #endif
    
    #if ACTIVATE_HEATING == 1
        ///Initialize heating table
        init_heating_table();
    #endif
    
    #if IC_FROM_FILE == 1
        init_from_file(grid,argv[1],NCELLS);
    #else
        t = set_ambient(grid,argv[1],NCELLS);
    #endif

    Cell initial_ambient = grid[N-1];

    //Let's create a data output for t=0
    #if DEBUG == 1
        grid_diagnostics(grid,N);
        cout<<"time = "<<t/year<<" kyr"<<endl;
        int test_number;
    #endif
    write_output(grid,initial_ambient,t/year);

    ///Get first timestep
    next_dt = cstm_dt;
    for(int i=0;i<N;i++){
        tmp_dt = set_timestep(grid[i],cstm_dt);

        next_dt = (tmp_dt > next_dt) ? next_dt : tmp_dt;
    }
    dt = next_dt;

    ///And pointers
    Cell* curr;
    Cell* prd_curr;
    Cell* cor_curr;
    curr = grid;
    
    prd_curr = prd_grid;
    cor_curr = cor_grid;
    
    #if CREATE_MAP == 1
        //Create a map of pressure, temperature, density, etc
        //This data is not used for the simulation
        create_map();
    #endif
    
    ///MAIN LOOP
    while(t < tf){
        /**
            PREDICTOR PHASE
        **/
        
        ///Init ghost cells
        create_boundary(first,second,penult,last,grid,N);
        
        //Create links with ghost cells
        make_a_link(&second,&first);
        make_a_link(&last,&penult);
        last.prev  = find_last(curr);
        first.next = find_first(curr);
        //^Note that I deliberately leave grid[0].prev = null and grid[N-1].next = null on purpose to know the chain borders!
        
        ///FIRST CELL
            //Get fluxes
            fluxes[0] = first.right_flux(second,*curr);
            fluxes[1] = curr->left_flux(first,*curr);
            fluxes[2] = curr->right_flux(first,*(curr->get_next()) );
            fluxes[3] = curr->next->left_flux(*curr,*(curr->next->get_next()) );

            //Compute the predictor

            prd_curr = find_first(prd_curr);
            *prd_curr = curr->predict(fluxes,dt);

            //Move to next cell
            curr = curr->next;

            prd_curr->Link_next(*(prd_curr+1));
            prd_curr = prd_curr->next;

        ///SECOND CELL
            //Get fluxes
            fluxes[0] = fluxes[2];
            fluxes[1] = fluxes[3];
            fluxes[2] = curr->right_flux( *(curr->get_prev()) , *(curr->get_next()) );
            fluxes[3] = curr->next->left_flux(*curr,*(curr->next->get_next()) );

            //Compute the predictor
            *prd_curr = curr->predict(fluxes,dt);
            
            //Move to next cell
            curr = curr->next;

            prd_curr->Link_prev(*(prd_curr-1));
            prd_curr->Link_next(*(prd_curr+1));
            prd_curr = prd_curr->next;
            
        ///COMMON CELLS
            while(curr->next->next != nullptr){
                //Get fluxes
                fluxes[0] = fluxes[2];
                fluxes[1] = fluxes[3];
                fluxes[2] = curr->right_flux( *(curr->get_prev()) , *(curr->get_next()) );
                fluxes[3] = curr->next->left_flux(*curr,*(curr->next->get_next()));

                //Compute the predictor
                *prd_curr = curr->predict(fluxes,dt);
                //Move to next cell
                curr = curr->next;

                prd_curr->Link_prev(*(prd_curr-1));
                prd_curr->Link_next(*(prd_curr+1));
                prd_curr = prd_curr->next;
            }

        ///PENULTIMATE CELL
            //Get fluxes
            fluxes[0] = fluxes[2];
            fluxes[1] = fluxes[3];
            fluxes[2] = curr->right_flux( *(curr->get_prev()) , *(curr->get_next()) );
            fluxes[3] = curr->next->left_flux(*curr,last);

            //Compute the predictor
            *prd_curr = curr->predict(fluxes,dt);
            //Move to next cell
            curr = curr->next;

            prd_curr->Link_prev(*(prd_curr-1));
            prd_curr->Link_next(*(prd_curr+1));
            prd_curr = prd_curr->next;

        ///LAST CELL
            //Get fluxes
            fluxes[0] = fluxes[2];
            fluxes[1] = fluxes[3];
            fluxes[2] = curr->right_flux(*(curr->get_prev()),last);
            fluxes[3] = last.left_flux(*curr,penult);

            //Compute the predictor
            *prd_curr = curr->predict(fluxes,dt);

            prd_curr->Link_prev(*(prd_curr-1));

            //Reset the 'curr' chain to the first cell
            curr = find_first(curr);
            
        /**
            CORRECTOR PHASE
        **/

        ///Init ghost cells
        create_boundary(first,second,penult,last,prd_grid,N);
        //Create links with ghost cells
        make_a_link(&second,&first);
        make_a_link(&last,&penult);
        last.prev  = prd_curr; //That's why I didn't reset 'prd_curr' pointer
        prd_curr = find_first(prd_curr);
        first.next = prd_curr; //And now it is
        cor_curr = find_first(cor_curr);

        ///FIRST CELL
            //Get fluxes
            fluxes[0] = first.right_flux(second,*prd_curr);
            fluxes[1] = prd_curr->left_flux(first, *(prd_curr->get_next()) );
            fluxes[2] = prd_curr->right_flux(first, *(prd_curr->get_next()) );
            fluxes[3] = prd_curr->next->left_flux(*prd_curr, *(prd_curr->next->get_next()) );

            //Compute the predictor
            *cor_curr = prd_curr->correct(*curr,fluxes,dt);

            //Get preferred timestep
            tmp_dt = set_timestep(*cor_curr,cstm_dt);
            next_dt = (tmp_dt > next_dt) ? next_dt : tmp_dt;

            //Update grid for next timestep
            curr->soft_update(*cor_curr); //soft_update does preserve pointers
            curr = curr->next;
            prd_curr = prd_curr->next;

            cor_curr->Link_next(*(cor_curr+1));
            cor_curr = cor_curr->next;

        ///SECOND CELL
            //Get fluxes
            fluxes[0] = fluxes[2];
            fluxes[1] = fluxes[3];
            fluxes[2] = prd_curr->right_flux(*(prd_curr->get_prev()) , *(prd_curr->get_next()));
            fluxes[3] = prd_curr->next->left_flux(*prd_curr , *(prd_curr->next->get_next()) );

            //Compute the predictor
            *cor_curr = prd_curr->correct(*curr,fluxes,dt);

            //Get preferred timestep
            tmp_dt = set_timestep(*cor_curr,cstm_dt);
            next_dt = (tmp_dt > next_dt) ? next_dt : tmp_dt;

            //Update grid for next timestep
            curr->soft_update(*cor_curr); //soft_update does preserve pointers
            curr = curr->next;
            prd_curr = prd_curr->next;

            cor_curr->Link_prev(*(cor_curr-1));
            cor_curr->Link_next(*(cor_curr+1));
            cor_curr = cor_curr->next;

        ///COMMON CELLS
            while(curr->next->next != nullptr){
                //Get fluxes
                fluxes[0] = fluxes[2];
                fluxes[1] = fluxes[3];
                fluxes[2] = prd_curr->right_flux(*(prd_curr->get_prev()) , *(prd_curr->get_next()));
                fluxes[3] = prd_curr->next->left_flux(*prd_curr , *(prd_curr->next->get_next()) );

                //Compute the predictor
                *cor_curr = prd_curr->correct(*curr,fluxes,dt);

                //Get preferred timestep
                tmp_dt = set_timestep(*cor_curr,cstm_dt);
                next_dt = (tmp_dt > next_dt) ? next_dt : tmp_dt;

                //Update grid for next timestep
                curr->soft_update(*cor_curr); //soft_update does preserve pointers
                curr = curr->next;
                prd_curr = prd_curr->next;

                cor_curr->Link_prev(*(cor_curr-1));
                cor_curr->Link_next(*(cor_curr+1));
                cor_curr = cor_curr->next;
            }
            

        ///PENULTIMATE CELL
            //Get fluxes
            fluxes[0] = fluxes[2];
            fluxes[1] = fluxes[3];
            fluxes[2] = prd_curr->right_flux(*(prd_curr->get_prev()) , *(prd_curr->get_next()));
            fluxes[3] = prd_curr->left_flux(*prd_curr , last);

            //Compute the predictor
            *cor_curr = prd_curr->correct(*curr,fluxes,dt);

            //Get preferred timestep
            tmp_dt = set_timestep(*cor_curr,cstm_dt);
            next_dt = (tmp_dt > next_dt) ? next_dt : tmp_dt;

            //Update grid for next timestep
            curr->soft_update(*cor_curr); //soft_update does preserve pointers
            curr = curr->next;
            prd_curr = prd_curr->next;

            cor_curr->Link_prev(*(cor_curr-1));
            cor_curr->Link_next(*(cor_curr+1));
            cor_curr = cor_curr->next;

        ///LAST CELL
            //Get fluxes
            fluxes[0] = fluxes[2];
            fluxes[1] = fluxes[3];
            fluxes[2] = prd_curr->right_flux(*(prd_curr->get_prev()) , last);
            fluxes[3] = last.left_flux(*prd_curr,penult);

            //Compute the predictor
            *cor_curr = prd_curr->correct(*curr,fluxes,dt);

            //Get preferred timestep
            tmp_dt = set_timestep(*cor_curr,cstm_dt);
            next_dt = (tmp_dt > next_dt) ? next_dt : tmp_dt;

            //Update grid for next timestep
            curr->soft_update(*cor_curr); //soft_update does preserve pointers
            cor_curr->Link_prev(*(cor_curr-1));

        /**
            UPDATE NEXT TIMESTEP
        **/

        if(AMR_delay_t > t){ refine_grid(grid,N,true); }
        else{ refine_grid(grid,N); }
        
        curr = find_first(grid); //find_first(curr) is dangerous
        dt = next_dt;
        t += dt;

        if(dt < SOFT_ZERO){
            throw runtime_error("Infinite loop! Timestep is too small");
        }
        
        //Reset next timestep
        next_dt = cstm_dt; 
        //^This is necessary, because when it encounters the smallest timestep in the whole run, this simulation will keep it until the end.

        /**
            IF APPLICABLE, GENERATE A FILE
        **/
        if( t >= next_output_t && curr_arg < argc){

            write_output(grid,initial_ambient,next_output_t/year); //Output file will give time in yr

            if(curr_arg != argc-1){ //Is not the last argument
                curr_arg++;
                next_output_t = todouble(argv[curr_arg])*year;
            }
            else{
                next_output_t = 2.*tf; //No more times
            }
            //If I want to account times...
            cout<<"time elapsed: "<<difftime(time(NULL),real_tstart)<<" s"<<endl;
        }

        #if DEBUG == 1
            #if STEP_BY_STEP ==1 //Put before everything else to allow me to stop with the initial condition
            cin>>test_number;

            if(test_number == 0){
                throw runtime_error("Test!");
            }
            #endif // STEP_BY_STEP
            grid_diagnostics(grid,N);
            cout<<"time = "<<t/year<<" kyr"<<endl;
            write_output(grid,initial_ambient,next_output_t/year);
        #endif
    }

    //Generate final output file
    write_output(grid,initial_ambient,t/year); //Output file will give time in yr

    //Destroy everything
    delete[] grid;
    delete[] prd_grid;
    delete[] cor_grid;

    delete[] fluxes;

    time_t real_tend = time(NULL);
    cout<<"time elapsed: "<<difftime(real_tend,real_tstart)<<" s"<<endl;

    return 0;
}
//That's all!
