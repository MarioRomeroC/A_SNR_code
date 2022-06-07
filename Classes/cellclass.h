#ifndef CELLCLASS_H_INCLUDE
#define CELLCLASS_H_INCLUDE

/******
    Cell class stores the 'thermohydrodynamics' of the code

    There are also 'special' cells created for some calculations, such 'get_derivative' or 'get_sources'.
    They appear behind the curtains (that's why there's some private functions), and not in main to avoid confusion.
*******/

class Cell{

  //friends:
        //these are friend classes and functions, and have access to private values
        friend void update_thermodynamics(Cell*, const int); //Updates P, T, gamma and t_cool using grackle for the whole grid
        friend void initialize_fields(Cell*, const int);
        friend class Flux;
        friend class Source;

    private:
        //These are the magnitudes that follow a differential equation
        my_float y[NUM_EQUATIONS]; // = {mass density, TOTAL energy density, momentum density, entropy}

        //These are other magnitudes related with termodynammics, and do not enter into the differential equations explicitly
        my_float P;      //Pressure
        my_float T;      //Temperature
        my_float gamma;  //Adiabatic constant
        my_float Z;      //Metallicity (Unused)
        my_float t_cool; //Cooling time. It has a sign, negative means 'gas cools', and possitive 'gas heats'
        my_float meanweight; //Atomic mean weight
        
        bool system_E; //Needed for dual energy formalism: we use system_E by default, and system_S (entropy) otherwise

        //Finally we save location of said cells
        my_float r; //Location of the cell
        my_float dr; //Width of such cell

        
        int LoR = INI_REF; //Level of Refinement. All cells will be initialized with INI_REF refinement. Be aware of that when creating cells
        //Booleans to know if a cell need refinement or not
        bool flag_for_merge;
        bool flag_for_split;

        int flag_for_shock; //0 -> No shock; 1 -> I'm a shock; 2 -> I have a shock one cell in front/behind me; 3 -> I have a shock two cells in front/behind me; etc
        int flag_for_gradient; //Same as above, but with energy gradients.
        

        /// DO IMPORTANT STUFF
        //These functions appear inside the public functions that compute the predictor, average and corrector.
        void update_thermodynamics(); //Updates thermodynamic magnitudes for this cell. See 'Physics/thermodynamics.cpp' for the function
        void resync_thermodynamics();

        /*Next functions create memory, you should delete in the code (in predict and correct functions)*/
        my_float* get_derivative(Flux,Flux,Flux,Flux); //Get all y[elem] derivatives for this cell. See 'Physics/hydrodynamics.cpp' for the function
        my_float* get_derivative(Cell,Cell,Cell,Cell); //Get all y[elem] derivatives for this cell. See 'Physics/hydrodynamics.cpp' for the function


        //Cell get_sources(); //Get the sources of each equation. See 'Physics/hydrodynamics.cpp' for the function
        my_float* get_sources();


        void resync(){
            //Resyncronize E and S evolution, depending if system E is used or not
            if(system_E == true){
                //So entropy is resync
                y[3] = entropy(y[0],P,gamma);
            }
            else{
                //So energy density is resync
                y[1] = energy_density(y[0],this->get_velocity(),P,gamma);
                system_E = true;
            }
        }

        /// CELL RECONSTRUCTION
        //Spatially speaking
        Cell left_interface(Cell prev, Cell next){ //next is needed
            //In literature left interface is WEST (W)
            my_float data[NUM_EQUATIONS];
            my_float delta_u;
            for(int elem=0;elem<NUM_EQUATIONS;++elem){
                delta_u = slope_limiter(prev[elem],y[elem],next[elem]);
                data[elem] = y[elem] - delta_u; //Notice the sign. This is the change with 'right_interface'
            }
            //Thermodynamic variables must be updated accordingly
            delta_u = slope_limiter(prev.P,P,next.P);
            my_float face_P = P - delta_u;

            delta_u = slope_limiter(prev.T,T,next.T);
            my_float face_T = T - delta_u;

            delta_u = slope_limiter(prev.Z,Z,next.Z);
            my_float face_Z = Z - delta_u;

            delta_u = slope_limiter(prev.gamma,gamma,next.gamma);
            my_float face_gamma = gamma - delta_u;

            delta_u = slope_limiter(prev.t_cool,t_cool,next.t_cool);
            my_float face_t_cool = t_cool - delta_u;

            delta_u = slope_limiter(prev.meanweight,meanweight,next.meanweight);
            my_float face_mu = meanweight - delta_u;

            Cell result = Cell(data,face_P,face_T,face_gamma,face_Z,r-0.5*dr,r-0.5*dr,system_E);
            result.update_tcool(face_t_cool);
            result.update_meanweight(face_mu);

            return result;
        }
        
        Cell right_interface(Cell prev, Cell next){ //next is needed
            //In literature right interface is EAST (E)
            my_float data[NUM_EQUATIONS];
            my_float delta_u;
            for(int elem=0;elem<NUM_EQUATIONS;++elem){
                delta_u = slope_limiter(prev[elem],y[elem],next[elem]);
                data[elem] = y[elem] + delta_u; //Notice the sign. This is the change with 'left_interface'
            }
            //Thermodynamic variables must be updated accordingly
            delta_u = slope_limiter(prev.P,P,next.P);
            my_float face_P = P + delta_u;

            delta_u = slope_limiter(prev.T,T,next.T);
            my_float face_T = T + delta_u;

            delta_u = slope_limiter(prev.Z,Z,next.Z);
            my_float face_Z = Z + delta_u;

            delta_u = slope_limiter(prev.gamma,gamma,next.gamma);
            my_float face_gamma = gamma + delta_u;

            delta_u = slope_limiter(prev.t_cool,t_cool,next.t_cool);
            my_float face_t_cool = t_cool + delta_u;

            delta_u = slope_limiter(prev.meanweight,meanweight,next.meanweight);
            my_float face_mu = meanweight + delta_u;

            Cell result = Cell(data,face_P,face_T,face_gamma,face_Z,r+0.5*dr,r+0.5*dr,system_E);
            result.update_tcool(face_t_cool);
            result.update_meanweight(face_mu);

            return result;
        }

    public:
        ///MEMORY LOCATIONS OF NEIGHBOUR CELLS
        Cell* next;
        Cell* prev;
        /*
            Constructors does NOT link cells, you have to manually do that in the code.
            It's just writting 'current_cell.next = &next_cell.
            If you are not using pointers for the same cells!

            An unused cell is one that have both pointers as null. Remember that!
        */

        ///CONSTRUCTORS
        Cell(){
            next = nullptr;
            prev = nullptr;
        }
        Cell(const Cell& my_copy){ //Copy constructor
            for(int elem=0; elem<NUM_EQUATIONS; ++elem){
                y[elem] = my_copy.y[elem];
            }
            system_E = my_copy.system_E;

            P = my_copy.P;
            T = my_copy.T;
            gamma = my_copy.gamma;
            Z = my_copy.Z;
            t_cool = my_copy.t_cool;
            meanweight = my_copy.meanweight;

            r  = my_copy.r;
            dr = my_copy.dr;

            next = my_copy.next;
            prev = my_copy.prev;
        }
        
        //This is the main constructor, in which you only give data of mass, energy and momentum density; and entropy as well
        Cell(my_float data[NUM_EQUATIONS], my_float rmin, my_float rmax,bool check = true){
            //Init the equations
            for(int elem=0;elem<NUM_EQUATIONS;++elem){
                y[elem] = data[elem];
            }

            //Safety net: avoid a hard zero for energy density (and entropy)
            if(std::abs(y[1]) < SOFT_ZERO){
                y[1] = SOFT_ZERO;
                y[3] = entropy(y[0],P); //The equivalent for entropy
            }

            system_E = check;
            //To avoid funky issues, initialize metallicity first
            this->Z = INIT_METALLICITY;

            //Location of the cell (moved before 'update_thermodynamics' just for testing, order is irrelevant
            r  = (rmax+rmin)*0.5;
            dr = std::abs(rmax-rmin); //No std => result is an int

            //Init other variables
            this->update_thermodynamics();


            next = nullptr;
            prev = nullptr;

        }
        //This is an expanded constructor, in which you give ALL data minus cooling time and heating function! (that should be initialized with 'update_thermodynamics' routine)
        Cell(my_float data[NUM_EQUATIONS], my_float my_P, my_float my_T, my_float my_gamma, my_float my_Z, my_float rmin, my_float rmax, bool check = true){
            //Init the equations
            for(int elem=0;elem<NUM_EQUATIONS;++elem){
                y[elem] = data[elem];
            }

            //Init other variables
            this->P = my_P;
            this->T = my_T;
            this->gamma = my_gamma;
            this->Z = my_Z;

            //t_cool is defaulted to zero, it's only used to get a timestep in the main loop
            //And it's updated later
            this->t_cool = -INFINITY; //Yes, you read right! by default this is negative infinity! (physical meaning = no cooling)
            ///^^^WARNING HERE, cooling appear when computing sources and not fluxes. Fluxes uses this prescription and this value is unused.

            //Safety net: avoid a hard zero for energy density (and entropy)
            if(std::abs(y[1]) < SOFT_ZERO){
                y[1] = SOFT_ZERO;
                y[3] = entropy(y[0],P); //The equivalent for entropy
            }
            this->meanweight = compute_meanweight(y[0],P,T); //Is this a HUGE liability?

            r  = (rmax+rmin)*0.5;
            dr = std::abs(rmax-rmin); //No std => result is an int
            system_E = check;

            next = nullptr;
            prev = nullptr;
        }

        ///UPDATE SOME VALUE

        void update_tcool(my_float my_tcool){
            t_cool = my_tcool;
        }
        void update_meanweight(my_float my_meanweight){
            meanweight = my_meanweight;
        }
        void update_metallicity(my_float my_Z){
            Z = my_Z;
        }

        ///SHOW ME SOME VALUE
        my_float eq(int index){
            if(index >= NUM_EQUATIONS || index < 0){
                throw std::out_of_range("index out of range");
            }else{
                return y[index];
            }
        }
        //Operator overloading
        my_float& operator[](const int index){//HEADACHE
            if(index >= NUM_EQUATIONS || index < 0){
                throw std::out_of_range("index out of range");
            }else{
                return y[index];
            }
        }
        
        void change_system(){ system_E = !system_E; }

        my_float get_pressure(){return P;}
        my_float get_temperature(){return T;}
        my_float get_metallicity(){return Z;}
        my_float get_adiabatic_const(){return gamma;}
        my_float get_meanweight(){return meanweight;}

        my_float get_cooling_time(){
            //Carerful with this one! this return the absolute value of t_cool
            //Use cooling function for the sign instead
            return t_cool;
        }

        my_float location(){return r;}
        my_float width(){return dr;}

        bool get_system(){return system_E;}

        my_float get_velocity(){return y[2]/y[0];}
        my_float get_magnetosonic(){return magnetosonic_speed(y[0],P,gamma);}

        my_float get_number_density(){return number_density(y[0],meanweight);}

        my_float get_internal_energy(){return internal_energy(P,gamma);}
        my_float get_specific_internal_energy(){return specific_internal_energy(P,y[0],gamma);}

        my_float get_cooling_function(){
            //Techincally, is this evaluation.
            #if ACTIVATE_COOLING == 1
                return internal_energy(P,gamma) / ( number_density(y[0],meanweight)*number_density(y[0],meanweight) * t_cool );
                //If t_cool = -inf, you should get a zero. But I saw NaN as well, and that's a problem with post-processing
            #else
                return 0.0;
            #endif
        }
        my_float get_heating_function();

        my_float get_volume(){
            //Get cell volume
            my_float r_min = this->r - this->dr;
            my_float r_max = this->r + this->dr;
            #if SYMMETRY == 2
                //Spherical coordinates. You get a volume
                return (4.*PI/3.)*(r_max*r_max*r_max - r_min*r_min*r_min);
            #elif SYMMETRY == 1
                //Cylindrical coordinates. You get an area
                return PI*(r_max*r_max - r_min*r_min);
            #else
                //Cartesian coordinates. You get a length
                return this->dr;
            #endif
        }

        //Due to precedence issues with () and -> . I've been forced to write these two functions
        Cell* get_next(){ return next; }
        Cell* get_prev(){ return prev; }

        ///MOVE CELL
        Cell move_location(my_float new_r){
            //Moves the Cell from 'this->r' to 'new_r'.
            //Remaining magnitudes are conserved
            return Cell(this->y,this->P,this->T,this->gamma,this->Z,new_r-0.5*dr,new_r+0.5*dr,this->system_E);
            //Note that this->dr is the same as writting dr alone
        }

        void Link_next(Cell& next_cell){
            this->next = &next_cell;
            next_cell.prev = this;
        }
        void Link_prev(Cell& prev_cell){
            this->prev = &prev_cell;
            prev_cell.next = this;
        }
        bool is_unused(){
            if( next == nullptr && prev == nullptr ){ return true; }
            else{ return false; }
        }

        ///AMR RELATED FUNCTIONS
        int level_of_refinement(){return LoR;}

        int compare_LoRs(Cell* second){
            //Positive values means that *second has a higher LoR than this
            //Negative values means the opposite

            if(second != nullptr){
                return (second->LoR - this->LoR);
            }else{
                return 0;
            }
        }

        void needs_merge(Cell*,bool);
        void needs_split(Cell*,bool);
        void near_shock(); //Also catches shocks
        void high_gradient();
        void check_refinement(Cell*,bool);

        bool unacceptable_value(){
            for(int elem=0;elem<NUM_EQUATIONS;++elem){
                //A value is unacceptable if is negative (except momentum)
                if(y[elem] < 0. && elem != 2){ return true; }
            }
            return false;
        }

        void split(Cell*,const int);
        void merge(Cell*,Cell*,Cell*,const int);
        //void restore(Cell*,Cell*,const int,const int);

        void split(Cell*,const int,INTERPOLATOR*);
        void merge(Cell*,INTERPOLATOR*);

        void prepare_for_refinement(){
            this->flag_for_merge = false;
            this->flag_for_split = false;
            this->flag_for_shock = 0;
            this->flag_for_gradient = 0;
        }
        bool should_be_split(){ return this->flag_for_split; }
        bool should_be_merged(){ return this->flag_for_merge; }

        int shocked_cell(){ return this->flag_for_shock; }
        int high_gradient_cell(){ return this->flag_for_gradient; }

        void force_split(){ this->flag_for_split = true; }
        void deny_merge(){ this->flag_for_merge = false; }

        ///SOFT COPY CELL
        //A.k.a keep the pointers!
        void soft_update(Cell my_copy){
            //Same as soft_copy, except for the return type.
            //Update the values of the cell, but preserve the pointers
            for(int elem=0; elem<NUM_EQUATIONS; ++elem){
                y[elem] = my_copy.y[elem];
            }
            system_E = my_copy.system_E;

            P = my_copy.P;
            T = my_copy.T;
            gamma = my_copy.gamma;
            Z = my_copy.Z;
            t_cool = my_copy.t_cool;
            meanweight = my_copy.meanweight;
            //heat_e = my_copy.heat_e;

            r  = my_copy.r;
            dr = my_copy.dr;
        }

        ///GENERATE INTERFACE FLUXES
        //I would like to write the function here, but compiler complains.
        //Code is after class Flux code, inside 'fluxclass.h'
        Flux left_flux(Cell&, Cell&);
        Flux right_flux(Cell&, Cell&);

        ///SIMULATION PUBLIC FUNCTIONS
        Cell predict(Flux*,my_float);

        Cell correct(Cell,Flux*,my_float);

};


#endif
