#ifndef THERMODYNAMICS_CPP
#define THERMODYNAMICS_CPP

//As title says, solves thermodynamics
bool thermodynamics_negative_pressure_warning_given = false;

///NEEDED TABLES
Bilinear EoS_table;
void init_EoS_table(){
    EoS_table = load_table(EOS_FILE);
}

#if ACTIVATE_COOLING == 1 
    Bilinear Cooling_table;
    void init_cooling_table(){
        Cooling_table = load_table(COOLING_FILE);
    }
#endif

#if ACTIVATE_HEATING == 1
    Bilinear heating_table;    
    void init_heating_table(){
        heating_table = load_table(HEATING_FILE);
    }
#endif

bool use_dual_energy(Cell* curr){
    //Dual energy formalism uses entropy instead of energy to track values.
    //Using energy allows energy conservation, and using entropy allows accurate temperatures.

    ///CONDITION 1: Thermal energy is too low.
    my_float v_i = curr->eq(2)/curr->eq(0);
    my_float curr_P = pressure(curr->eq(0),curr->eq(1), v_i, ADIAB_CONST);
    if(curr_P == FAILURE){
        //Check if we are not in a shock (entropy is NOT conserved there)
        my_float v_r; //Velocity_right (interface)
        my_float r_r; //Radius_right
        if(curr->next != nullptr){
            v_r = (curr->next->eq(2) )/(curr->next->eq(0) );
            r_r = curr->next->location();
        }else{
            v_r = v_i;
            r_r = curr->location();
        }
        my_float v_l; //Velocity_left (interface)
        my_float r_l; //Radius_left
        if(curr->prev != nullptr){
            v_l = (curr->prev->eq(2) )/(curr->prev->eq(0) );
            r_l = curr->prev->location();
        }else{
            v_l = v_i;
            r_l = curr->location();
        }
        
        my_float rq_r = pow(r_r,SYMMETRY);
        my_float rq_l = pow(r_l,SYMMETRY); //rq = r^q = r^SYMMETRY (0,1,2 for Cartesian,Cylindrical,Spherical)
        my_float condition = rq_r*v_r - rq_l*v_l;

        if(condition >= 0.){ 
            if(thermodynamics_negative_pressure_warning_given == false){
                std::cout<<"Warning: Entropy conservation used instead of energy's! This may give unexpected results."<<std::endl;
                thermodynamics_negative_pressure_warning_given = true;
            }
            return true;
        }
        else{ return false; }
    }
    else{ return false; }
}
/// MAIN ROUTINE
void Cell::update_thermodynamics(){ 

    ///Step 1: Get pressure with equation variables (Dual energy is checked here!)
    if(!use_dual_energy(this)){ //If entropy system is NOT used!
        this->system_E = true;
        this->P = pressure(this->y[0],this->y[1],this->y[2]/this->y[0],ADIAB_CONST);
    }else{
        this->system_E = false;
        this->P = pressure(this->y[0],this->y[3],ADIAB_CONST);
    }
    ///Step 2: Read the EoS table to get the temperature and mean weight
    #if EOS_INTERPOLATION == 1
        this->T = EoS_table.get_value(this->y[0],this->P);
    #else
        this->T = EoS_table.get_log(this->y[0],this->P);
    #endif
    
    this->meanweight = compute_meanweight(this->y[0],this->P,this->T);
    
    ///Step 3: Get the cooling time and thus the cooling
        //Note that, in my code, cooling rate has negative sign (in the table is positive sign. This is to avoid computing a logarithm of a negative value).
    this->gamma = ADIAB_CONST;
    #if ACTIVATE_COOLING == 1
        if(this->T < MAX_COOLING_TEMPERATURE){
            #if COOLING_INTERPOLATION == 1 
                my_float dot_u = -Cooling_table.get_value(this->y[0],this->P);
            #else
                my_float dot_u = -Cooling_table.get_log(this->y[0],this->P);
            #endif
                this->t_cool = this->get_internal_energy() / dot_u;
            }else{
                this->t_cool = -1e100; //If T>MAX_COOLING_TEMPERATURE I remove cooling to avoid overcooling when the shock is forming at the very first timesteps.
                /*For some reason, writing -INFINITY here gives a NaN at some point that is not here (see below, I don't get a 'Bad Cooling' error)
                 Looks like that some kind of "soft infinity" (i.e.: very high number) fixes this issue*/
            }
    #else
        this->t_cool = - INFINITY;//Yes, you read right! by default this is negative infinity! (physical meaning = no cooling)
    #endif
    
    //thow an error if there's something wrong
    if(this->P != this->P || this->t_cool != this->t_cool){ // || this->meanweight < 0 || this->meanweight > 2){ //You got an error!
        std::cout<<this->location()<<" "<<this->y[0]<<" "<<this->y[1]<<" "<<this->y[2]<<" "<<this->y[3]<<std::endl;
        //Note to future self, if the output here is a nan in both y[1] and y[3], check the source terms in correct function (timeEvolution.cpp)
        if(this->P != this->P){
            throw std::runtime_error("Bad Pressure!");
        }
        else{
            throw std::runtime_error("Bad Cooling");
        }
    }
}

my_float Cell::get_heating_function(){
    //Get the value of the heating function (in erg/cm3/s) by reading a table
    #if ACTIVATE_HEATING == 1
        //my_float number_density = this->get_number_density(); //To avoid two function calls
        if(this->T < MAX_HEATING_TEMPERATURE){
            #if HEATING_INTERPOLATION == 1
                return this->get_number_density() * heating_table.get_value(this->y[0],this->P); //erg/cm3/s
            #else
                return this->get_number_density() * heating_table.get_log(this->y[0],this->P); //erg/cm3/s
            #endif
        }else{
            return 0.0;
        }
    #else
        return 0.0;
    #endif
}
#endif
