#ifndef TIME_EVOLUTION_CPP
#define TIME_EVOLUTION_CPP

/****

    Here is stored the predictor-corrector scheme used for the equations, and the choice of the correct timestep

****/

/// SET SIMULATION TIMESTEP

my_float basic_stability(Cell my_cell){
    /*The timestep is stable if
        v*(dt/dx) < CFL_NUMBER
        CFL_NUMBER is 0.5 or less
    */

    double v = std::abs(my_cell.get_velocity());
    double c = my_cell.get_magnetosonic(); //Sound speed
    double a = (v > c) ? v : c ; //If true, take v. If false take c.
    double result = std::abs(my_cell.width()/a); //= dx/v

    return result*CFL_NUMBER;
}

my_float set_timestep(Cell my_cell,my_float max_dt = INFINITY){
    //This routine compares the maximum timestep of each stability function (below)
    //and compares with max_dt.
    //Returns the minimum timestep

    my_float next_dt = max_dt;                      //Custom timestep
    my_float ba_dt   = basic_stability(my_cell);    //Hydrodynamic timestep
    my_float cool_dt = std::abs(my_cell.get_cooling_time());  //Cooling timestep
    cool_dt *= CFL_NUMBER2; //Note that basic_stability reduces its timestep by the same amount

    //Take the minimum of all timesteps
    next_dt = (next_dt > ba_dt) ? ba_dt : next_dt;
    next_dt = (next_dt > cool_dt) ? cool_dt : next_dt;
    
    if(next_dt < SOFT_ZERO){
        std::cout<<"Gotcha!"<<std::endl;
        std::cout<<ba_dt<<' '<<cool_dt<<' '<<next_dt<<std::endl;
        throw std::runtime_error("Test stop!");
    }
    
    return next_dt;
}


/// PREDICTOR-CORRECTOR SCHEMES

Cell Cell::predict(Flux* fluxes, my_float dt){
    /**
        Allow me to explain each argument

        '*fluxes' is an array of four Cells: [prev2, prev, next, next2]. Caller is not needed because is '*this'
        'dt' is the current simulation timestep.

        For this function, the return object is the time evolution of *this, following an Euler scheme without sources:
        result = (*this) + dt*this->get_derivative(...)

        Previous line is this function in a nutshell, but devil is in details since there are lots of safety nets...
    **/

    ///-----------------
    //Compute the derivative of *this
    my_float* derivatives = this->get_derivative(fluxes[0], fluxes[1], fluxes[2], fluxes[3]);

    ///-----------------
    //Compute hydrodynamic magnitudes

    my_float data[NUM_EQUATIONS];

    for(int elem=0;elem<NUM_EQUATIONS;++elem){
        data[elem] = this->eq(elem) + dt*derivatives[elem];
    }

    delete[] derivatives; //get_derivative generates a dynamic array, so we have to delete it.

    ///-----------------
    //Compute thermodynamic magnitudes

    Cell result = Cell(data,r-0.5*dr,r+0.5*dr,this->system_E); //Adding my_fields and my_units as arguments, thermodynamics and dual energy formalism are updated as well
    //result will inherit caller system_E. However, during construction of 'result' this may change, specially if system_E was true.
    //result WON'T inherit pointers, though!!

    ///-----------------
    //resync energy and entropy

    result.resync();//my_fields,my_units);

    ///-----------------
    //give result
    return result;
}


Cell Cell::correct(Cell original, Flux* fluxes, my_float dt){
    /**
        Allow me to explain each argument

        'original' is the object from *this comes from. i.e.: *this = original.predict(...)
        '*fluxes' is an array of four Cells: [prev2, prev, next, next2]. However, these comes from the predictor values.
        'dt' is the current simulation timestep.

        For this function, we make two operations:
        -Compute the average of *this and original, and use the result to compute source terms
        -Compute:
            average + 0.5*dt*this->get_derivative(...) + dt*average.get_sources()

        Again, devil is here as well
    **/

    ///-----------------
    //Compute the average

    my_float av_data[NUM_EQUATIONS];

    //Cell average;

    for(int elem=0;elem<NUM_EQUATIONS;++elem){
        av_data[elem]  = 0.5 * (original.eq(elem) + this->eq(elem));
    }
    //Thermodynamics should be updated with the same formula
    my_float av_P = 0.5 * (this->P + original.P);
    my_float av_T = 0.5 * (this->T + original.T);
    my_float av_gamma = 0.5 * (this->gamma + original.gamma);
    my_float av_Z = 0.5 * (this->Z + original.Z);

    bool av_system = this->system_E & original.system_E; //=logic AND

    //Call the extended constructor
    Cell average = Cell(av_data,av_P,av_T,av_gamma,av_Z,r-0.5*dr,r+0.5*dr,av_system);
    
    //Update t_cool to avoid a default t_cool = -inf
    //Do the same with heating energy, as it is not created in the constructor
    average.t_cool =  0.5 * (this->t_cool + original.t_cool);

    ///-----------------
    //Get sources and derivative
    my_float* derivatives = this->get_derivative(fluxes[0], fluxes[1], fluxes[2], fluxes[3]); //derivatives is a dynamic array
    my_float* sources = average.get_sources(); //Another dynamic array

    ///-----------------
    //Compute hydrodynamic magnitudes

    my_float data[NUM_EQUATIONS];
    //Cell result;

    for(int elem=0;elem<NUM_EQUATIONS;++elem){
        data[elem] = average.eq(elem) + 0.5*dt*derivatives[elem] + dt*sources[elem];
    }

    delete[] derivatives; //No needed anymore
    delete[] sources;

    ///-----------------
    //Compute thermodynamic magnitudes

    Cell result = Cell(data,r-0.5*dr,r+0.5*dr,this->system_E); //Adding my_fields and my_units as arguments, thermodynamics and dual energy formalism are updated as well
    //result will inherit caller system_E. However, during construction of 'result' this may change, specially if system_E was true.
    //result WON'T inherit pointers, though!!

    ///-----------------
    //resync energy and entropy

    result.resync();//my_fields,my_units);

    ///-----------------
    //give result

    return result;

}

#endif
