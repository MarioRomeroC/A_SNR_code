#ifndef HYDRODYNAMICS_CPP
#define HYDRODYNAMICS_CPP

//As title says, solves hydrodynamics
//All functions related to that are included here


/// MAIN ROUTINES

my_float* Cell::get_derivative(Flux flux_p_E,Flux flux_c_W,Flux flux_c_E,Flux flux_n_W){
    ///---------------------------------------------
    //Next step is to get local velocities, labeled as a^+ or a^-
    //You need a pair for each interface (total = two a^+ and two a^-)

    //West max local speed
    my_float a_max_W = max_value(flux_c_W.local_speed_max(),flux_p_E.local_speed_max(), 0.);

    //West min local speed
    my_float a_min_W = min_value(flux_c_W.local_speed_min(),flux_p_E.local_speed_min(), 0.);

    //East max local speed
    my_float a_max_E = max_value(flux_c_E.local_speed_max(),flux_n_W.local_speed_max(), 0.);

    //East min local speed
    my_float a_min_E = min_value(flux_c_E.local_speed_min(),flux_n_W.local_speed_min(), 0.);


    ///---------------------------------------------
    //Now get KNP fluxes with all the values we have computed
    //You need one KNP flux for each interface

    my_float F_E[NUM_EQUATIONS];
    my_float F_W[NUM_EQUATIONS];

    my_float* data = new my_float[NUM_EQUATIONS];

    //For using Cartesian, cylindrical and spherical coordinates, I'm using method 2 of Wang & Johnsen 2017
    my_float r_E = flux_c_E.location();//curr_E.location();
    my_float r_W = flux_c_W.location();//curr_W.location();
    my_float super_r_E = std::pow(r_E,SYMMETRY);
    my_float super_r_W = std::pow(r_W,SYMMETRY);

    for(int elem=0;elem<NUM_EQUATIONS;++elem){

        //Get West KNP flux (You will find East KNP flux, I had problems finding this version...)
        F_W[elem] = a_max_W * flux_p_E.get_flux(elem) - a_min_W * flux_c_W.get_flux(elem)
            + a_max_W * a_min_W * ( flux_c_W.eq(elem) - flux_p_E.eq(elem) ); //Numerator
        F_W[elem] = F_W[elem] / (a_max_W - a_min_W); //Denominator

        //Get East KNP flux
        F_E[elem] = a_max_E * flux_c_E.get_flux(elem) - a_min_E * flux_n_W.get_flux(elem)
            + a_max_E * a_min_E * ( flux_n_W.eq(elem) - flux_c_E.eq(elem)) ; //Numerator
        F_E[elem] = F_E[elem] / (a_max_E - a_min_E) ; //Denominator

        //Get the derivative for 'elem' equation
        //SYMMETRY = {0,1,2} for {cartesian, cylindrical, spherical}
        data[elem] = -(SYMMETRY + 1.)*(super_r_E * F_E[elem] - super_r_W * F_W[elem])/(super_r_E*r_E - super_r_W*r_W);

    }
    ///Add a geometric source term in momentum equation, or method 2 doesn't work for cylindrical/spherical coordinates.
        data[2] += ( (SYMMETRY+1.)*(super_r_E*flux_c_E.P - super_r_W*flux_c_W.P)/(super_r_E*r_E - super_r_W*r_W) ) - ( (flux_c_E.P - flux_c_W.P)/(r_E - r_W) );
        //Note that this source is zero if we remain in cartesian coordinates (SYMMETRY = 0).

    return data;
}

//Get source
my_float* Cell::get_sources(){
    ///Self-explanatory: get the source term of Euler equations
    //Only source not taken into consideration is geometric source, which is added in get_derivative

    ///----------------
    //Construct *this Flux

    Source my_sources = Source(*this);

    ///----------------
    //Create the source

    my_float* data = new my_float[NUM_EQUATIONS];

    for(int elem=0;elem<NUM_EQUATIONS;++elem){
        data[elem] = my_sources.get_source(elem);
    }

    return data;
}

/// CELL RECONSTRUCTION

//Max-Min of three values
my_float max_value(my_float value1, my_float value2, my_float value3){

    my_float tmp = (value1>value2) ? value1 : value2;
    return (tmp>value3) ? tmp : value3;
}

my_float min_value(my_float value1, my_float value2, my_float value3){

    my_float tmp = (value1<value2) ? value1 : value2;
    return (tmp<value3) ? tmp : value3;
}

//Limiters
my_float slope_limiter(my_float prev, my_float curr, my_float next){
    //Computes a correction based on van Leer 1997
    //This can be used for left and right interface of current cell

    my_float value = (next - curr) * (curr - prev);
    value = (value>0) ? value : 0.;
    my_float result = value / (next-prev);
    if( result != result){return 0.;} //If result = NaN, give me a zero. (Happens if you have constant density across a fluid section, for example)
    else{return result;}
}

#endif
