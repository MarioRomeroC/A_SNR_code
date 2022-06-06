#ifndef GAS_PHYSICS_CPP
#define GAS_PHYSICS_CPP

//Gas physics related functions are stored here
//Unlike thermodynamics, these functions are the easy ones found in textbooks (e.g.: Pressure of an ideal gas)

/// MEANWEIGHT

my_float compute_meanweight(my_float rho, my_float P, my_float T){

    my_float num = BOLTZMANN_CTE*T*rho;
    my_float den = HYDROGEN_MASS*P;

    return num/den;
}


/// PRESSURE

my_float pressure(my_float rho, my_float e, my_float v, my_float gamma, my_float B){
    //Compute pressure according to an ideal gas equation.
    //Arguments are mass and energy density, velocity, adiabatic constant, and magnetic field (default value is zero)

    my_float kin_e = 0.5*rho*v*v;
    my_float mag_e = 0.5*B*B/VACUUM_PERM;
    my_float thm_e = e - kin_e - mag_e;

    //Dual energy check (see 'use_dual_energy' in 'thermodynamics.cpp')
    if(thm_e / e > ETA_E){ return (gamma - 1.)*thm_e; }
    else{ return FAILURE; }
}

my_float pressure(my_float rho, my_float S, my_float gamma){
    return S*pow(rho,gamma -1.);
}

/// TEMPERATURE

my_float temperature(my_float rho, my_float P, my_float mu){
    //Following an ideal gas equation of state

    my_float n = number_density(rho,mu);
    my_float result = P/(n*BOLTZMANN_CTE);

    if( result < 0.0 ){ return FAILURE;  }
    else{ return result; }
}

/// NUMBER DENSITY

my_float number_density(my_float rho,my_float mu){
    return rho / (mu*HYDROGEN_MASS);
}

/// ENERGY
//Units of erg/cm^3 unless 'specific' is in the name of the function. In that case, units are erg/g

my_float internal_energy(my_float P, my_float gamma){
    return P / (gamma - 1.);
}

my_float specific_internal_energy(my_float P, my_float rho, my_float gamma){
    return internal_energy(P,gamma)/rho;
}

my_float energy_density(my_float rho, my_float v, my_float P, my_float gamma = ADIAB_CONST){

    my_float thm_e = P/(gamma - 1.);
    my_float kin_e = 0.5*rho*v*v;

    return kin_e + thm_e;
}

/// ENTROPY

my_float entropy(my_float rho, my_float P, my_float gamma = ADIAB_CONST){
    return P/pow(rho,gamma-1.);
}

/// OTHERS

my_float magnetosonic_speed(my_float rho, my_float P, my_float gamma, my_float B){
    //Compute SQUARED sound speed
    my_float cs2 = (gamma * P) / rho;
    //Compute SQUARED Alfvén speed
    my_float ca2 = B*B / (VACUUM_PERM*rho);

    //now give result, in velocity units
    return std::sqrt(cs2+ca2);
}

#endif
