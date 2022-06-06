#ifndef SET_AMBIENT_CPP
#define SET_AMBIENT_CPP

///Auxiliary functions in this file
//Note that these are written as such to try different initial conditions for explosions
#if MASS_ON == 1
    my_float structure_function(my_float r){
        //This is \Delta \rho (r), important!
        //Do not worry about normalization, code will do that for you
        return 1.;
    }
#else
    my_float structure_function(my_float w){
        //Kept because it worked really well.
        return (1./std::sqrt(PI))*std::exp(-w*w); //Just a Gaussian
    }
#endif

my_float w_profile(my_float i, my_float n){
    //This defines w. In some cases, x=r and X=R_shock
    //In other cases, x=v and X=v_shock;
    //For a general use, x=variable and X=Maximum value. w must remain dimensionless!
    return 1.;
}

void normalize(my_float* my_function, int N, my_float const_num=1.){
    my_float norm = 0.0;
    for(int i=0;i<N;++i){
        norm += my_function[i];
    }
    for(int i=0;i<N;++i){
        my_function[i] *= (const_num/norm);
    }
}

#if MASS_ON == 1
my_float create_central_explosion(Cell* grid, const int N_cells){
    //Returns starting time, but the important part is to get IC profiles right
    //In order to work properly, all cells must be at the same refinement level (you only have to care about that at starting time)
    //Note that constant velocity profile is a special case.

    //(1.1) Define some constants related with define.h
    const my_float shock_radius = SHOCK_RADIUS*PARSEC; //In define.h it was in pc, here is now in code units (cm)

    const my_float ejected_mass = EJECTED_MASS*SOLAR_MASS; //In define.h it was in solar masses, here is now in code units (g)
    const my_float ejected_energy = EXPLOSION_ENERGY;

    const my_float explosion_meanweight = 0.62; //That number is the value for a fully ionized gas

    //(1.2) Define other constants related with this simulation
    const my_float simulation_radius = grid[N_cells-1].location();

    const int n = std::round(N_cells*(shock_radius/simulation_radius));

    //(2.1) Create needed arrays.
    my_float w[n]; //w profile (=velocity profile)
    my_float V[n]; //Cell volumes
    my_float f[n]; //Structure functions
    my_float normalization = 0.;

    for(int i=0; i<n; ++i){
        my_float ri = grid[i].location();
        my_float dr = grid[i].width();

        w[i] = w_profile(my_float(i),my_float(n));
        V[i] = (4.*PI/3.) * ( std::pow(ri+0.5*dr,3.) - std::pow(ri-0.5*dr,3.) );
        f[i] = structure_function(grid[i].location());

        normalization += f[i]*V[i];
    }
    //(2.2) Repeat again the loop to get normalized f[i] [ Sum(f[i])=1. ]
    my_float energy_sum = 0.; //This will be the Sum(f[i]*w[i]*w[i]). This is needed for shock velocity and energy injection
                                //It's necessary if radial velocity profile is not constant
    for(int i=0; i<n; ++i){
        f[i] *= ejected_mass/normalization; //No phytonic way to do that. Too bad!
        energy_sum += f[i]*w[i]*w[i]*V[i];
    }


    //(3.1) Get shock velocity and starting time
    my_float v_shock = std::sqrt( 2.*ejected_energy / ( energy_sum )  );
    my_float t_0 = shock_radius / v_shock; //This is the return value.

    //(3.2) Get Mass and energy injection
    my_float explosion_data[NUM_EQUATIONS];
    my_float explosion_P;
    my_float explosion_T;
    bool explosion_system_E = true;

    for(int i=0;i<n;++i){
        my_float ri = grid[i].location();
        my_float dr = grid[i].width();
        my_float vel = w[i]*v_shock; //Velocity profile.

        my_float added_mass   = f[i];
        my_float added_energy = 0.5*f[i]*vel*vel;

        //Hydrodynamic variables
        explosion_data[0] = grid[i].eq(0) + added_mass;

        explosion_data[1] = grid[i].eq(1) + added_energy;

        explosion_data[2] = grid[i].eq(2) + added_mass * vel;
        explosion_data[3] = explosion_data[1]*(ADIAB_CONST - 1.)/std::pow(explosion_data[0], ADIAB_CONST - 1.);

        //Thermodynamic variables
        explosion_P = pressure(explosion_data[0],explosion_data[1],explosion_data[2]/explosion_data[0],grid[i].get_adiabatic_const());
        if(explosion_P == FAILURE){ //Use dual energy formalism
            explosion_P = pressure(explosion_data[0],explosion_data[3],grid[i].get_adiabatic_const());
            explosion_system_E = false;
        }
        explosion_T = explosion_meanweight*(explosion_P/explosion_data[0])*(HYDROGEN_MASS/BOLTZMANN_CTE);

        //Recreate the cell again
        Cell explosion_cell = Cell(explosion_data,ri-0.5*dr,ri+0.5*dr);

        grid[i].soft_update(explosion_cell); //Soft copy preserves the pointers
    }

    return t_0;
}

#else
my_float create_central_explosion(Cell* grid, const int N_cells, my_float my_explosion, my_float my_ejected_mass){
    //Following Doumler & Knebe 10, follow their prescription for a Sedov blastwave.
    //You only need a few cells

    //Note: If you have to reproduce the same IC with higher/lower resolution. Multiply/divide the next two const int.
    //For example, if you want double resolution, multiply by 2 'stop_at_radius' and 'shock_position'
    const int stop_at_radius = 4*100; //=3.5*dx
    const int shock_position = 5;
    const my_float ejected_mass = SOLAR_MASS*my_ejected_mass;
    const my_float shock_velocity = std::sqrt( 2.*my_explosion/ejected_mass );
    const my_float shock_radius = shock_position*grid[0].width();
    
    my_float t0 = shock_radius/shock_velocity;

    //Normalize structure function:
    my_float* f = new my_float[stop_at_radius];
    my_float norm = 0.;


    for(int i=0;i<stop_at_radius;++i){
        f[i] = structure_function(grid[i].location()/shock_radius);
        norm += f[i];
    }
    for(int i=0;i<stop_at_radius;++i){
        f[i] /= norm;
    }

    //Do things
    my_float explosion_data[NUM_EQUATIONS];

    for(int i=0;i < stop_at_radius; ++i){
        my_float xi = grid[i].location() ;
        my_float dx = grid[i].width() ;
        my_float Vi = (4./3.)*PI*( std::pow(xi+0.5*dx,3.) - std::pow(xi-0.5*dx,3.) );

        my_float injected_energy = (my_explosion/Vi)*f[i];

        explosion_data[0] = grid[i].eq(0);// + MASS_ON*(ejected_mass/Vi)*structure_function(wi);
        explosion_data[1] = grid[i].eq(1) + injected_energy;
        explosion_data[3] = explosion_data[1]*(ADIAB_CONST - 1.)/std::pow(explosion_data[0], ADIAB_CONST - 1.);
        explosion_data[2] = grid[i].eq(2);

        //Recreate the cell again
        
        Cell explosion_cell = Cell(explosion_data,xi-0.5*dx,xi+0.5*dx);

        grid[i].soft_update(explosion_cell); //Soft copy preserves the pointers
    }

    delete[] f;

    return t0;
}
#endif


my_float set_ambient(Cell* grid, string name, const int N_cells){
    //Will create your ambient. For the first step, cells are ordered (i.e grid[i] is guaranteed to be i-th cell, which is no longer true when refinement activates later)

    //Open file and read it.
    ifstream inFile;

    inFile.open(name);

    my_float ambient_data[NUM_EQUATIONS];
    my_float ejected_mass;
    my_float explosion_energy;
    my_float simulation_box;
    my_float ambient_Z;

    inFile>>ambient_data[0];
    inFile>>ambient_data[1];
    inFile>>ambient_data[2];

    ambient_data[3] = ambient_data[1]*(ADIAB_CONST - 1.)/std::pow(ambient_data[0], ADIAB_CONST - 1.);

    inFile>>ambient_Z;

    inFile>>simulation_box;

    #if CENTRAL_EXPLOSION == 1
        inFile>>explosion_energy;
        inFile>>ejected_mass;
    #endif // CENTRAL_EXPLOSION

    inFile.close();

    //Now fill grid
    my_float dx = simulation_box/N_cells;
    my_float xi;
    for(int i = 0; i<N_cells; ++i){
        xi = (my_float(i) + 0.5)*dx; //x_0 = 0

        grid[i] = Cell(ambient_data,xi-0.5*dx,xi+0.5*dx);
        if(i>0){ grid[i].Link_prev(grid[i-1]); }
    }

    //Now set extra things
    #if CENTRAL_EXPLOSION == 1
        #if MASS_ON ==1
        return create_central_explosion(grid,N_cells, my_fields,my_units);
        #else
        return create_central_explosion(grid, N_cells,explosion_energy,ejected_mass);
        #endif
    #else
        return 0.0;
    #endif
}


void create_map(){
    /*
        This function is used to create useful info, such as values for T, Lambda and mu for each pair (rho,P)

        It is NOT used for the simulation, it only stores useful info
    */
    
    const int N_rho = 18*20; //Before: 12*20
    const int N_P = 28*10; //+1;
    //Uncomment these to have ridiculous resolution (But the file is going to be REALLY big!)
    //const int N_rho = 18*200;
    //const int N_P = 28*100;
    
    const int N = N_rho * N_P; //Subdivisions

    //Create our limits
    const my_float max_rho = 10000.;// * HYDROGEN_MASS; //In amu/cm^3
    const my_float min_rho = 0.00001;// * HYDROGEN_MASS;
    const int max_rho_i = std::log10(max_rho);
    const int min_rho_i = std::log10(min_rho);

    const my_float max_P = 1.e13;// * BOLTZMANN_CTE; //In K/cm^3
    const my_float min_P = 1.e-1;// * BOLTZMANN_CTE;
    const int max_P_j = std::log10(max_P);
    const int min_P_j = std::log10(min_P);

    //Some dummy variables [remember: dumb values because we only want a map, not a simulation. it is expected to fail at first timestep]
    const my_float  x = -1.0;
    const my_float dx = -1.0;

    //Create the main loop
    const my_float di = (max_rho_i - min_rho_i)/my_float(N_rho);
    const my_float dj = (max_P_j - min_P_j)/my_float(N_P);
    
    ofstream mapfile;
    string path_to_file = OUTPUT_PATH;
    string name1 = "Map_Z";
    string name2 = tostring(INIT_METALLICITY);
    string name3 = ".dat";
    mapfile.open(path_to_file+name1+name2+name3);
    
    mapfile<<" density | "<<"pressure | "<<"temperature | "<<"number density | "<<"cooling function | "<<"mean weight | "<<"cooling time | ";
    mapfile<<"heating rate "<<endl;
    mapfile<<std::setprecision(NUMBER_PRECISION);
    
    my_float heating_value = 0.0;
    my_float net_cooling;

    for(my_float i=min_rho_i; i<=max_rho_i; i+=di){
        for(my_float j=min_P_j; j<=max_P_j; j+=dj){
            //Create the dummy cell
            my_float map_data[NUM_EQUATIONS];
            my_float rho = HYDROGEN_MASS * std::pow(10.,i);
            my_float P = BOLTZMANN_CTE * std::pow(10.,j);

            map_data[0] = rho; //Density, this is easy
            map_data[1] = P / (ADIAB_CONST - 1.); // Energy density
            map_data[2] = 0. ; //No velocity so the line above is true
            map_data[3] = map_data[1]*(ADIAB_CONST - 1.)/std::pow(map_data[0], ADIAB_CONST - 1.);

            Cell dummy_cell;

            dummy_cell = Cell(map_data,x-0.5*dx,x+0.5*dx);

            if( dummy_cell.get_temperature() >= TEMPERATURE_CUTOFF ){
                net_cooling = -dummy_cell.get_cooling_function(); //Cells store cooling values with negative sign.
                heating_value = dummy_cell.get_heating_function();
            }
            else{ net_cooling = 0.0; }
            
            //Update the file
            mapfile<<dummy_cell.eq(0)<<" "<<dummy_cell.get_pressure();
            mapfile<<" "<<dummy_cell.get_temperature()<<" "<<dummy_cell.get_number_density()<<" "<<net_cooling;
            mapfile<<" "<<dummy_cell.get_meanweight()<<" "<<-dummy_cell.get_cooling_time();
            
            mapfile<<" "<<heating_value<<endl;
        }
        //Blank line for gnuplot!
        mapfile<<endl;
    }

    mapfile.close();

}

#endif
