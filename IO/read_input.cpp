#ifndef READ_INPUT_CPP
#define READ_INPUT_CPP

int count_lines(string name){
    //Given the name of a file, open it and count the number of lines

    int number_of_lines = 0;
    ifstream file;
    file.open(name);

    if(!file.is_open()){
        throw runtime_error("File not found");
    }else{
        //Count lines otherwise
        string line;
        while(!file.eof()){
            getline(file,line);
            if(line!=""){ number_of_lines++; } //Don't count blank lines
        }
    }

    file.close();
    return number_of_lines;
}

my_float* append(my_float* old_array,my_float new_value, int& size){ //(!) Why is this here?
    
    my_float* new_array = new my_float [size+1];
    
    //Copy array.
    for (int i=0;i<size;++i){
        new_array[i] = old_array[i];
    }
    //Add the new value
    new_array[size] = new_value;
    
    //Delete old array
    delete[] old_array;
    //Change size
    size++;
    
    return new_array;
    
}

void smooth_shock(Cell* grid, const int N_cells){ //(!) When is this called?
    //If this is used before the first refinement, cells in grid are ordered
    /* Ideally, a shock front should be two cells: an upstream and downstream values.
       Numerically, a shock front is composed by three/four cells (at least from the results I made from this code)
       If you start with a shock made of two cells, the computer interprets a different shock compared to a three/four-cells' shock.
       This will give a numerical oscillation that can be easily ignored as long as
       the ambient temperature (most evident quantity to change) is several orders of magnitude below shocked temperature.
       This 'overtemperature' is a factor 2 (or less) the expected temperature, and can be a huge deal if we have a shocked and ambient temperatures of 1e7 and 5e6 K, respectively.
    */

    const int shock_cells = 3;

    for(int i=0;i<N_cells-shock_cells - 1;++i){
        //Check if there is a shock between two cells
        my_float v_der = (grid[i+1].get_velocity() - grid[i].get_velocity())/(grid[i+1].location() - grid[i].location());
        my_float T_der = (grid[i+1].get_temperature() - grid[i].get_temperature())/(grid[i+1].location() - grid[i].location());
        my_float d_der = (grid[i+1].eq(0) - grid[i].eq(0))/(grid[i+1].location() - grid[i].location());
        if(v_der < 0. && T_der*d_der > 0.){
            //Create new data
            int i_0 = i;
            int i_f = i+shock_cells;
            for(int j=i;j<i_f;++j){
                Cell new_cell;
                my_float new_data[NUM_EQUATIONS];
                for(int elem=0;elem<NUM_EQUATIONS;++elem){
                    //Select interpolation is defined in AMR_main
                    new_data[elem] = grid[j].eq(elem);
                }
                my_float new_T = select_interpolation(grid[j].location(), grid[i_0].location(), grid[i_0].get_temperature(), grid[i_f].location(), grid[i_f].get_temperature());
                my_float new_gamma = select_interpolation(grid[j].location(), grid[i_0].location(), grid[i_0].get_adiabatic_const(), grid[i_f].location(), grid[i_f].get_adiabatic_const());
                my_float new_Z = select_interpolation(grid[j].location(), grid[i_0].location(), grid[i_0].get_metallicity(), grid[i_f].location(), grid[i_f].get_metallicity());

                my_float new_P = pressure(new_data[0],new_data[1],new_data[2]/new_data[0],new_gamma);
                if( grid[j].get_system() == false || new_P == FAILURE){
                    //Dual energy formalism
                    new_P = pressure(new_data[0],new_data[3],new_gamma);
                    new_cell = Cell(new_data,new_P,new_T,new_gamma,new_Z,grid[j].location()-0.5*grid[j].width(),grid[j].location()+0.5*grid[j].width(),false);
                }else{
                    new_cell = Cell(new_data,new_P,new_T,new_gamma,new_Z,grid[j].location()-0.5*grid[j].width(),grid[j].location()+0.5*grid[j].width(),true);
                }

                my_float new_tcool = -select_interpolation(grid[j].location(), grid[i_0].location(), -grid[i_0].get_cooling_time(), grid[i_f].location(), -grid[i_f].get_cooling_time());
                my_float new_meanweight = select_interpolation(grid[j].location(), grid[i_0].location(), grid[i_0].get_meanweight(), grid[i_f].location(), grid[i_f].get_meanweight());

                new_cell.update_tcool(new_tcool);
                new_cell.update_meanweight(new_meanweight);
                grid[j].soft_update(new_cell); //Preserves the pointers

            }
            i = i_f;
        }
    }
}

void init_from_file(Cell* grid, string name, const int N_cells){
    /*
       Generate the initial conditions from a file given as argument. The order of columns must be

       x   rho   e   rho*v  S [ignored]

       Separated by tabs and without any initial text
       S = modified entropy = pressure / rho^(gamma-1)
       gamma = ADIABAT_CONST in define.h
       Units are in cgs.
    */

    ///Count number of lines of input file
    ifstream infile;
    const int N_lines = count_lines(name); //Throws an exception if you have given a bad filename


    ///Read input file
    /**CHANGE THIS IF ADDICTIONAL VARIABLES ARE ADDED**/
    enum Readables { Position, Mass_density, Energy_density, Momentum_density, Modified_entropy, N_variables};
    my_float** table = new my_float*[N_variables];
    for(int column=0;column<N_variables;++column){
        table[column] = new my_float[N_lines];
    }

    infile.open(name);
    string trash;
    for(int line=0;line<N_lines;++line){
        infile>>table[Position][line];
        infile>>table[Mass_density][line];
        infile>>table[Energy_density][line];
        infile>>table[Momentum_density][line];
        infile>>table[Modified_entropy][line];
        //And columns left are trash
        getline(infile,trash);
    }

    infile.close(); //No longer needed

    ///Generate grid of cells
    my_float x0 = table[Position][0];
    my_float xf = table[Position][N_lines-1];
    my_float dx = (xf - x0)/(N_cells-1.);

    //Let's define interpolating functions
    Slopes rho_x = Slopes(table[Position],table[Mass_density],N_lines,SYMMETRY);
    Slopes ene_x = Slopes(table[Position],table[Energy_density],N_lines,SYMMETRY);
    Slopes mom_x = Slopes(table[Position],table[Momentum_density],N_lines,SYMMETRY);
    Slopes ent_x = Slopes(table[Position],table[Modified_entropy],N_lines,SYMMETRY);

    for(int i = 0; i < N_cells; ++i){
        my_float data[NUM_EQUATIONS]; //Mass, energy, momentum densities; and modified entropy
        my_float xi = x0 + i*dx;

        //Fill data array
        //Ambient values are guaranteed to be constants
        //And non-ambient values are interpolated to be energy-conserved.
        data[0] = rho_x.get_value(xi);
        data[1] = ene_x.get_value(xi);
        data[2] = mom_x.get_value(xi);
        data[3] = ent_x.get_value(xi);

        //Create cell
        grid[i] = Cell(data,xi-0.5*dx,xi+0.5*dx);
        if(i>0){ grid[i].Link_prev(grid[i-1]); }
    }

    ///Deleting dynamic arrays
    for(int column=0;column<N_variables;++column){
        delete[] table[column];
    }
    delete[] table;

}

Bilinear load_table(string name){
    /*
    This function does what its name says: it loads a table.
    In the code, the function is called to create the Equation of state, cooling and heating tables.
    These tables are z(x,y), structured as
    x0 y00 z00
    x0 y01 z01
    x0 y02 z02
    x1 y10 z10
    x1 y11 z11
    etc
    Thus you need two sizes for y (size_y and size_2)
     */
    
    //Count number of lines of input file
    ifstream infile;
    const int N_lines = count_lines(name); //Throws an exception if you have given a bad filename
    
    //Define your arrays (I will need to append, so this will be difficult...)
    int size_x = 1; //Will change over time
    int size_y = 1;
    int size_z = 1;
    int size_2 = 1; //Same size as size_z. This is to avoid problems with the parameter by reference in append function
    my_float *x = new my_float[size_x];
    my_float *y_all = new my_float[size_2];
    my_float *z_all = new my_float[size_z];
    
    infile.open(name);
    
    //Read file
    my_float x_temp, y_temp, z_temp;
    //First line
    infile>>x_temp;
    infile>>y_temp;
    infile>>z_temp;
    x[0] = x_temp;
    y_all[0] = y_temp;
    z_all[0] = z_temp;
    bool y_filled = false; //This is used to know size_y
    my_float epsilon;
    
    for(int line=1;line<N_lines;++line){
        infile>>x_temp;
        infile>>y_temp;
        infile>>z_temp;
        
        //Fill x
        epsilon = abs((x_temp - x[size_x-1])/x[size_x-1]);
        
        if(epsilon > 1e-4){
            x = append(x,x_temp,size_x);
            y_filled = true;
        }
        
        //Fill y
        if( y_filled == false ){
            size_y++;
        }
        y_all = append(y_all,y_temp,size_2);
        
        z_all = append(z_all,z_temp,size_z);
        
    }
    
    infile.close();
    
    //Now fill the z, because we have x,y and their sizes. Not the most elegant, but...
    my_float** z = new my_float* [size_x];
    my_float** y = new my_float* [size_x];
    for(int i=0;i<size_x;++i){
        z[i] = new my_float [size_y];
        y[i] = new my_float [size_y];
    }
    
    int all_count = 0;
    for(int i=0;i<size_x;++i){
        for(int j=0;j<size_y;++j){
            z[i][j] = z_all[all_count];
            y[i][j] = y_all[all_count];
            all_count++;
        }
    }
    
    delete[] y_all; //Not needed
    delete[] z_all;
    
    Bilinear result = Bilinear(x,y,z,size_x,size_y);
    
    delete[] x;
    for(int i=0;i<size_x;++i){
        delete[] y[i];
        delete[] z[i];
    }
    delete[] y;
    delete[] z;
    
    return result;
    
}

#endif
