#ifndef WRITE_FILE_CPP
#define WRITE_FILE_CPP

//Filename is self-explanatory

void write_output(Cell* grid, Cell initial_ambient, const my_float my_time = 0.0,const int N = NCELLS){

    ///Creating auxiliary pointers
    Cell* curr;
    Cell* last;
    curr = grid;
    last = grid;
    last = find_last(last);

    ///Creating the path
    string path_to_file = OUTPUT_PATH;
    string first_name = "d"+tostring(initial_ambient.eq(0)/HYDROGEN_MASS)+"_Z"+tostring(INIT_METALLICITY)+"_P"+tostring(initial_ambient.get_pressure()/BOLTZMANN_CTE);//OUTPUT_FILENAME;
    string second_name = "_t=";
    string third_name = tostring(my_time);
    string fourth_name = "kyr.txt";

    #if DEBUG == 0
        string name = path_to_file + first_name + second_name + third_name + fourth_name;
    #else
        string name = "DEBUG.txt";
    #endif

    std::cout<<"Creating file "<<name<<std::endl;

    ///Creating the file
    ofstream file;
    file.open(name.c_str());

    ///Creating the intro
    #if KEY_ON == 1
        //Location and hydrodynamic variables first
        file<<" location | "<<"mass den | "<<"energy den | "<<"momentum den | "<<"mod entropy || ";
        //thermodynamic variables next
        file<<"pressure | "<<"temperature | "<<" number den | "<<"internal energy || ";
        //Other variables then
        file<<" velocity | "<<"cooling function | "<<" metallicity | "<<" gamma | "<<" meanweight ||";
        #if AMR_ON == 1
            //Adaptive Mesh Refinement variables last
            file<<" LoR | "<<" dx";
        #endif
        file<<std::endl;
    #else
        //Generate a separate file with the key
        ofstream keyfile;
        keyfile.open(path_to_file+"key.txt");

        keyfile<<"1  : location"<<std::endl;
        keyfile<<"2  : mass density"<<std::endl;
        keyfile<<"3  : energy density"<<std::endl;
        keyfile<<"4  : momentum density"<<std::endl;
        keyfile<<"5  : modified entropy"<<std::endl;
        keyfile<<"6  : pressure"<<std::endl;
        keyfile<<"7  : temperature"<<std::endl;
        keyfile<<"8  : number density"<<std::endl;
        keyfile<<"9  : internal energy"<<std::endl;
        keyfile<<"10 : velocity"<<std::endl;
        keyfile<<"11 : cooling function"<<std::endl;
        keyfile<<"12 : metallicity"<<std::endl;
        keyfile<<"13 : gamma"<<std::endl;
        keyfile<<"14 : mean weight"<<std::endl;
        #if AMR_ON == 1
            keyfile<<"15 : Level of Refinement"<<std::endl;
            keyfile<<"16 : Cell width"<<std::endl;
            keyfile<<"17 : Shocked cell"<<std::endl;
            keyfile<<"18 : Mach number"<<std::endl;
            keyfile<<"19 : heating function (erg/cm3/s)"<<std::endl;
        #else
            keyfile<<"15 : Mach number"<<std::endl;
            keyfile<<"16 : heating function (erg/cm3/s)"<<std::endl;
        #endif

        keyfile<<"--------------------"<<std::endl;
        keyfile<<"All units are in cgs, outputs are given with "<<NUMBER_PRECISION<<" digits"<<std::endl;
        keyfile<<"And all calculations are done with "<<(FLOAT_PRECISION+1)*8<<" digits"<<std::endl;

        keyfile.close();
    #endif

    ///Writing the file
    file<<std::setprecision(NUMBER_PRECISION);
    while(curr->next != nullptr){

        file<<curr->location()<<" "<<curr->eq(0)<<" "<<curr->eq(1)<<" "<<curr->eq(2)<<" "<<curr->eq(3)<<" ";
        file<<curr->get_pressure()<<" "<<curr->get_temperature()<<" "<<curr->get_number_density()<<" "<<curr->get_internal_energy()<<" ";
        file<<curr->get_velocity()<<" "<<-curr->get_cooling_function()<<" "<<curr->get_metallicity()<<" "<<curr->get_adiabatic_const()<<" "<<curr->get_meanweight()<<" ";
        #if AMR_ON == 1
            file<<curr->level_of_refinement()<<" "<<curr->width()<<" "<<curr->shocked_cell();
        #endif
        file<<" "<<curr->get_velocity()/curr->get_magnetosonic()<<" "<<curr->get_heating_function();

        file<<std::endl;

        curr = curr->next;
    }
    //Due to while condition, last cell is not written, we can solve this adding these lines
    file<<curr->location()<<" "<<curr->eq(0)<<" "<<curr->eq(1)<<" "<<curr->eq(2)<<" "<<curr->eq(3)<<" ";
    file<<curr->get_pressure()<<" "<<curr->get_temperature()<<" "<<curr->get_number_density()<<" "<<curr->get_internal_energy()<<" ";
    file<<curr->get_velocity()<<" "<<-curr->get_cooling_function()<<" "<<curr->get_metallicity()<<" "<<curr->get_adiabatic_const()<<" "<<curr->get_meanweight()<<" ";
    #if AMR_ON == 1
        file<<curr->level_of_refinement()<<" "<<curr->width()<<" "<<curr->shocked_cell();
    #endif
    file<<" "<<curr->get_velocity()/curr->get_magnetosonic()<<" "<<curr->get_heating_function();

    file<<std::endl;

    ///Closing the file
    file.close();

    ///That's all
}

#endif
