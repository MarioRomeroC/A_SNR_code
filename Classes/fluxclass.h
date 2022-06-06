#ifndef FLUXCLASS_H_INCLUDE
#define FLUXCLASS_H_INCLUDE

/******

    Flux class stores the fluxes of each differential (Euler) equations, as well as their source/sink terms.
        The only source/sink that is not included here is a geometric source needed for cylindrical/spherical coordinates

    It's not an inherited class, but need access to Cell members. It's mean to be used behind the curtains

******/

class Flux{

  //friends:
        friend class Cell;

    private:

        my_float f[NUM_EQUATIONS]; //= fluxes

        //Auxiliary variables
        my_float a_max;
        my_float a_min;
        my_float r; //Usually interface location
        my_float P;
        my_float y[NUM_EQUATIONS];

    public:
        ///CONSTRUCTORS
        Flux(){}
        Flux(const Flux& my_copy){
            for(int elem=0;elem<NUM_EQUATIONS;++elem){
                f[elem] = my_copy.f[elem];
                y[elem] = my_copy.y[elem];
            }

            a_max = my_copy.a_max;
            a_min = my_copy.a_min;

            P = my_copy.P;
            r = my_copy.r;
        }
        
        Flux(Cell& my_cell){

            ///Appoint maker
            Cell* maker;
            maker = &my_cell;

            ///Declare useful variables
            my_float v = maker->get_velocity();

            ///Construct auxiliary variables
            a_max = v + maker->get_magnetosonic();
            a_min = v - maker->get_magnetosonic();
            r = maker->location();
            P = maker->P;

            for(int elem=0;elem<NUM_EQUATIONS;++elem){
                y[elem] = maker->eq(elem);
            }

            ///Construct the fluxes
            //-----------------
            //Mass conservation equation
            f[0] = maker->eq(2);

            //-----------------
            //Energy conservation equation
            f[1] = v * (maker->eq(1) + P);

            //-----------------
            //Momentum conservation equation
            f[2] = maker->eq(0) * v * v + maker->P;

            //-----------------
            //Entropy conservation equation
            f[3] = maker->eq(3) * v;


            ///That's all!
        }

        ///ACCESS
        my_float get_flux(int index){
            if(index >= NUM_EQUATIONS || index < 0){
                throw std::out_of_range("index out of range");
            }else{
                return f[index];
            }
        }
        my_float eq(int index){
            if(index >= NUM_EQUATIONS || index < 0){
                throw std::out_of_range("index out of range");
            }else{
                return y[index];
            }
        }

        ///GIVE ME SOME VALUES

        my_float local_speed_max(){return a_max;}
        my_float local_speed_min(){return a_min;}
        my_float location(){return r;}
};

///Some cell public functions that have Flux classes, due to compiler complaining

Flux Cell::right_flux(Cell& prev, Cell& next){
    //Returns a flux at the right interface
    Cell curr_E = this->right_interface(prev,next);

    return Flux(curr_E);
}

Flux Cell::left_flux(Cell& prev, Cell& next){
    //Returns a flux at the left interface
    Cell curr_W = this->left_interface(prev,next);

    return Flux(curr_W);
}

///Class Source

class Source{
    //friends:
        friend class Cell;

    private:

        my_float s[NUM_EQUATIONS];

    public:
        //Constructors
        Source(){}
        Source(Cell& my_cell){
            ///Create and appoint maker
            Cell* maker;
            maker = &my_cell;

            ///Construct the sources
            //-----------------
            //Mass conservation sources
            s[0] = 0.; //Nothing at the moment

            //-----------------
            //Energy conservation sources
            s[1] = 0.0;
            if(maker->get_temperature() > TEMPERATURE_CUTOFF){
                //Add cooling
                #if ACTIVATE_COOLING == 1 //This directive will be moved to updated_thermodynamics when I remove grackle
                s[1]  = maker->get_internal_energy()/maker->get_cooling_time();
                //^^ Has a plus sign because cooling gives a negative cooling function (and heating has positive values).
                #endif
                //Add heating
                s[1] += maker->get_heating_function();
                    
            }else{
                s[1] = 0.0;
            }
                

            //-----------------
            //Momentum conservation sources
            s[2] = 0.; //Remember that geometric terms are included in 'Cell::get_derivative', not here

            //-----------------
            //Entropy conservation sources
            s[3] = s[1]*(maker->gamma - 1.)/(std::pow(maker->eq(0), maker->gamma - 1.)); //Radiative cooling plus heating
            

        }
        ///ACCESS
        my_float get_source(int index){
            if(index >= NUM_EQUATIONS || index < 0){
                throw std::out_of_range("index out of range");
            }else{
                return s[index];
            }
        }
};

#endif
