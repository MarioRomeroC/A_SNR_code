#ifndef INTERPOLATIONS_CPP
#define INTERPOLATIONS_CPP

///Auxiliary functions for interpolation
my_float lineal_interpolation(my_float X, my_float x0, my_float y0, my_float x1, my_float y1){

    if(x0 != x1){
        return y0+(X-x0)*(y1-y0)/(x1-x0);
    }else{
        return y0; //For lineal interpolation this works, because X=x0=x1. For extrapolation this is not formally defined
    }

}

my_float semilog_interpolation(my_float X, my_float x0, my_float y0, my_float x1, my_float y1){

    if(x0 != x1){
        my_float Log_y0   = std::log10(y0);
        my_float Log_y1y0 = std::log10(y1/y0);

        my_float Log_Y = Log_y0 + (X-x0)*Log_y1y0/(x1-x0);

        return std::pow(10.,Log_Y);
    }else{
        return y0;
    }
}

my_float log_interpolation(my_float X, my_float x0, my_float y0, my_float x1, my_float y1){

    if(x0 != x1){
        my_float Log_y0   = std::log10(y0);
        my_float Log_Xx1  = std::log10(X/x0);
        my_float Log_y1y0 = std::log10(y1/y0);
        my_float Log_x1x0 = std::log10(x1/x0);

        my_float Log_Y = Log_y0 + Log_Xx1*Log_y1y0/Log_x1x0;

        if(Log_Y != Log_Y){
            std::cout<<y1<<" "<<y0<<std::endl;
            throw std::runtime_error("Nan found in logarithm.");
        }
        return std::pow(10.,Log_Y);
    }else{
        return y0;
    }
}

my_float select_interpolation(my_float X, my_float x0, my_float y0, my_float x1, my_float y1){
//Originally, it selected between linear and log interpolation. Now is kinda rebundant...
    return lineal_interpolation(X,x0,y0,x1,y1);
}

///Some useful functions
my_float relative_error(my_float real_value, my_float approximation){
    //Written to avoid repetition
    my_float result = std::abs((real_value - approximation) / real_value );
    return result;
}

my_float slope_angle(my_float x0, my_float y0, my_float x1, my_float y1){
    //As title suggests, computes the angle with x axis for a straight line between x0 and x1.
    //Gives the result in degrees.

    my_float slope = (y1-y0)/(x1-x0);
    my_float angle = std::atan(slope);

    return angle*180./PI;
}

bool turning_point(my_float x_l, my_float y_l, my_float x_c, my_float y_c, my_float x_r, my_float y_r){
    //So do you want to find a maximum-minimum. That's easy, compute the slopes and see if changes sign
    //If true, then x_c is the closest point you have to the turning point of your function.

    my_float left_slope  = (y_c - y_l)/( x_c - x_l );
    my_float right_slope = (y_r - y_c)/( x_r - x_c );

    if(left_slope*right_slope < 0.){ return true; }
    else{ return false; }

}


bool in_front(Cell* prev_cell, Cell* next_cell){
    //You are in shock or cold front if divergence of velocity is negative
    my_float dx = next_cell->location()-prev_cell->location();
    my_float grad_v = (next_cell->get_velocity()-prev_cell->get_velocity())/dx;
    if(grad_v < 0.0){ return true; }
    else{ return false; }
}

bool in_shock(Cell* prev_cell, Cell* next_cell, Cell* curr_cell){

    const my_float mach_limit = 0.; //Not 1 because I want to track shockwaves degraded to waves but not 0 to protect against numerical oscillations
    bool is_shocked;
    my_float dx = next_cell->location()-prev_cell->location();
    my_float grad_v = (next_cell->get_velocity()-prev_cell->get_velocity())/dx;
    if(grad_v < 0.0){
        //OK, now check density and temperature
        #if ENTROPY_FOR_SHOCKS == 0
            my_float grad_d = (next_cell->eq(0)-prev_cell->eq(0))/dx;
        #else
            my_float next_S = next_cell->eq(3) / next_cell->eq(0);
            my_float prev_S = prev_cell->eq(3) / prev_cell->eq(0);
            //my_float grad_d = (next_cell->eq(3)-prev_cell->eq(3))/dx;
            my_float grad_d = (next_S - prev_S)/dx;
        #endif // ENTROPY_FOR_SHOCKS
        my_float grad_T = (next_cell->get_temperature()-prev_cell->get_temperature())/dx;
        if(grad_d*grad_T > 0.0){
            //Finally, mach numbers
            my_float prev_M = std::abs(prev_cell->get_velocity())/prev_cell->get_magnetosonic();
            my_float next_M = std::abs(next_cell->get_velocity())/next_cell->get_magnetosonic();
            my_float curr_M = std::abs(curr_cell->get_velocity())/curr_cell->get_magnetosonic();
            if( curr_M > mach_limit || prev_M > mach_limit || next_M > mach_limit){
                is_shocked = true;
            }else{
                is_shocked = false;
            }
        }else{
            is_shocked = false;
        }
    }else{
        is_shocked = false;
    }

    return is_shocked;
}

bool in_shock_or_ambient(Cell* prev_cell, Cell* next_cell){
    //Checks if first and second cells are in shock or ambient.

    /*
        Shock conditions:
        -Divergence of velocity is negative
        -The product between gradients of temperature an density is positive.
        -Mach number is higher than 1.

        In practice, these conditions is check the sign of a derivative
    */

    bool is_shocked;
    my_float dx = next_cell->location()-prev_cell->location();
    my_float grad_v = (next_cell->get_velocity()-prev_cell->get_velocity())/dx;
    if(grad_v < 0.0){
        //OK, now check density and temperature
        #if ENTROPY_FOR_SHOCKS == 0
            my_float grad_d = (next_cell->eq(0)-prev_cell->eq(0))/dx;
        #else
            my_float grad_d = (next_cell->eq(3)-prev_cell->eq(3))/dx;
        #endif
        my_float grad_T = (next_cell->get_temperature()-prev_cell->get_temperature())/dx;
        if(grad_d*grad_T > 0.0){
            //Finally, mach numbers
            my_float prev_M = prev_cell->get_velocity()/prev_cell->get_magnetosonic();
            my_float next_M = next_cell->get_velocity()/next_cell->get_magnetosonic();
            if(prev_M > 1. || next_M > 1.){
                is_shocked = true;
            }else{
                is_shocked = false;
            }
        }else{
            is_shocked = false;
        }
    }else{
        is_shocked = false;
    }

    /*
        Ambient conditions:
        -Both cells have same values. (except entropy)
    */
    bool is_ambient = true;
    for(int elem=0;elem<NUM_EQUATIONS-1;++elem){
        if(prev_cell->eq(elem) != next_cell->eq(elem) ){
            is_ambient = false;
            break;
        }
    }

    return is_ambient | is_shocked; //'|' is logic OR.
}

#endif
