#ifndef INCLUDE_INTERNOMIALS
#define INCLUDE_INTERNOMIALS

//These classes are for 1D interpolation. Think as a f(x).
//Not only gives you an interpolation of points, but also integrated values (in case you want the special condition of integrals must be conserved)
//These classes are slower in perfomance but more powerful for Eulerian fluid codes or AMR routines in which you need cell-conservation values
//They can work at their own, just be sure to add '#define my_float double', '#define <cmath>' and '#define <algorithm>' at the beginning

class Slopes{
    /*
        Given a set of points x, and their function values y(x). It computes the interpolations

        But also computes a set of integrated points of the function:

        f(x) = f[i] = int_{x[0]}^{x[i]} y(x) dV   //Computed by trapezoidal rule. This code does not assume that you have equispaced points

        Is 1D symmetry dependent (i.e: Cartesian, cylindrical or spherical coordinates), by default I assume cartesian coordinates.
        Therefore I suggest defining a macro variable called SYMMETRY in which 0 is cartesian, 1 cylindrical, and 2 spherical, and enter as argument in 'integrated_value(...)' function.
        Furthermore, this class only works with cmath library activated.

        The values this class gives is an integral conserved quantity of y(x), between two values x_L and x_R (x_R > x_L):
        result = 1/(V(x_R)-V(x_L)) * (f(x_R) - f(x_L))

        Note that if x_R and/or x_L are not points in *x array, a linear interpolation between the values in which x_R > x[i] and x_R < x[i-1]  (same for x_L).
    */

    private:

        my_float *x; //Points
        my_float *y; //Values
        my_float *f; //Integrated values
        int N;     //Their number of points
        int geom;  //Geometry

        my_float linear_interpolation(my_float X, my_float x1, my_float x2, my_float y1, my_float y2){
            //Self-explainatory!
            if(x2 == x1){
                if(X==x1){
                    //If you use X for an interpolation, then it should return y1 (X=x1=x2). Unless you have something that is not a function, y1=y2 should be true
                    return y1;
                }else{
                    //But it you're extrapolating, you have y(x) = 0. You can prove it computing 1/y(x) in order to remove any (x2-x1) denominator
                    std::cout<<"Warning: Bad extrapolation defined!"<<std::endl;
                    return NAN; //However, computer will give you an indetermination. I'm thinking of giving you a nan or 0.
                }
            }else{
                return y1 + (X - x1)*( (y2-y1)/(x2-x1) ); //Vile copy from wikipedia
            }
        }

        my_float linear_conservation(my_float r_R, my_float r_L, my_float fR, my_float fL){
            //Gives the integral described above
            return (geom+1)*(fR - fL) / ( pow(r_R,geom+1) - pow(r_L, geom+1) );
        }

    protected:
        my_float trapezoidal_rule(my_float x_0, my_float x_f, my_float f_0, my_float f_f){
            //Compute the trapezoidal rule of one interval [x_0,x_f] in order to compute the integral
            return (x_f - x_0)*(f_f + f_0)/2.;
        }
        void sort_values(my_float* X, my_float* Y){ //Both arrays' size is N, and their result is *y created
            //Sort *x is easy:
            // std::sort(x, x+N);
            // However we have to rearrange the positions of Y according to the sorted x array
            // so... I don't know if it's better to call sort x, and after rearrange y
            // or sorting with my code (slower), but at the same time assign the y.
            // I'm taking first option

            y = new my_float[N];

            bool *assigned = new bool[N];
            for(int j=0;j<N;++j){
                assigned[j] = false;
            }
            /* Imagine that I give you these arrays:
                x = [ 1, 3, 5, 6, 6, 6, 6, 6]
                y = [ 2, 4, 0, 1, 1, 3, 5, 7]
               You are interpolating something that is not a function.
               Having x[i] == X[j] does not guarantee saving latter points (If you indeed have a function, it doesn't matter, all y(x) is uniquely defined).
            */

            for(int i=0;i<N;++i){
                for(int j=0;j<N;++j){
                    if(x[i] == X[j] && assigned[j] == false){ //Shoudn't give an error of 1 != 1 if is called afterwards initialization and sort, no mathematical operation shoudn't have been done!
                        y[i] = Y[j];
                        assigned[j] = true;
                        break; //Right position found! No need to continue
                    }
                }
            }

            delete[] assigned;

            //Finished! *y has been initialized!
        }

        void find_interval(my_float X, int& i1, int& i2){
            //Find the subinterval [x[i1], x[i2]] in which X is located

            if(X < x[0]){//First case is X < x[0]
                i1 = 0;
                i2 = 1;
            }
            else if(X > x[N-1]){ //Second case is the opposite: X>x[N-1]
                i1 = N-2;
                i2 = N-1;
            }
            else{//Third, and most common case
                for(int i=1;i<N;++i){
                    if( x[i-1]<=X && X<=x[i] ){ //Got them!
                        i1 = i-1;
                        i2 = i;
                        break;
                    }
                }
            }
        }

        void find_nearest(my_float X, int& i1, int& i2){
            //Given x find nearest values (by their indexes i1, i2) in x*
            //i1 and i2 are your outputs, as long they are declared variables in the scope in which this function is called this works fine
            /*
                This generates a side effect if values are not equispaced:
                for example if you want the linear interpolation at x=3.
                and your array of points is x[] = {1,2, 6}. This routine chooses x=1 and x=2 for the linear interpolation, not x=2 and x=6

                This was coded before the correct implementation above, but I think that's a nice feature to keep.
            */

            my_float dist =  1e6 * std::abs(x[N-1]-x[0]); //Arbitrary large distance
            my_float low1 = dist;
            my_float low2 = dist;
            i1 = 0;
            i2 = 1;
            //First nearest
            for(int i=0;i<N;++i){
                dist = std::abs(X-x[i]);
                if( dist < low1 ){
                    low1 = dist;
                    i1   = i;
                }
            }

            //Second nearest;
            for(int i=0;i<N;++i){
                if(i == i1){continue;} //Exclude first nearest for this loop
                dist = std::abs(X-x[i]);
                if( dist < low2 ){
                    low2 = dist;
                    i2   = i;
                }
            }

            //return i1 and i2
        }

    public:
        Slopes(){
            x = nullptr;
            y = nullptr;
            f = nullptr;
        }
        Slopes(my_float* X, my_float* Y, const int n, const int geometry = 0){
            //Init the x
            N = n;
            x = new my_float[N];
            for(int i=0;i<N;i++){
                x[i] = X[i];
            }
            std::sort(x, x+N); //And sort it
            
            //Init the y
            sort_values(X,Y); //*y is 'new'ed there

            //Init the f
            geom = geometry;
            my_float suma = 0.0;
            f = new my_float[N];
            f[0] = suma;
            my_float r0,r1;
            for(int i=1;i<N;++i){
                //In 1d: int{ f dV } = int { f * r^geom dr }, ignoring constants because they will be removed when values are computed
                r0 = std::pow(x[i-1],geom);
                r1 = std::pow(x[i],geom);
                suma += trapezoidal_rule(x[i-1],x[i],y[i-1]*r0,y[i]*r1);
                f[i] = suma;
            }

            //All done!
        }
        ~Slopes(){
            if(x != nullptr){
                delete[] x;
                x = nullptr;
            }
            if(y != nullptr){
                delete[] y;
                y = nullptr;
            }
            if(f != nullptr){
                delete[] f;
                f = nullptr;
            }
        }
        Slopes(const Slopes& other){ //Copy constructor
            //C++ rule of 3: If you have pointer members to dynamic memory, you need a destructor, a copy constructor, and assigment operator.
            //This is to avoid shallow copy (copy of pointer members points to the same address as the original) and awful deletes

            //To generate a deep copy, we have to assign new memory locations to the copy!
            x = new my_float[other.N];
            y = new my_float[other.N];
            f = new my_float[other.N];
            N = other.N;
            geom = other.geom;

            for(int i=0;i<N;i++){
                x[i] = other.x[i];
                y[i] = other.y[i];
                f[i] = other.f[i];
            }

        }
        Slopes& operator=(const Slopes& other){
            //This represent the operation *this = other. WE ARE OVERRIDING *this, with potential memory leaks if we are not careful!

            //Prepare new arrays
            my_float *new_x = new my_float[other.N];
            my_float *new_y = new my_float[other.N];
            my_float *new_f = new my_float[other.N];

            for(int i=0;i<other.N;++i){
                new_x[i] = other.x[i];
                new_y[i] = other.y[i];
                new_f[i] = other.f[i];
            }

            //Now override *this
            if( x != nullptr ){ delete[] x; }
            if( y != nullptr ){ delete[] y; }
            if( f != nullptr ){ delete[] f; }
            //Unlike destructors, {x,y,f} won't be null pointers, but...
            x = new_x;
            y = new_y;
            f = new_f;

            N = other.N;
            geom = other.geom;

            return *this;
        }


        my_float get_value(my_float X, bool use_nearest = false){
            //Find the two nearest points
            int i1;
            int i2;

            if(use_nearest){
                find_nearest(X,i1,i2);
            }else{
                find_interval(X,i1,i2);
            }

            //Get the interpolation
            return linear_interpolation(X, x[i1],x[i2],y[i1],y[i2]);

        }

        my_float get_integrated(my_float X_R, my_float X_L, bool use_nearest = false){
            //We have to do the procedure of get value twice, but our objetive here is get f_R and f_L

            //Right value
            int i1,i2;
            if(use_nearest){
                find_nearest(X_R,i1,i2);
            }else{
                find_interval(X_R,i1,i2);
            }
            my_float f_R = linear_interpolation(X_R, x[i1],x[i2],f[i1],f[i2]);

            //Left value
            int j1,j2;
            if(use_nearest){
                find_nearest(X_L,j1,j2);
            }else{
                find_interval(X_L,j1,j2);
            }

            my_float f_L = linear_interpolation(X_L, x[j1],x[j2],f[j1],f[j2]);

            //Now give the result
            return linear_conservation(X_R,X_L,f_R,f_L);
        }

        ///Get derivatives
        //due to 'geom' being a variable, we need to make a gradient and a divergence

        int stores(my_float X){
            for(int i=0;i<N;++i){
                if(abs(X-x[i]) < EPSILON_FLOAT){ return i; }
            }
            return -1;
        }

        my_float gradient(my_float X, bool use_nearest = false){
            //Regardless of coordinate system, it's always a standard derivative
            //Furthermore, because this is 'connect points by straight line', you return the slope

            //Unlike getting values, for derivatives is important to know if X is stored
            int i = this->stores(X);//-1 is the return value if X is not stored

            
            int i1;
            int i2;

            if(use_nearest){
                find_nearest(X,i1,i2);
            }else{
                find_interval(X,i1,i2);
            }


            //Return the slope between points
            return (y[i2]-y[i1])/(x[i2]-x[i1]);

        }

        my_float divergence(my_float X, bool use_nearest = false){
            //This is x^{-geom}*derivative(x^{geom}*y)
            //Similar but not equal!

            //Unlike getting values, for derivatives is important to know if X is stored
            int i = this->stores(X);//-1 is the return value if X is not stored

            my_float super_x1,super_x2,slope;

            
            int i1;
            int i2;

            if(use_nearest){
                find_nearest(X,i1,i2);
            }else{
                find_interval(X,i1,i2);
            }
            //Compute the 'slope'
            super_x1 = std::pow(x[i1],geom);
            super_x2 = std::pow(x[i2],geom);

            slope = (super_x2*y[i2] - super_x1*y[i1])/(x[i2]-x[i1]);

            my_float super_X = std::pow(X,geom);
            return slope/super_X;
        }

        ///Get integral
        //Literally! not a conserved-values interpolation
        my_float integrate(my_float X_i, my_float X_f, my_float dz = 1.){
            //This is trickier than it looks.
            //We are doing this integral -> int_V y(r) dV
            //There are differences between coordinates, and an ambiguity if geom != 2. This is resolved with dz

            //Find the right points
            int i1; //for X_f
            int i2;

            int j1; //for X_i
            int j2;

            //Finding_nearest makes no sense here!
            find_interval(X_f,i1,i2);
            find_interval(X_i,j1,j2);


            //Make an integral from x[0] to X_f
            //Remember that f[i] = int_{x[0]}^{x[i]} y(x) x^{geom} dx

            my_float r0 = std::pow(x[i1],geom);
            my_float r1 = std::pow(X_f,geom);
            my_float Y  = get_value(X_f);

            my_float Int_f = f[i1] + trapezoidal_rule(x[i1],X_f,y[i1]*r0,Y*r1);
            //^We have to account if X_f is in the middle of the interval, isn't it?

            //Make an integral from x[0] to X_i
            r0 = std::pow(x[j1],geom);
            r1 = std::pow(X_i,geom);
            Y  = get_value(X_i);

            my_float Int_i = f[j1] + trapezoidal_rule(x[j1],X_i,y[j1]*r0,Y*r1);

            //And give the result (remember the correction from geom!)
            my_float result = (Int_f - Int_i);

            switch(geom){

                case 2:
                    return 4.*std::acos(-1)*result;
                case 1:
                    return 2.*std::acos(-1)*dz*result;
                case 0:
                    return dz*dz*result;
            }
        }
};

#endif
