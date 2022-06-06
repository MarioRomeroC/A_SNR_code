#ifndef INCLUDE_BILINEAR
#define INCLUDE_BILINEAR

/*
 *Same bilinear class as MarioRomeroC/NaiveRT (although there is called 'Table'). In order to work, you need
 *-A 1d array of x of length Nx
 *-A 2d arrays of y and z, having dimensions of Nx x N2.
 *All arrays MUST be ordered, from smallest to highest values.
 *Each x[i] has an array of y, whose values can differ from other x[j], although each y[i] must have same lengths (N2)
 *This way you can work with structured tables in x but not making a perfect square with values in y (as long as each x have the same number of y values attached).
 *I'm using binary seach instead of 'finding intervals', making a more efficient (and hopefully faster) code.
 *
 *As standalone file, you need "<cstdlib>" and "<cmath>" libraries for std::abs() and std::log10() functions.
 *Also add #define my_float double at the beginning
 *
 *Mario Romero           November 2020
 */

/**Search functions
 * I wrote them inside a class in order to encapsulate binary_search
 * Also as static functions. I only need one copy
 **/
class Search{
    private:
        static void binary_search(const my_float fx, const my_float* f, int& i1, int& i2){//That const my_float* may give an error
            //This only works if fx is between *f, if you want the full version (using extrapolation), call find_interval instead.
            
            //Check if we already found the interval
            int diff = i2-i1;
            switch(diff){
                case 1: //interval found
                    break;
                case 0:
                    //Should not happen, but means that i1=i2, so get the closest next point.
                    if( i1 != 0 && std::abs(f[i1-1]-fx) <= std::abs(f[i1+1] - fx) ){
                        i1--;
                    }
                    else{
                        i2++;
                    }
                    break;
                default:
                    //Get a middle index
                    int im = i1 + diff/2;
                    
                    my_float factor = (f[i1]-fx)*(f[im]-fx);
                    if(factor <= 0){ //fx is inside the interval i1 to im
                        i2 = im;
                    }else{ //In this case, between im to i2
                        i1 = im;
                    }
                    binary_search(fx,f,i1,i2);
            }
        }
    public:
        static void interval(const my_float fx, const my_float* f, const int len, int& i1, int& i2, bool& warn){ //That const my_float* may give an error
            //Given an (ordered) array f and a point fx that are read-only (thus the 'const'), find the nearest indexes in order to interpolate (or extrapolate)
            //The return of this function is i1 and i2, as they are called by reference
            
            if(fx <= f[0]){
                i1 = 1;
                i2 = 10;
                warn = true;
            }
            else if(fx >= f[len-1]){
                i1 = len-10;//2
                i2 = len-2;
                warn = true;
            }
            else{
                //Prepare for binary search
                i1 = 0;
                i2 = len-1;
                binary_search(fx,f,i1,i2);
            }
            
        }
};

/**Proper class**/
class Bilinear{
    
    private:
        my_float *x;
        int Nx;
        
        my_float **y;
        my_float **z;
        int N2;
        
        bool Xextrapolation_warning_given = false;
        bool Yextrapolation_warning_given = false;
        
        my_float bilinear_interpolation(my_float X, my_float Y){
            //Find the indexes
            int i1,i2,j1,j2,k1,k2;
            find_indexes(X,Y,i1,i2,j1,j2,k1,k2);
            
            //Get coordinates and points
            my_float x1  = x[i1];
            my_float x2  = x[i2];
            
            my_float y11 = y[i1][j1];
            my_float y12 = y[i1][j2];
            my_float y21 = y[i2][k1];
            my_float y22 = y[i2][k2];
            
            my_float z11 = z[i1][j1];
            my_float z12 = z[i1][j2];
            my_float z21 = z[i2][k1];
            my_float z22 = z[i2][k2];
            
            //Bilinear interpolation is to use z = a*x+b*y+c*x*y+d. We need 4 points to solve the system
            //I write here the solutions
            my_float a = ( (z11-z22)*(y11-y12)*(y21-y22) + y22*(z21-z22)*(y11-y12) - y11*(z11-z12)*(y21-y22) ) / ( (x1-x2)*(y11-y12)*(y21-y22) );
            my_float b = ( (z11-z12)*(y21-y22)*x2 - (z21-z22)*(y11-y12)*x1 ) / ( (x2-x1)*(y11-y12)*(y21-y22) );
            my_float c = ( (z11-z12)*(y21-y22) - (z21-z22)*(y11-y12) ) / ( (x1-x2)*(y11-y12)*(y21-y22) );
            my_float d = z11 - a*x1 - b*y11 - c*x1*y11; //Total laziness
            
            return a*X+b*Y+c*X*Y+d;
        }
        
        my_float bilog_interpolation(my_float X, my_float Y){
            //Find the indexes
            int i1,i2,j1,j2,k1,k2;
            find_indexes(X,Y,i1,i2,j1,j2,k1,k2);
            
            //Get coordinates and points
            my_float x1  = std::log10(x[i1]);
            my_float x2  = std::log10(x[i2]);
            
            my_float y11 = std::log10(y[i1][j1]);
            my_float y12 = std::log10(y[i1][j2]);
            my_float y21 = std::log10(y[i2][k1]);
            my_float y22 = std::log10(y[i2][k2]);
            
            my_float z11 = std::log10(z[i1][j1]);
            my_float z12 = std::log10(z[i1][j2]);
            my_float z21 = std::log10(z[i2][k1]);
            my_float z22 = std::log10(z[i2][k2]);
            
            //Bilinear interpolation is to use z = a*x+b*y+c*x*y+d. We need 4 points to solve the system
            //I write here the solutions
            my_float a = ( (z11-z22)*(y11-y12)*(y21-y22) + y22*(z21-z22)*(y11-y12) - y11*(z11-z12)*(y21-y22) ) / ( (x1-x2)*(y11-y12)*(y21-y22) );
            my_float b = ( (z11-z12)*(y21-y22)*x2 - (z21-z22)*(y11-y12)*x1 ) / ( (x2-x1)*(y11-y12)*(y21-y22) );
            my_float c = ( (z11-z12)*(y21-y22) - (z21-z22)*(y11-y12) ) / ( (x1-x2)*(y11-y12)*(y21-y22) );
            my_float d = z11 - a*x1 - b*y11 - c*x1*y11; //Total laziness
            
            my_float logX = std::log10(X);
            my_float logY = std::log10(Y);
            
            return a*logX+b*logY+c*logX*logY+d;
        }
        
        static void warning_message(my_float point,bool is_X){
            std::cout<<"Warning: I extrapolated in";
            if(is_X){
                std::cout<<" X ";
            }
            else{
                std::cout<<" Y ";
            }
            std::cout<<"values!"<<std::endl;
            std::cout<<"Value used: "<<point<<std::endl;
        }
        
    protected:
        void find_indexes(my_float X, my_float Y, int& i1, int& i2, int& j1, int& j2, int& k1, int& k2){
            /*In order to do a bilinear interpolation, we need four points:
             * (1) (x[i1],y[i1][j1]) -> z[i1][j1]
             * (2) (x[i1],y[i1][j2]) -> z[i1][j2]
             * (3) (x[i2],y[i2][k1]) -> z[i2][k1]
             * (4) (x[i2],y[i2][k2]) -> z[i2][k2]
             * 
             * This function finds the 6 indexes needed for that
             */
            
            //Finding the indexes of x
            //find_interval(X,x,Nx,i1,i2);
            bool X_warning = false;
            Search::interval(X,x,Nx,i1,i2,X_warning);
            if(!Xextrapolation_warning_given && X_warning){
                warning_message(X,true);
                Xextrapolation_warning_given = true;
            }
            
            bool Y_warning = false;
            //Each x has an y[i] attached
            Search::interval(Y,y[i1],N2,j1,j2,Y_warning);
            Search::interval(Y,y[i2],N2,k1,k2,Y_warning);
            if(!Yextrapolation_warning_given && Y_warning){
                warning_message(Y,false);
                Yextrapolation_warning_given = true;
            }
            
        }
    
    public:
        Bilinear(){
            x = nullptr;
            y = nullptr;
            z = nullptr;
        }
        Bilinear(my_float* X, my_float** Y, my_float** Z, int n_x, int n_2){ //Remember, X,Y and Z must be ordered!
            //Init the x
            Nx = n_x;
            x = new my_float[Nx];
            for(int i=0;i<Nx;i++){
                x[i] = X[i];
            }
            
            //Init y and z
            N2 = n_2;
            y = new my_float*[Nx];
            z = new my_float*[Nx];
            for(int i=0;i<Nx;++i){
                y[i] = new my_float[N2];
                z[i] = new my_float[N2];
            }
            for(int i=0;i<Nx;++i){
                for(int j=0;j<N2;++j){
                    y[i][j] = Y[i][j];
                    z[i][j] = Z[i][j];
                }
            }
        }
        //C++ rule of 3: If you have pointer members to dynamic memory, you need a destructor, a copy constructor, and assigment operator.
        ~Bilinear(){
            if( x != nullptr){
                delete[] x;
                x = nullptr;
            }
            if( y != nullptr){
                for(int i=0; i<Nx;i++){
                    delete[] y[i];
                }
                delete[] y;
                y = nullptr;
            }
            if( z != nullptr){
                for(int i=0; i<Nx;i++){
                    delete[] z[i];
                }
                delete[] z;
                z = nullptr;
            }
        }
        Bilinear(const Bilinear& other){
            //This is to avoid shallow copy (copy of pointer members points to the same address as the original) and awful deletes
            
            //Copy x
            x = new my_float[other.Nx];
            Nx = other.Nx;
            for(int i=0;i<Nx;i++){
                x[i] = other.x[i];
            }
            
            //Copy y and z
            N2 = other.N2;
            y = new my_float*[Nx];
            z = new my_float*[Nx];
            for(int i=0;i<Nx;++i){
                y[i] = new my_float[N2];
                z[i] = new my_float[N2];
            }
            for(int i=0;i<Nx;i++){
                for(int j=0;j<N2;j++){
                    y[i][j] = other.y[i][j];
                    z[i][j] = other.z[i][j];
                }
            }
            
        }
        Bilinear& operator=(const Bilinear& other){
            //This represent the operation *this = other. WE ARE OVERRIDING *this, with potential memory leaks if we are not careful!
            
            //Override *this first (but not assign to null)
            if( x != nullptr){delete[] x;}
            if( y != nullptr){
                for(int i=0;i<Nx;++i){ delete[] y[i];}
                delete[] y;
            }
            if( z != nullptr){
                for(int i=0;i<Nx;++i){ delete[] z[i]; }
                delete[] z;
            }
            
            //Repeat the copy constructor code
            //Copy x
            x = new my_float[other.Nx];
            Nx = other.Nx;
            for(int i=0;i<Nx;i++){
                x[i] = other.x[i];
            }
            
            //Copy y and z
            N2 = other.N2;
            y = new my_float*[Nx];
            z = new my_float*[Nx];
            for(int i=0;i<Nx;++i){
                y[i] = new my_float[N2];
                z[i] = new my_float[N2];
            }
            for(int i=0;i<Nx;i++){
                for(int j=0;j<N2;j++){
                    y[i][j] = other.y[i][j];
                    z[i][j] = other.z[i][j];
                }
            }
            
            return *this;
        }
        
        my_float get_value(my_float X, my_float Y){
            return bilinear_interpolation(X,Y);
        }
        
        my_float get_log(my_float X, my_float Y){
            return std::pow(10.0,bilog_interpolation(X,Y));
        }
        
};

#endif
