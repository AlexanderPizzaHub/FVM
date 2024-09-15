#include <valarray>
#include <vector>
#include <string>
#include <cmath>

// Typedefs

using label   = long long int;
using scalar  = double;
using tensor1 = std::valarray<scalar>;
using tensor2 = std::valarray<std::valarray<scalar>>;
using std::vector;


namespace Const
{

#define PI 3.14159265
inline namespace Case1
{
    // specify domain and mesh
    inline constexpr label meshsize = 100;
    inline constexpr scalar dx = 0.01;

    // mobility and diffusivity
    inline const scalar compute_velocity(scalar x)
    {
        scalar vel;
        //vel= 1.0*sin(PI*x);
        vel = 1.0;

        return vel;
    };

    // inline constexpr scalar velocity = 0.0;

    inline constexpr scalar D = 0.01;

    // initial condition
    inline const scalar initial_condition(scalar x)
    {
        return (x>0.3 & x<0.7) ? 1.0 : 0.2 ;
        //return x*(1-x);
        //return 1+x;
    };

    inline constexpr scalar dt=0.00001;
    inline constexpr scalar stopTime = 1.0;

    inline constexpr label writeInterval = 1000;
    inline constexpr bool First_Order_Reconstruction=false;


}

}

namespace FVM
{

class cell
{
    public:
    
    explicit cell(const label i);
    ~cell();

    label cellID; 
    scalar centrod;
    scalar volume;

    scalar phi;  // flow variable 
    scalar sphi; // slope of flow variable


    scalar velocity;
    scalar svelocity; 
    scalar diff_const; 
    scalar sdiff_const;

};

class face
{
    public:

    explicit face(const cell& cellL,const cell& cellR);
    ~face();

    const cell& cellL;
    const cell& cellR;

    scalar phiL,sphiL,velocityL,svelocityL,diff_constL,sdiff_constL; // interpolated variables
    scalar phiR,sphiR,velocityR,svelocityR,diff_constR,sdiff_constR; // interpolated variables
    scalar phi0,sphi0,velocity0,svelocity0,diff_const0,sdiff_const0; // interpolated variables

    scalar flux;

    void interpolate(); // interpoalte variables to face
    void compute_flux(const scalar dt); // c\ompute flux at face

};

class Solver
{

    public:
    vector<cell> cells; 
    vector<face> faces;

    scalar dt;

    label step;
    scalar runTime;

    void initialization();

    void getDt();

    void reconstruct();

    void evolve();

    void update();

    void write(std::string str);

    void iterate();

    void info();

    
    Solver();
    ~Solver();
};

}