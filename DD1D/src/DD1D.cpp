#include <iostream>
#include <cmath>

// Output part
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "DD1D.hpp"

namespace Tools
{
    scalar sign(const scalar x)
    {
        return (x > 0) - (x < 0);
    }
}

FVM::cell::cell(const label i):
    cellID(i)
{}

FVM::cell::~cell()
{}

FVM::face::face(const cell& Lcell, const cell& Rcell):
    cellR(Rcell),cellL(Lcell)
{}

FVM::face::~face()
{}

void FVM::face::interpolate()
{
    using namespace Const;

    const scalar delta = 0.5; // coef for upwinding(really?)

    #define interpolateDisc(F)             \ 
    F##L = cellL.F + 0.5*dx*cellL.s##F;    \ 
    F##R = cellR.F - 0.5*dx*cellR.s##F;    \
    F##0 = delta*F##L + (1.0-delta)*F##R;  \
    s##F##0 = (cellR.F-cellL.F)/dx;
    //s##F##0 = delta*cellL.s##F + (1.0-delta)*cellR.s##F; //this seems not working well



    /*
    #define interpolateDisc(F)             \   
    F##L = cellL.F;                        \
    F##R = cellR.F;                        \
    F##0 = delta*F##L + (1.0-delta)*F##R;  \
    s##F##0 = (cellR.F-cellL.F)/dx; */

    interpolateDisc(phi);
    interpolateDisc(velocity);
    interpolateDisc(diff_const);

}

void FVM::face::compute_flux(const scalar dt)
{
    using namespace Const;

    interpolate();

    flux = (velocity0*phi0 - diff_const0*sphi0)*dt;
    //std::cout<< " "<< cellL.phi << ' ' << cellR.phi << ' ' << phi0 << ' ' << sphi0 <<std::endl;
    //std::cout<< cellL.diff_const << ' ' << cellR.diff_const << ' ' << diff_const0 << ' ' << sdiff_const0 <<std::endl;
}

FVM::Solver::Solver():
    step(0),runTime(0),dt(Const::dt)
{
using Const::meshsize;
for(label i=0;i<=meshsize+1;i++)
{
    cells.emplace_back(i); //cell with index 0 and meshsize+1 will be ghost cell
}

for(label i=0;i<=meshsize;i++)
{
    faces.emplace_back(cells[i],cells[i+1]);
}
}

FVM::Solver::~Solver()
{}

void FVM::Solver::initialization()
{
    using namespace Const;
    for(label i=1;i<=meshsize;i++)
    {
        cells[i].phi = initial_condition((i-0.5)*dx);
        cells[i].velocity = compute_velocity((i-0.5)*dx);
        cells[i].svelocity = 0.0;
        cells[i].diff_const = D;
        cells[i].sdiff_const = 0.0;

        std::cout << "init: " << cells[i].phi << std::endl;
    }

    // set values of ghost cells
    cells[0].phi = cells[1].phi;
    cells[meshsize+1].phi = cells[meshsize].phi;

}

void FVM::Solver::reconstruct()
{
    using namespace Const;
    using namespace Tools;

    cells[1].sphi = (cells[2].phi-cells[1].phi)/dx;
    cells[meshsize].sphi = (cells[meshsize].phi-cells[meshsize-1].phi)/dx;
    cells[1].svelocity = (cells[2].velocity-cells[1].velocity)/dx;
    cells[meshsize].svelocity = (cells[meshsize].velocity-cells[meshsize-1].velocity)/dx;
    cells[1].sdiff_const = (cells[2].diff_const-cells[1].diff_const)/dx;
    cells[meshsize].sdiff_const = (cells[meshsize].diff_const-cells[meshsize-1].diff_const)/dx;

    for(label i=2;i<=meshsize-1;i++)
    {
        if(First_Order_Reconstruction) break;
        cell *cellL = &cells[i-1];
        cell *cellC = &cells[i];
        cell *cellR = &cells[i+1];

        scalar sL,sR;
        #define reconstructDisc(F) \
        sL = (cellC->F - cellL->F)/dx;  \
        sR = (cellR->F - cellC->F)/dx;  \
        cellC->s##F = (sL+sR)/2.0; 
        //cellC->s##F = ( sign(sL) + sign(sR) ) * abs(sR) * abs(sL)/ ( abs(sR) + abs(sL) + 1e-50 ); // minmod limiter, this makes figure fluffy

        reconstructDisc(phi);
        reconstructDisc(velocity);
        reconstructDisc(diff_const);
    }
    cells[0].sphi = -cells[1].sphi;
    cells[meshsize+1].sphi = -cells[meshsize].sphi;
    cells[0].svelocity = -cells[1].svelocity;
    cells[meshsize+1].svelocity = -cells[meshsize].svelocity;
    cells[0].sdiff_const = -cells[1].sdiff_const;
    cells[meshsize+1].sdiff_const = -cells[meshsize].sdiff_const;

}

void FVM::Solver::evolve()
{
    for(auto &face : faces)
    {
        face.compute_flux(dt);
    }
}

void FVM::Solver::update()
{
    using namespace Const;
    for(label i=1;i<=meshsize;i++)
    {
        cells[i].phi = cells[i].phi - (faces[i].flux - faces[i-1].flux)/dx;
    }
    // set values of ghost cells
    cells[0].phi = cells[1].phi + 2;
    cells[meshsize+1].phi = cells[meshsize].phi-2;

    runTime += dt;
    step += 1;
}

void FVM::Solver::write(std::string str)
{
    std::stringstream fileName;

    // File name like "prim1000.dat"
    fileName << "./result/phi" << str << ".dat" << std::endl;

    std::string s;

    fileName >> s;

    std::ofstream out(s.c_str());

    // Write as Tecplot format
    out << "X Density"
        << std::endl;

    for (label i = 1; i <= Const::meshsize + 1; i++ )
    {
        const scalar cellCenter = (i-0.5)*Const::dx;

        scalar density = cells[i].phi;

        out << std::setprecision(10)
            << cellCenter           << " "
            << density                  
            << std::endl;
    }
    out.close();
}

void FVM::Solver::info()
{
    using std::cout;
    using std::endl;

    cout << "step = "    << step     << "\t" 
         << "runTime = " << runTime  << "\t" 
         << "dt = "      << dt      
         << endl;
}

void FVM::Solver::iterate()
{

    initialization();
    write("0");

    while (runTime < Const::stopTime)
    {
        reconstruct();
        evolve();
        update();

        if (step % Const::writeInterval == 0)
        {
            write( std::to_string(step) );
            info();
        }

        
    }
    
    write("End");
}

int main(int argc, char* argv[])
{
using namespace Const;
FVM::Solver DD1D;
DD1D.iterate();
return 0;
}