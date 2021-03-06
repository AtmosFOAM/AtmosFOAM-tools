#include "sinusoidalTracerField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(sinusoidalTracerField, 0);
addToRunTimeSelectionTable(tracerField, sinusoidalTracerField, dict);

sinusoidalTracerField::sinusoidalTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
    tracerField(velocityField),
    
    mag(dict.lookupOrDefault<scalar>("mag", scalar(0))),
    
    a_x(dict.lookupOrDefault<scalar>("a_x", scalar(0))),
    a_y(dict.lookupOrDefault<scalar>("a_y", scalar(0))),
    a_z(dict.lookupOrDefault<scalar>("a_z", scalar(0))),
    
    b_x(dict.lookupOrDefault<scalar>("b_x", scalar(1))),
    b_y(dict.lookupOrDefault<scalar>("b_y", scalar(1))),
    b_z(dict.lookupOrDefault<scalar>("b_z", scalar(1))),
    
    k_x(dict.lookupOrDefault<scalar>("k_x", scalar(0))),
    k_y(dict.lookupOrDefault<scalar>("k_y", scalar(0))),
    k_z(dict.lookupOrDefault<scalar>("k_z", scalar(0))),
    
    L_x(dict.lookupOrDefault<scalar>("L_x", scalar(1))),
    L_y(dict.lookupOrDefault<scalar>("L_y", scalar(1))),
    L_z(dict.lookupOrDefault<scalar>("L_z", scalar(1))),
    
    phi_x(dict.lookupOrDefault<scalar>("phi_x", scalar(0))),
    phi_y(dict.lookupOrDefault<scalar>("phi_y", scalar(0))),
    phi_z(dict.lookupOrDefault<scalar>("phi_z", scalar(0))),
    
    funcMax(dict.lookupOrDefault<scalar>("funcMax", -GREAT)),
    funcMin(dict.lookupOrDefault<scalar>("funcMin", -GREAT))
{};

scalar sinusoidalTracerField::tracerAt
(
        const point& p,
        const Time& t
) const
{
    scalar x = p.x();
    scalar y = p.y();
    scalar z = p.z();
    scalar func = mag * (a_x + b_x*Foam::cos(2*M_PI*k_x*x/L_x - phi_x))
                      * (a_y + b_y*Foam::cos(2*M_PI*k_y*y/L_y - phi_y))
                      * (a_z + b_z*Foam::cos(2*M_PI*k_z*z/L_z - phi_z));
    return max(min(func, funcMax), funcMin);
}

