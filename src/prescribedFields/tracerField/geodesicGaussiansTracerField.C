#include "geodesicGaussiansTracerField.H"
#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(geodesicGaussiansTracerField, 0);
addToRunTimeSelectionTable(tracerField, geodesicGaussiansTracerField, dict);

geodesicGaussiansTracerField::geodesicGaussiansTracerField
(
    const dictionary& dict,
    const advectable& velocityField
)
:
    tracerField(velocityField),
    R("radius", dimLength, readScalar(dict.lookup("radius"))),
    hmax("hmax", dimless, readScalar(dict.lookup("hmax"))),
    b("b", dimless, readScalar(dict.lookup("b")))
{};

scalar geodesicGaussiansTracerField::tracerAt(const point& p, const Time& t) const
{
    scalar lon1 = 5.0*M_PI/6.0;
    scalar lat1 = 0;

    scalar lon2 = 7.0*M_PI/6.0;
	scalar lat2 = 0;

	const dimensionedVector dimensionedP("p", dimLength, p);

    const dimensionedVector centre1("centre1", dimLength, point(
            R.value() * Foam::cos(lat1) * Foam::cos(lon1),
            R.value() * Foam::cos(lat1) * Foam::sin(lon1),
            R.value() * Foam::sin(lat1)
    ));

    const dimensionedVector centre2("centre2", dimLength, point(
            R.value() * Foam::cos(lat2) * Foam::cos(lon2),
            R.value() * Foam::cos(lat2) * Foam::sin(lon2),
            R.value() * Foam::sin(lat2)
    ));

    return (hmax * exp(-b*magSqr(dimensionedP - centre1)/sqr(R)) +
        hmax * exp(-b*magSqr(dimensionedP - centre2)/sqr(R))).value();
}
