#include "catch.hpp"

#include "sphericalVector.H"
#include "checks.H"

using namespace Foam;

TEST_CASE("zeroCartesianAtCentreOfEarth")
{
    sphericalVector velocity(3, 4, -2);
    sphericalVector p(0, 0, 0);

    check(velocity.toCartesian(p), vector(0, 0, 0));
}

TEST_CASE("cartesianVelocityAtEquatorOnPrimeMeridian")
{
    sphericalVector velocity(3, 4, -2);
    sphericalVector p(5, 0, 0);

    check(velocity.toCartesian(p), vector(-2, 3, 4));
}

TEST_CASE("cartesianVelocityAtNorthPole")
{
    sphericalVector velocity(3, 4, -2);
    sphericalVector p(0, 0, 5);

    check(velocity.toCartesian(p), vector(4, -3, -2));
}

TEST_CASE("cartesianVelocityAtSouthPole")
{
    sphericalVector velocity(3, 4, -2);
    sphericalVector p(0, 0, -5);

    check(velocity.toCartesian(p), vector(-4, -3, 2));
}

TEST_CASE("cartesianVelocityAt45N")
{
    sphericalVector velocity(0, 1, 0);
    sphericalVector p(1, 0, 1);

    check(velocity.toCartesian(p), vector(-Foam::sqrt(0.5), 0, Foam::sqrt(0.5)));
}

TEST_CASE("unitTensorAtEquatorOnPrimeMeridian")
{
    sphericalVector p(5, 0, 0);

    check(p.unitTensor(), tensor(vector(0, 1, 0), vector(0, 0, 1), vector(1, 0, 0)));
}

