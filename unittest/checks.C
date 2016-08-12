#include "checks.H"

void check(Foam::vector actual, Foam::vector expected)
{
    CHECK(actual.x() == approx(expected.x()));
    CHECK(actual.y() == approx(expected.y()));
    CHECK(actual.z() == approx(expected.z()));
}

void check(Foam::tensor actual, Foam::tensor expected)
{
    CHECK(actual.xx() == approx(expected.xx()));
    CHECK(actual.xy() == approx(expected.xy()));
    CHECK(actual.xz() == approx(expected.xz()));

    CHECK(actual.yx() == approx(expected.yx()));
    CHECK(actual.yy() == approx(expected.yy()));
    CHECK(actual.yz() == approx(expected.yz()));

    CHECK(actual.zx() == approx(expected.zx()));
    CHECK(actual.zy() == approx(expected.zy()));
    CHECK(actual.zz() == approx(expected.zz()));
}
