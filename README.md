# AtmosFOAM-tools
[![DOI](https://zenodo.org/badge/64257768.svg)](https://zenodo.org/badge/latestdoi/64257768)

A repository  which contains generic libraries and utilities for
    https://github.com/hertzsprung/AtmosFOAM
    and
    https://github.com/pbrowne/AMMM
    Upgraded for OpenFOAM-4.x or OpenFOAM-dev:
    https://github.com/OpenFOAM/OpenFOAM-dev
    46d69e1 commit a6056b7329d0ac93403beef3a1048e97dbf010b1


## Installation
First, install [OpenFOAM 4 or dev](http://www.openfoam.org/download/).

Second, set the variables the ATMOSFOAM_TOOLS_SRC and GMTU in the .bashrc file:

export ATMOSFOAM_TOOLS_SRC=$HOME/$WM_PROJECT/$USER-$WM_PROJECT_VERSION/AtmosFOAM-tools/src
export GMTU=$HOME/$WM_PROJECT/$USER-$WM_PROJECT_VERSION/AtmosFOAM-tools/gmtUser

Third compile AtmosFOAM-tools using Allwmake

