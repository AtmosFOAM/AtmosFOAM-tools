# AtmosFOAM-tools
AtmosFOAM-tools contains generic libraries and utilities that support atmospheric simulations with [OpenFOAM](https://openfoam.org/).  These generic tools can be combined with [AtmosFOAM](https://github.com/AtmosFOAM/AtmosFOAM) and [AMMM](https://github.com/AtmosFOAM/AMMM) repositories.

## Source installation
1. Install a recent version of [openfoam-dev](http://www.openfoam.org/download/).
2. Export environment variables `ATMOSFOAM_TOOLS_SRC` and `GMTU` in your `~/.bashrc` file:

       export ATMOSFOAM_TOOLS_SRC=/path/to/AtmosFOAM-tools/src
       export GMTU=/path/to/AtmosFOAM-tools/gmtUser
      
3. Compile AtmosFOAM-tools using `./Allwmake`
4. gmtFoam plots OpenFOAM output using [gmt](http://gmt.soest.hawaii.edu/) so to use gmtFoam you will also need to install gmt.
