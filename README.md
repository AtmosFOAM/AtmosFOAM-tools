# AtmosFOAM-tools
AtmosFOAM-tools contains generic libraries and utilities that support atmospheric simulations with [OpenFOAM](https://openfoam.org/).  These generic tools can be combined with [AtmosFOAM](https://github.com/AtmosFOAM/AtmosFOAM) and [AMMM](https://github.com/AtmosFOAM/AMMM) repositories.

## Source installation
* Install a recent version of [openfoam-dev](http://www.openfoam.org/download/).
* Install libgdal. On ubuntu this is done with `apt-get install libgdal-dev`
* Go to directory
cd $WM_PROJECT_USER_DIR
and download AtmosFOAM-tools using:
git clone https://github.com/AtmosFOAM/AtmosFOAM-tools.git
* Export environment variables `ATMOSFOAM_TOOLS_SRC` and `GMTU` in your `~/.bashrc` file:

       export ATMOSFOAM_TOOLS_SRC=/path/to/AtmosFOAM-tools/src
       export GMTU=/path/to/AtmosFOAM-tools/gmtUser
      
* Compile AtmosFOAM-tools:
cd AtmosFOAM-tools
./Allwmake

