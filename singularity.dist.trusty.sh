#!/bin/bash
set -e

VERSION=14.04
CODENAME=trusty

./singularity.bootstrap.sh $VERSION $CODENAME
./dist.sh $CODENAME
sudo singularity exec -e -w $CODENAME.img apt-get update -qq
sudo singularity exec -e -w $CODENAME.img apt-get install atmosfoam-tools=$VERSION -y --allow-unauthenticated --no-install-recommends
git clean -xfd .
