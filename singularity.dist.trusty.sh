#!/bin/bash
set -e

VERSION=14.04
CODENAME=trusty

PACKAGE_VERSION=$(./singularity.bootstrap.sh $VERSION $CODENAME)
./dist.sh $CODENAME
sudo singularity exec -e -w $CODENAME.img apt-get update -qq
sudo singularity exec -e -w $CODENAME.img apt-get install atmosfoam-tools=$PACKAGE_VERSION -y --allow-unauthenticated --no-install-recommends
git clean -xfd .
