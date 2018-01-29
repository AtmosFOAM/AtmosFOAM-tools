#!/bin/bash
set -e

VERSION=17.10
CODENAME=artful

PACKAGE_VERSION=$(./singularity.bootstrap.sh $VERSION $CODENAME)
./dist.sh $CODENAME
sudo singularity exec -e -w $CODENAME.img apt-get update -qq --allow-insecure-repositories
sudo singularity exec -e -w $CODENAME.img apt-get install atmosfoam-tools=$PACKAGE_VERSION -y --allow-unauthenticated --no-install-recommends
git clean -xfd .
