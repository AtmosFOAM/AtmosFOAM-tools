#!/bin/bash
set -e

display_usage() {
	echo -e "Usage: singularity.bootstrap.sh <version> <codename>\n"
}

if [ $# -le 1 ]
then
	display_usage
	exit 1
fi

export VERSION=14.04
export CODENAME=trusty

./singularity.bootstrap.sh $VERSION $CODENAME
./dist.sh $CODENAME
sudo singularity exec -e -w $CODENAME.img apt-get update -qq
sudo singularity exec -e -w $CODENAME.img apt-get install atmosfoam-tools=$VERSION -y --allow-unauthenticated --no-install-recommends
git clean -xfd .
