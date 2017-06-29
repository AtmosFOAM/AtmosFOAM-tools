#!/bin/bash
set -e

display_usage() {
	echo -e "Usage: dist.sh <codename>\n"
}

if [ $# -lt 1 ]
then
	display_usage
	exit 1
fi

VERSION=$(date +"%Y%m%d%H%M%S")
CODENAME=$1

# create a skeletal debian/changelog
SINGULARITYENV_DEBFULLNAME=$DEBFULLNAME \
SINGULARITYENV_DEBEMAIL=$DEBEMAIL \
singularity exec -e $CODENAME.img dch --create --package "atmosfoam-tools" --distribution $CODENAME --newversion=$VERSION "Build for $(git rev-parse HEAD)"

# create an amd64 binary package
# using the Makefile 'debian/rules'
SINGULARITYENV_DEBFULLNAME=$DEBFULLNAME \
SINGULARITYENV_DEBEMAIL=$DEBEMAIL \
SINGULARITYENV_AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID \
SINGULARITYENV_AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY \
singularity exec -e $CODENAME.img debuild -i -us -uc -b

# upload to the debian apt repository located in the
# Amazon S3 bucket.  deb-s3 expects AWS_ACCESS_KEY_ID and
# AWS_SECRET_ACCESS_KEY environment variables to be set.
deb-s3 upload --bucket atmosfoam-apt --codename=$CODENAME --component=dev ../atmosfoam-tools_${VERSION}_amd64.deb

# test the new .deb can be installed
sudo apt-get install atmosfoam-tools=$VERSION -y --allow-unauthenticated --no-install-recommends

git clean -xfd .
