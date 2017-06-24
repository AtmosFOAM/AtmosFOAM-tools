#!/bin/bash
set -e
VERSION=$(date +"%Y%m%d%H%M%S")
CODENAME="trusty"

# create a skeletal debian/changelog
dch --create --package "atmosfoam-tools" --distribution $CODENAME --newversion=$VERSION "Build for $(git rev-parse HEAD)"

# create an amd64 binary package
# using the Makefile 'debian/rules'
debuild -i -us -uc -b

# upload to the debian apt repository located in the
# Amazon S3 bucket.  deb-s3 expects AWS_ACCESS_KEY_ID and
# AWS_SECRET_ACCESS_KEY environment variables to be set.
deb-s3 upload --bucket atmosfoam-apt --codename=$CODENAME ../atmosfoam-tools_${VERSION}_amd64.deb
