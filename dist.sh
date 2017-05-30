#!/bin/bash
set -e
VERSION=$(date +"%Y%m%d%H%M%S")

dch --create --package "atmosfoam-tools" --distribution "trusty" --newversion=$VERSION "Build for $(git rev-parse HEAD)"
debuild -i -us -uc -b
deb-s3 upload --bucket atmosfoam-apt --codename=dev ../atmosfoam-tools_${VERSION}_amd64.deb
