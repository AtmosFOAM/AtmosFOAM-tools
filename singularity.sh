#!/bin/bash
sudo singularity create --size 2048 trusty.img
sudo singularity bootstrap trusty.img Singularity

SINGULARITYENV_DEBFULLNAME=$DEBFULLNAME \
SINGULARITYENV_DEBEMAIL=$DEBEMAIL \
SINGULARITYENV_AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID \
SINGULARITYENV_AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY \
	singularity exec -e trusty.img ./dist.sh
