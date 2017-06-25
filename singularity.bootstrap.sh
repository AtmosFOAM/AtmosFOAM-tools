#!/bin/bash
sudo singularity create --size 2048 trusty.img
sudo singularity bootstrap trusty.img Singularity
