#!/bin/bash
set -e
#TODO: make sure src libraries can see headers from other src libraries
#wmake -j -a src
#wmake -j all applications/utilities/postProcessing
wmake -j applications/utilities/postProcessing/globalSum
wmake -j applications/utilities/postProcessing/sumFields
