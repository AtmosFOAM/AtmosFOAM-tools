# Debian packaging
The continuous integration system, [Travis-CI](https://travis-ci.org), is triggered when a commit is pushed to the repository.  The `.travis.yml` configuration installs the necessary build dependencies before invoking `dist.sh` to compile and package AtmosFOAM-tools.  This script creates a binary `.deb` package that is uploaded to an APT package repository hosted in an [Amazon S3 bucket](http://docs.aws.amazon.com/AmazonS3/latest/dev/UsingBucket.html).

* `debian/control` describes the package and its dependencies
* `debian/rules` is a Makefile that is used by `debuild` to prepare the package file structure, with files installed into `debian/atmosfoam-tools/usr/`
