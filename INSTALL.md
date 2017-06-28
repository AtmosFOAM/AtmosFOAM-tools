# Debian packaging
The continuous integration system, [Travis CI](https://travis-ci.org), is triggered when a commit is pushed to the repository.  Separate binary `.deb` packages are created for different Ubuntu releases.  The `.travis.yml` configuration coordinates the 4-step process:

1. [Singularity](http://singularity.lbl.gov/) is installed to host the different Ubuntu releases
2. [deb-s3](https://github.com/krobertson/deb-s3) is installed, allowing `.deb` packages to be uploaded to an APT package repository hosted in an [Amazon S3 bucket](http://docs.aws.amazon.com/AmazonS3/latest/dev/UsingBucket.html)
3. For each Ubuntu release, a Singularity container is bootstrapped and a `.deb` package is built within the container
4. On the Travis CI host, the `.deb` package is uploaded to the APT package repository using deb-s3

## Debian packaging files

* `debian/control` describes the package and its dependencies
* `debian/rules` is a Makefile that is used by `debuild` to prepare the package file structure, with files installed into `debian/atmosfoam-tools/usr/`

## Configuring the Amazon S3 bucket

1. Create an S3 bucket called `atmosfoam-apt` [in the AWS S3 console](https://console.aws.amazon.com/s3)
2. Enable static web hosting
3. Create an IAM user with `AmazonS3FullAccess` [in the AWS IAM console](https://console.aws.amazon.com/iam/home)
4. Create an access key pair for the new IAM user
5. Modify [.travis.yml](https://github.com/AtmosFOAM/AtmosFOAM-tools/blob/master/.travis.yml) to use the new `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY`, using [Travis CI encryption](https://docs.travis-ci.com/user/environment-variables/#Defining-encrypted-variables-in-.travis.yml) to secure the secret key
