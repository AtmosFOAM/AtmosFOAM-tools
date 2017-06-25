Bootstrap:docker
From:ubuntu:14.04

%post
	apt-get update -qq
	apt-get install wget software-properties-common apt-transport-https -y --no-install-recommends
	sh -c "wget -O - http://dl.openfoam.org/gpg.key | apt-key add -"
	add-apt-repository "http://dl.openfoam.org/ubuntu dev" -y
	apt-get -qq update
	apt-get install git openfoam-dev devscripts debhelper ruby-dev -y
	gem install deb-s3
