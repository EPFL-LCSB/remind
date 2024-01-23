#!/bin/sh


#--cpuset-cpus 6-10
docker run 	--rm -it 			\
		-v $(pwd)/work:/home/remind/work \
		-v $(pwd)/..:/remind	\
		remind_docker



#docker run 	--rm -it 			\
#		-v $(pwd)/work:/home/coptik/work \
#		-v $(pwd)/..:/coptik	\
#		-v $(pwd)/..:/openbread	\
#		coptik_docker

