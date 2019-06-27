#!/bin/bash 
HN=`hostname`
if [ "$HN" == "rf73u" ]; then
    echo "Running docker on the rf73u machine"
    docker run -it -v /home/rf73/neiproj/public:/root/public \
                    -v /home/rf73/rafgh/GWAStutorial:/root/dockspace \
                    bioc0 /bin/bash
elif [ "$HN" == "biotime" ]; then
    echo "Running docker on the rf73u machine"
    docker run -it -v /home/nutria/neiproj/public:/root/public \
                    -v /home/nutria/rafgh/GWAStutorial:/root/dockspace \
                    bioc0 /bin/bash
fi
