#!/bin/bash

for A in Lp Wb Pd Br Gr Ma Om BF LF Hm
do
    for B in Lp Wb Pd Br Gr Ma Om BF LF Hm
    do
        if [[ $A < $B ]] ; then
            echo """
sh 01.create_scripts_for_SAMap_one2one.sh $A $B
""" >run.${A}.${B}.sh
chmod u+x run.${A}.${B}.sh
        fi
    done
done
