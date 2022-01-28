#!/usr/bin/bash

set -e 

# grep ERROR LOG | grep ready | sed -e 's/), Node(/\n/g' -e 's/.*(//'  | grep runs | parallel -j 8 bash ~/RGP/scripts/parallel_falcon_run_daligner.sh

direc=$@

if [ ! -e "${direc}/run.sh.done" ]; then
    if [ -d "${direc}" ]; then
        cd ${direc};
        echo ${direc};
        bash user_script.sh && touch run.sh.done
    fi
fi

