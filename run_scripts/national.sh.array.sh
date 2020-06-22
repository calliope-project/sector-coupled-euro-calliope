#!/bin/sh

function process_case () {
    case "$1" in

    1) calliope run ../build/model/national/model.yaml --scenario directional-rooftop-pv --save_netcdf out_1_directional-rooftop-pv.nc --save_plots plots_1_directional-rooftop-pv.html ;;

    esac
}

if [[ $# -eq 0 ]] ; then
    echo No parameter given, running all runs sequentially...
    for i in $(seq 1 1); do process_case $i; done
else
    echo Running run $1
    process_case $1
fi
