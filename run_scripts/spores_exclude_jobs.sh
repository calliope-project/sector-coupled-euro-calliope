#!/bin/sh

case "$1" in

    1) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs roof_mounted_pv ;;
    2) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs open_field_pv ;;
    3) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs wind_onshore_monopoly ;;
    4) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs wind_onshore_competing ;;
    5) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs wind_offshore ;;
    6) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs biofuel_supply ;;
    7) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs electrolysis ;;
    8) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs dac ;;
    9) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs heavy_transport_ice ;;
    10) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs heavy_transport_ev ;;
    11) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs light_transport_ice ;;
    12) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs light_transport_ev ;;
    13) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs hydrogen_storage ;;
    14) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs battery ;;
    15) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs dc_ohl_transmission ;;
    16) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs dc_subsea_transmission ;;
    17) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs dc_underground_transmission ;;
    18) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs ac_ohl_transmission ;;
    19) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs ac_ohl_mountain_transmission ;;
    20) python "${3}/spores_exclude.py" $2 --n_spores 10 --tech_groups heat_storage_small ;;
    21) python "${3}/spores_exclude.py" $2 --n_spores 10 --tech_groups heat_storage_big ;;
    22) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs wind_onshore_monopoly,wind_onshore_competing ;;
    23) python "${3}/spores_exclude.py" $2 --n_spores 10 --tech_groups wind ;;
    24) python "${3}/spores_exclude.py" $2 --n_spores 10 --tech_groups pv ;;
    25) python "${3}/spores_exclude.py" $2 --n_spores 10 --techs electrolysis,dac ;;
    26) python "${3}/spores_exclude.py" $2 --n_spores 10 --tech_groups transport_ice ;;
    27) python "${3}/spores_exclude.py" $2 --n_spores 10 --tech_groups transport_ev ;;
    28) python "${3}/spores_exclude.py" $2 --n_spores 10 --tech_groups transmission ;;
    29) python "${3}/spores_exclude.py" $2 --n_spores 10 --tech_groups storage ;;

esac
