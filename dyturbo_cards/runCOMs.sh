#!/bin/bash
oshname="runAll.sh"

if [ -f $oshname ]; then
    rm -rf $oshname
fi

# pp
declare -a channels=(Wp Wm)
declare -a coms=(2 3 4 5 5.02 6 7 8 9 10 13 14 20)

for chan in "${channels[@]}"; do
    fileorg="cards/${chan}-13tev_inc.in"
    for com in "${coms[@]}"; do
        out="${chan}-${com}tev_inc_pp"
        ifile="cards/coms/${out}.in"
        sed "s/sroot        = 13e3/sroot        = ${com}e3/g" $fileorg > $ifile
        echo "./bin/dyturbo ${ifile} >& logs/${out}.log" >> $oshname
    done
done

# ppbar
declare -a channels=(Wp Wm Z)
declare -a coms=(0.5 0.7 1 1.4 1.96 2 2.5 3)

for chan in ${channels[@]}; do
    fileorg="cards/${chan}-13tev_inc.in"
    for com in ${coms[@]}; do
        out="${chan}-${com}tev_inc_ppbar"
        ifile="cards/coms/${out}.in"
        sed -e "s/sroot        = 13e3/sroot        = ${com}e3/g" -e "s/ih2          = 1/ih2          = -1/g" $fileorg > $ifile
        echo "./bin/dyturbo ${ifile} >& logs/${out}.log" >> $oshname
    done
done

chmod +x $oshname