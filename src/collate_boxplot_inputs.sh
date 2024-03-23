#!/bin/bash

basedir="$1"
Nc="$2"
Representations="F AS S"
massChannels="vector axialvector scalar tensor axialtensor"
decayChannels="pseudoscalar vector axialvector"

cd "${basedir}"
for representation in ${Representations}
do
  cd $i
  [[ -f ${representation}_decayconsts.txt ]] && rm ${representation}_decayconsts.txt
  [[ -f ${representation}_masses.txt ]] && rm ${representation}_masses.txt
  for channel in ${massChannels}
  do
    cat mass_${channel}_${representation}_Sp${Nc}.dat | head -1 | awk '{print $1}' >> ${representation}_masses.txt;
    cat mass_${channel}_${representation}_Sp${Nc}.dat | head -1 | awk '{print $2}' >> ${representation}_masses.txt;
  done
  for channel in ${decayChannels}
  do
    cat decayconst_${channel}_${representation}_Sp${Nc}.dat | head -1 | awk '{print $1}' >> ${representation}_decayconsts.txt;
    cat decayconst_${channel}_${representation}_Sp${Nc}.dat | head -1 | awk '{print $2}' >> ${representation}_decayconsts.txt;
  done
  cd ..
done
