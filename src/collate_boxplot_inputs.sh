#!/bin/bash

basedir="$1"
Nc="$2"
Representations="F AS S"
massChannels="vector axialvector scalar tensor axialtensor"
decayChannels="pseudoscalar vector axialvector"

cd "${basedir}"
for representation in ${Representations}
do
  masses_filename="${representation}_masses.txt"
  decayconsts_filename="${representation}_decayconsts.txt"

  cd "${representation}"
  [[ -f "${masses_filename}" ]] && rm "${masses_filename}"
  [[ -f "${decayconsts_filename}" ]] && rm "${decayconsts_filename}"
  for channel in ${massChannels}
  do
    cat ${channel}_masses_${representation}_Sp${Nc}.dat 2>/dev/null | head -1 | awk '{print $1} END {if(NR==0) {print "--"}}' >> "${masses_filename}"
    cat ${channel}_masses_${representation}_Sp${Nc}.dat 2>/dev/null | head -1 | awk '{print $2} END {if(NR==0) {print "--"}}' >> "${masses_filename}"
  done


  for channel in ${decayChannels}
  do
    cat ${channel}_decayconsts_${representation}_Sp${Nc}.dat 2>/dev/null | head -1 | awk '{print $1} END {if(NR==0) {print "--"}}' >> "${decayconsts_filename}"
    cat ${channel}_decayconsts_${representation}_Sp${Nc}.dat 2>/dev/null | head -1 | awk '{print $2} END {if(NR==0) {print "--"}}' >> "${decayconsts_filename}"
  done
  cd ..
done
