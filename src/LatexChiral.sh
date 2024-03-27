#!/bin/bash

#set -eux

N=$1
basedir="$2"
output_file="${PWD}/$3"

cd "${basedir}"

Reps="F AS S"
[ -f "${output_file}" ] && rm "${output_file}"

echo "\\begin{tabular}{|c|c|c|c|}" >> "${output_file}"
echo "\\hline" >> "${output_file}"
echo "\\multicolumn{4}{|c|}{\$Sp(${N})\$}\\\\" >> "${output_file}"
echo "\\hline" >> "${output_file}"
echo "Representation & Channel & Chiral limit & \$\\chi^2/{\\rm d.o.f.}\$\\\\" >> "${output_file}"
echo "\\hline" >> "${output_file}"

for rep in ${Reps}
do
    channel="pseudoscalar"
    if [ "${rep}" == "F" ]
    then
        representation="Fundamental"
        echo -n "\\multirow{8}{*}{${representation}} & \$\\hat{f}_{\\rm PS}\$ & " >> "${output_file}"
    elif [ "${rep}" == "AS" ]
    then
        representation="Antisymmetric"
        echo -n "\\multirow{8}{*}{${representation}} & \$\\hat{f}_{\\rm ps}\$ & " >> "${output_file}"
    else
        representation="Symmetric"
        echo -n "\\multirow{8}{*}{${representation}} & \$\\hat{f}_{\\mathcal{PS}}\$ & " >> "${output_file}"
    fi

    o2=$(cat ${rep}/${channel}_decayconsts_${rep}_Sp${N}.dat | head -1 | awk '{print $1}')
    oError2=$(cat ${rep}/${channel}_decayconsts_${rep}_Sp${N}.dat | head -1 | awk '{print $2}')

    if [ ${o2} == "--" ]
    then
        echo "-- & -- \\\\" >> "${output_file}"
    else
        o=$(echo ${o2} | awk '{print sqrt($1)}')
        oError=$(echo ${oError2} ${o} | awk '{print $1 /2 /$2}')
        num=1
        tmp=${oError}
        while (( $(echo "$tmp < 1.0" | bc -l) ))
        do
            tmp=$(echo $tmp | awk '{print $1 * 10}')
            num=`expr $(($num +1 ))`
        done
        oErrRounded=$(printf "%.${num}f" ${oError})
        chi2=$(cat ${rep}/${channel}_decayconsts_${rep}_Sp${N}.dat | tail -1 | awk '{print $1}' | xargs printf "%1.2f")

        echo -n $(printf "%.${num}f" $o) >> "${output_file}"
        echo -n "(" >> "${output_file}"
        echo -n $(printf "%s" "${oErrRounded}" | cut -c `expr $((${num}+1))`-) >> "${output_file}"
        echo -n ") & " >> "${output_file}"
        if (( $(echo "${chi2} > 3.0" |bc -l) ));
        then
            echo -n "\\textcolor{red}{${chi2}}" >> "${output_file}"
        else
            echo -n "${chi2}" >> "${output_file}"
        fi
        echo "\\\\" >> "${output_file}"
    fi

    echo "\\cline{2-4}" >> "${output_file}"

    channel="vector"

    if [ "${rep}" == "F" ]
    then
        echo -n " & \$\\hat{f}_{\\rm V}\$ & " >> "${output_file}"
    elif [ "${rep}" == "AS" ]
    then
        echo -n " & \$\\hat{f}_{\\rm v}\$ & " >> "${output_file}"
    else
        echo -n " & \$\\hat{f}_{\\mathcal{V}}\$ & " >> "${output_file}"
    fi

    o2=$(cat ${rep}/${channel}_decayconsts_${rep}_Sp${N}.dat | head -1 | awk '{print $1}')
    oError2=$(cat ${rep}/${channel}_decayconsts_${rep}_Sp${N}.dat | head -1 | awk '{print $2}')

    if [ ${o2} == "--" ]
    then
        echo "-- & -- \\\\" >> "${output_file}"
    else
        o=$(echo ${o2} | awk '{print sqrt($1)}')
        oError=$(echo ${oError2} ${o} | awk '{print $1 /2 /$2}')
        num=1
        tmp=${oError}
        while (( $(echo "$tmp < 1.0" | bc -l) ))
        do
            tmp=$(echo $tmp | awk '{print $1 * 10}')
            num=`expr $(($num +1 ))`
        done
        oErrRounded=$(printf "%.${num}f" ${oError})
        chi2=$(cat ${rep}/${channel}_decayconsts_${rep}_Sp${N}.dat | tail -1 | awk '{print $1}' | xargs printf "%1.2f")

        echo -n $(printf "%.${num}f" $o) >> "${output_file}"
        echo -n "(" >> "${output_file}"
        echo -n $(printf "%s" "${oErrRounded}" | cut -c `expr $((${num}+1))`-) >> "${output_file}"
        echo -n ") & " >> "${output_file}"
        if (( $(echo "${chi2} > 3.0" |bc -l) ));
        then
            echo -n "\\textcolor{red}{${chi2}}" >> "${output_file}"
        else
            echo -n "${chi2}" >> "${output_file}"
        fi
        echo "\\\\" >> "${output_file}"
    fi

    echo "\\cline{2-4}" >> "${output_file}"

    channel="axialvector"

    if [ "${rep}" == "F" ]
    then
        echo -n " & \$\\hat{f}_{\\rm AV}\$ & " >> "${output_file}"
    elif [ "${rep}" == "AS" ]
    then
        echo -n " & \$\\hat{f}_{\\rm av}\$ & " >> "${output_file}"
    else
        echo -n " & \$\\hat{f}_{\\mathcal{AV}}\$ & " >> "${output_file}"
    fi

    o2=$(cat ${rep}/${channel}_decayconsts_${rep}_Sp${N}.dat | head -1 | awk '{print $1}')
    oError2=$(cat ${rep}/${channel}_decayconsts_${rep}_Sp${N}.dat | head -1 | awk '{print $2}')

    if [ ${o2} == "--" ]
    then
        echo "-- & -- \\\\" >> "${output_file}"
    else
        o=$(echo ${o2} | awk '{print sqrt($1)}')
        oError=$(echo ${oError2} ${o} | awk '{print $1 /2 /$2}')
        num=1
        tmp=${oError}
        while (( $(echo "$tmp < 1.0" | bc -l) ))
        do
            tmp=$(echo $tmp | awk '{print $1 * 10}')
            num=`expr $(($num +1 ))`
        done
        oErrRounded=$(printf "%.${num}f" ${oError})
        chi2=$(cat ${rep}/${channel}_decayconsts_${rep}_Sp${N}.dat | tail -1 | awk '{print $1}' | xargs printf "%1.2f")

        echo -n $(printf "%.${num}f" $o) >> "${output_file}"
        echo -n "(" >> "${output_file}"
        echo -n $(printf "%s" "${oErrRounded}" | cut -c `expr $((${num}+1))`-) >> "${output_file}"
        echo -n ") & " >> "${output_file}"
        if (( $(echo "${chi2} > 3.0" |bc -l) ));
        then
            echo -n "\\textcolor{red}{${chi2}}" >> "${output_file}"
        else
            echo -n "${chi2}" >> "${output_file}"
        fi
        echo "\\\\" >> "${output_file}"
    fi

    echo "\\cline{2-4}" >> "${output_file}"

    channel="vector"

    if [ "${rep}" == "F" ]
    then
        echo -n " & \$\\hat{m}_{\\rm V}\$ & " >> "${output_file}"
    elif [ "${rep}" == "AS" ]
    then
        echo -n " & \$\\hat{m}_{\\rm v}\$ & " >> "${output_file}"
    else
        echo -n " & \$\\hat{m}_{\\mathcal{V}}\$ & " >> "${output_file}"
    fi

    o2=$(cat ${rep}/${channel}_masses_${rep}_Sp${N}.dat | head -1 | awk '{print $1}')
    oError2=$(cat ${rep}/${channel}_masses_${rep}_Sp${N}.dat | head -1 | awk '{print $2}')

    if [ ${o2} == "--" ]
    then
        echo "-- & -- \\\\" >> "${output_file}"
    else
        o=$(echo ${o2} | awk '{print sqrt($1)}')
        oError=$(echo ${oError2} ${o} | awk '{print $1 /2 /$2}')
        num=1
        tmp=${oError}
        while (( $(echo "$tmp < 1.0" | bc -l) ))
        do
            tmp=$(echo $tmp | awk '{print $1 * 10}')
            num=`expr $(($num +1 ))`
        done
        oErrRounded=$(printf "%.${num}f" ${oError})
        chi2=$(cat ${rep}/${channel}_masses_${rep}_Sp${N}.dat | tail -1 | awk '{print $1}' | xargs printf "%1.2f")

        echo -n $(printf "%.${num}f" $o) >> "${output_file}"
        echo -n "(" >> "${output_file}"
        echo -n $(printf "%s" "${oErrRounded}" | cut -c `expr $((${num}+1))`-) >> "${output_file}"
        echo -n ") & " >> "${output_file}"
        if (( $(echo "${chi2} > 3.0" |bc -l) ));
        then
            echo -n "\\textcolor{red}{${chi2}}" >> "${output_file}"
        else
            echo -n "${chi2}" >> "${output_file}"
        fi
        echo "\\\\" >> "${output_file}"
    fi

    echo "\\cline{2-4}" >> "${output_file}"

    channel="axialvector"

    if [ "${rep}" == "F" ]
    then
        echo -n " & \$\\hat{m}_{\\rm AV}\$ & " >> "${output_file}"
    elif [ "${rep}" == "AS" ]
    then
        echo -n " & \$\\hat{m}_{\\rm av}\$ & " >> "${output_file}"
    else
        echo -n " & \$\\hat{m}_{\\mathcal{AV}}\$ & " >> "${output_file}"
    fi

    o2=$(cat ${rep}/${channel}_masses_${rep}_Sp${N}.dat | head -1 | awk '{print $1}')
    oError2=$(cat ${rep}/${channel}_masses_${rep}_Sp${N}.dat | head -1 | awk '{print $2}')

    if [ ${o2} == "--" ]
    then
        echo "-- & -- \\\\" >> "${output_file}"
    else
        o=$(echo ${o2} | awk '{print sqrt($1)}')
        oError=$(echo ${oError2} ${o} | awk '{print $1 /2 /$2}')
        num=1
        tmp=${oError}
        while (( $(echo "$tmp < 1.0" | bc -l) ))
        do
            tmp=$(echo $tmp | awk '{print $1 * 10}')
            num=`expr $(($num +1 ))`
        done
        oErrRounded=$(printf "%.${num}f" ${oError})
        chi2=$(cat ${rep}/${channel}_masses_${rep}_Sp${N}.dat | tail -1 | awk '{print $1}' | xargs printf "%1.2f")

        echo -n $(printf "%.${num}f" $o) >> "${output_file}"
        echo -n "(" >> "${output_file}"
        echo -n $(printf "%s" "${oErrRounded}" | cut -c `expr $((${num}+1))`-) >> "${output_file}"
        echo -n ") & " >> "${output_file}"
        if (( $(echo "${chi2} > 3.0" |bc -l) ));
        then
            echo -n "\\textcolor{red}{${chi2}}" >> "${output_file}"
        else
            echo -n ${chi2} >> "${output_file}"
        fi
        echo "\\\\" >> "${output_file}"

    fi
    echo "\\cline{2-4}" >> "${output_file}"

    channel="scalar"

    if [ "${rep}" == "F" ]
    then
        echo -n " & \$\\hat{m}_{\\rm S}\$ & " >> "${output_file}"
    elif [ "${rep}" == "AS" ]
    then
        echo -n " & \$\\hat{m}_{\\rm s}\$ & " >> "${output_file}"
    else
        echo -n " & \$\\hat{m}_{\\mathcal{S}}\$ & " >> "${output_file}"
    fi

    o2=$(cat ${rep}/${channel}_masses_${rep}_Sp${N}.dat | head -1 | awk '{print $1}')
    oError2=$(cat ${rep}/${channel}_masses_${rep}_Sp${N}.dat | head -1 | awk '{print $2}')

    if [ ${o2} == "--" ]
    then
        echo "-- & -- \\\\" >> "${output_file}"
    else
        o=$(echo ${o2} | awk '{print sqrt($1)}')
        oError=$(echo ${oError2} ${o} | awk '{print $1 /2 /$2}')
        num=1
        tmp=${oError}
        while (( $(echo "$tmp < 1.0" | bc -l) ))
        do
            tmp=$(echo $tmp | awk '{print $1 * 10}')
            num=`expr $(($num +1 ))`
        done
        oErrRounded=$(printf "%.${num}f" ${oError})
        chi2=$(cat ${rep}/${channel}_masses_${rep}_Sp${N}.dat | tail -1 | awk '{print $1}' | xargs printf "%1.2f")

        echo -n $(printf "%.${num}f" $o) >> "${output_file}"
        echo -n "(" >> "${output_file}"
        echo -n $(printf "%s" "${oErrRounded}" | cut -c `expr $((${num}+1))`-) >> "${output_file}"
        echo -n ") & " >> "${output_file}"
        if (( $(echo "${chi2} > 3.0" |bc -l) ));
        then
            echo -n "\\textcolor{red}{${chi2}}" >> "${output_file}"
        else
            echo -n "${chi2}" >> "${output_file}"
        fi
        echo "\\\\" >> "${output_file}"
    fi

    echo "\\cline{2-4}" >> "${output_file}"

    channel="tensor"

    if [ "${rep}" == "F" ]
    then
        echo -n " & \$\\hat{m}_{\\rm T}\$ & " >> "${output_file}"
    elif [ "${rep}" == "AS" ]
    then
        echo -n " & \$\\hat{m}_{\\rm t}\$ & " >> "${output_file}"
    else
        echo -n " & \$\\hat{m}_{\\mathcal{T}}\$ & " >> "${output_file}"
    fi

    o2=$(cat ${rep}/${channel}_masses_${rep}_Sp${N}.dat | head -1 | awk '{print $1}')
    oError2=$(cat ${rep}/${channel}_masses_${rep}_Sp${N}.dat | head -1 | awk '{print $2}')

    if [ ${o2} == "--" ]
    then
        echo "-- & -- \\\\" >> "${output_file}"
    else
        o=$(echo ${o2} | awk '{print sqrt($1)}')
        oError=$(echo ${oError2} ${o} | awk '{print $1 /2 /$2}')
        num=1
        tmp=${oError}
        while (( $(echo "$tmp < 1.0" | bc -l) ))
        do
            tmp=$(echo $tmp | awk '{print $1 * 10}')
            num=`expr $(($num +1 ))`
        done
        oErrRounded=$(printf "%.${num}f" ${oError})
        chi2=$(cat ${rep}/${channel}_masses_${rep}_Sp${N}.dat | tail -1 | awk '{print $1}' | xargs printf "%1.2f")

        echo -n $(printf "%.${num}f" $o) >> "${output_file}"
        echo -n "(" >> "${output_file}"
        echo -n $(printf "%s" "${oErrRounded}" | cut -c `expr $((${num}+1))`-) >> "${output_file}"
        echo -n ") & " >> "${output_file}"
        if (( $(echo "${chi2} > 3.0" |bc -l) ));
        then
            echo -n "\\textcolor{red}{${chi2}}" >> "${output_file}"
        else
            echo -n "${chi2}" >> "${output_file}"
        fi
        echo "\\\\" >> "${output_file}"
    fi

    echo "\\cline{2-4}" >> "${output_file}"

    channel="axialtensor"

    if [ "${rep}" == "F" ]
    then
        echo -n " & \$\\hat{m}_{\\rm AT}\$ & " >> "${output_file}"
    elif [ "${rep}" == "AS" ]
    then
        echo -n " & \$\\hat{m}_{\\rm at}\$ & " >> "${output_file}"
    else
        echo -n " & \$\\hat{m}_{\\mathcal{AT}}\$ & " >> "${output_file}"
    fi

    o2=$(cat ${rep}/${channel}_masses_${rep}_Sp${N}.dat | head -1 | awk '{print $1}')
    oError2=$(cat ${rep}/${channel}_masses_${rep}_Sp${N}.dat | head -1 | awk '{print $2}')

    if [ ${o2} == "--" ]
    then
        echo "-- & -- \\\\" >> "${output_file}"
    else
        o=$(echo ${o2} | awk '{print sqrt($1)}')
        oError=$(echo ${oError2} ${o} | awk '{print $1 /2 /$2}')
        num=1
        tmp=${oError}
        while (( $(echo "$tmp < 1.0" | bc -l) ))
        do
            tmp=$(echo $tmp | awk '{print $1 * 10}')
            num=`expr $(($num +1 ))`
        done
        oErrRounded=$(printf "%.${num}f" ${oError})
        chi2=$(cat ${rep}/${channel}_masses_${rep}_Sp${N}.dat | tail -1 | awk '{print $1}' | xargs printf "%1.2f")

        echo -n $(printf "%.${num}f" $o) >> "${output_file}"
        echo -n "(" >> "${output_file}"
        echo -n $(printf "%s" "${oErrRounded}" | cut -c `expr $((${num}+1))`-) >> "${output_file}"
        echo -n ") & " >> "${output_file}"
        if (( $(echo "${chi2} > 3.0" |bc -l) ));
        then
            echo -n "\\textcolor{red}{${chi2}}" >> "${output_file}"
        else
            echo -n ${chi2} >> "${output_file}"
        fi

        echo "\\\\" >> "${output_file}"
    fi

    echo "\\hline" >> "${output_file}"
    echo "\\hline" >> "${output_file}"
done

sed -i '' -e '$ d' "${output_file}"

echo "\\end{tabular}" >> "${output_file}"
