#!/bin/bash

basedir="$1"
output_file="${PWD}/$2"
caption_file="${PWD}/$3"

cd "${basedir}"

caption_chisquare=""

Reps="F AS S"
[ -f "${output_file}" ] && rm "${output_file}"

# Write headers for the LaTeX table:
echo "\\begin{tabular}{|c|c|c|c|}" >> "${output_file}"
echo "\\hline" >> "${output_file}"
echo "\\multicolumn{4}{|c|}{\$Sp(\\infty)\$}\\\\" >> "${output_file}"
echo "\\hline" >> "${output_file}"
echo "Representation & Channel & Chiral limit & \$\\chi^2/{\\rm d.o.f.}\$\\\\" >> "${output_file}"
echo "\\hline" >> "${output_file}"

for rep in ${Reps}
do
    channel="PS"

    if [ "${rep}" == "F" ]
    then
        echo -n "\\multirow{8}{*}{Fundamental} & \$\\hat{f}^2_{\\rm PS}/N_c\$ & " >> "${output_file}"
    elif [ "${rep}" == "AS" ]
    then
        echo -n "\\multirow{8}{*}{Antisymmetric} & \$\\hat{f}^2_{\\rm ps}/N_c^2\$ & " >> "${output_file}"
    else
        echo -n "\\multirow{8}{*}{Symmetric} & \$\\hat{f}^2_{\\mathcal{PS}}/N_c^2\$ & " >> "${output_file}"
    fi

    o=$(cat ${rep}_decayconsts.txt | grep -w ${channel} | awk '{print $2}')
    oError=$(cat ${rep}_decayconsts.txt | grep -w ${channel} | awk '{print $3}')

    if [ ${o} == "--" ]
    then
        echo "\$\cdots\$ & \$\cdots\$ \\\\" >> "${output_file}"
    else
        num=1
        tmp=${oError}
        while (( $(echo "$tmp < 1.0" | bc -l) ))
        do
            tmp=$(echo $tmp | awk '{print $1 * 10}')
            num=`expr $(($num +1 ))`
        done
        oErrRounded=$(printf "%.${num}f" ${oError})
        chi2=$(cat ${rep}_decayconsts.txt | grep -w ${channel} | awk '{print $4}')
        if [[ "${chi2}" != "--" ]]
        then
            chi2=$(echo ${chi2} | xargs printf "%1.2f")
        fi

        echo -n $(printf "%.${num}f" $o) >> "${output_file}"
        echo -n "(" >> "${output_file}"
        echo -n $(printf "%s" "${oErrRounded}" | cut -c `expr $((${num}+1))`-) >> "${output_file}"
        echo -n ") & " >> "${output_file}"
        if [[ "${chi2}" != "--" ]] && (( $(echo "${chi2} > 3.0" |bc -l) ));
        then
            echo -n "\\textcolor{red}{${chi2}}" >> "${output_file}"
            if [[ "" == "${caption_chisquare}" ]]
            then
                caption_chisquare=${chi2}
            fi
        elif [[ "${chi2}" != "--" ]]
        then
            echo -n "${chi2}" >> "${output_file}"
        else
            echo -n "\$\cdots\$" >> "${output_file}"
        fi
        echo "\\\\" >> "${output_file}"
    fi

    echo "\\cline{2-4}" >> "${output_file}"

    channel="V"

    if [ "${rep}" == "F" ]
    then
        echo -n " & \$\\hat{f}^2_{\\rm V}/N_c\$ & " >> "${output_file}"
    elif [ "${rep}" == "AS" ]
    then
        echo -n " & \$\\hat{f}^2_{\\rm v}/N_c^2\$ & " >> "${output_file}"
    else
        echo -n " & \$\\hat{f}^2_{\\mathcal{V}}/N_c^2\$ & " >> "${output_file}"
    fi

    o=$(cat ${rep}_decayconsts.txt | grep -w ${channel} | awk '{print $2}')
    oError=$(cat ${rep}_decayconsts.txt | grep -w ${channel} | awk '{print $3}')

    if [ ${o} == "--" ]
    then
        echo "\$\cdots\$ & \$\cdots\$ \\\\" >> "${output_file}"
    else
        num=1
        tmp=${oError}
        while (( $(echo "$tmp < 1.0" | bc -l) ))
        do
            tmp=$(echo $tmp | awk '{print $1 * 10}')
            num=`expr $(($num +1 ))`
        done
        oErrRounded=$(printf "%.${num}f" ${oError})
        chi2=$(cat ${rep}_decayconsts.txt | grep -w ${channel} | awk '{print $4}')
        if [[ "${chi2}" != "--" ]]
        then
            chi2=$(echo ${chi2} | xargs printf "%1.2f")
        fi

        echo -n $(printf "%.${num}f" $o) >> "${output_file}"
        echo -n "(" >> "${output_file}"
        echo -n $(printf "%s" "${oErrRounded}" | cut -c `expr $((${num}+1))`-) >> "${output_file}"
        echo -n ") & " >> "${output_file}"
        if [[ "${chi2}" != "--" ]] && (( $(echo "${chi2} > 3.0" |bc -l) ));
        then
            echo -n "\\textcolor{red}{${chi2}}" >> "${output_file}"
            if [[ "" == "${caption_chisquare}" ]]
            then
                caption_chisquare=${chi2}
            fi
        elif [[ "${chi2}" != "--" ]]
        then
            echo -n "${chi2}" >> "${output_file}"
        else
            echo -n "\$\cdots\$" >> "${output_file}"
        fi
        echo "\\\\" >> "${output_file}"
    fi

    echo "\\cline{2-4}" >> "${output_file}"

    channel="AV"

    if [ "${rep}" == "F" ]
    then
        echo -n " & \$\\hat{f}^2_{\\rm AV}/N_c\$ & " >> "${output_file}"
    elif [ "${rep}" == "AS" ]
    then
        echo -n " & \$\\hat{f}^2_{\\rm av}/N_c^2\$ & " >> "${output_file}"
    else
        echo -n " & \$\\hat{f}^2_{\\mathcal{AV}}/N_c^2\$ & " >> "${output_file}"
    fi

    o=$(cat ${rep}_decayconsts.txt | grep -w ${channel} | awk '{print $2}')
    oError=$(cat ${rep}_decayconsts.txt | grep -w ${channel} | awk '{print $3}')

    if [ ${o} == "--" ]
    then
        echo "\$\cdots\$ & \$\cdots\$ \\\\" >> "${output_file}"
    else
        num=1
        tmp=${oError}
        while (( $(echo "$tmp < 1.0" | bc -l) ))
        do
            tmp=$(echo $tmp | awk '{print $1 * 10}')
            num=`expr $(($num +1 ))`
        done
        oErrRounded=$(printf "%.${num}f" ${oError})
        chi2=$(cat ${rep}_decayconsts.txt | grep -w ${channel} | awk '{print $4}')
        if [[ "${chi2}" != "--" ]]
        then
            chi2=$(echo ${chi2} | xargs printf "%1.2f")
        fi

        echo -n $(printf "%.${num}f" $o) >> "${output_file}"
        echo -n "(" >> "${output_file}"
        echo -n $(printf "%s" "${oErrRounded}" | cut -c `expr $((${num}+1))`-) >> "${output_file}"
        echo -n ") & " >> "${output_file}"
        if [[ "${chi2}" != "--" ]] && (( $(echo "${chi2} > 3.0" |bc -l) ));
        then
            echo -n "\\textcolor{red}{${chi2}}" >> "${output_file}"
            if [[ "" == "${caption_chisquare}" ]]
            then
                caption_chisquare=${chi2}
            fi
        elif [[ "${chi2}" != "--" ]]
        then
            echo -n "${chi2}" >> "${output_file}"
        else
            echo -n "\$\cdots\$" >> "${output_file}"
        fi
        echo "\\\\" >> "${output_file}"
    fi

    echo "\\cline{2-4}" >> "${output_file}"

    channel="V"

    if [ "${rep}" == "F" ]
    then
        echo -n " & \$\\hat{m}^2_{\\rm V}\$ & " >> "${output_file}"
    elif [ "${rep}" == "AS" ]
    then
        echo -n " & \$\\hat{m}^2_{\\rm v}\$ & " >> "${output_file}"
    else
        echo -n " & \$\\hat{m}^2_{\\mathcal{V}}\$ & " >> "${output_file}"
    fi

    o=$(cat ${rep}_masses.txt | grep -w ${channel} | awk '{print $2}')
    oError=$(cat ${rep}_masses.txt | grep -w ${channel} | awk '{print $3}')

    if [ ${o} == "--" ]
    then
        echo "\$\cdots\$ & \$\cdots\$ \\\\" >> "${output_file}"
    else
        num=1
        tmp=${oError}
        while (( $(echo "$tmp < 1.0" | bc -l) ))
        do
            tmp=$(echo $tmp | awk '{print $1 * 10}')
            num=`expr $(($num +1 ))`
        done
        oErrRounded=$(printf "%.${num}f" ${oError})
        chi2=$(cat ${rep}_masses.txt | grep -w ${channel} | awk '{print $4}')
        if [[ "${chi2}" != "--" ]]
        then
            chi2=$(echo ${chi2} | xargs printf "%1.2f")
        fi

        echo -n $(printf "%.${num}f" $o) >> "${output_file}"
        echo -n "(" >> "${output_file}"
        echo -n $(printf "%s" "${oErrRounded}" | cut -c `expr $((${num}+1))`-) >> "${output_file}"
        echo -n ") & " >> "${output_file}"
        if [[ "${chi2}" != "--" ]] && (( $(echo "${chi2} > 3.0" |bc -l) ));
        then
            echo -n "\\textcolor{red}{${chi2}}" >> "${output_file}"
            if [[ "" == "${caption_chisquare}" ]]
            then
                caption_chisquare=${chi2}
            fi
        elif [[ "${chi2}" != "--" ]]
        then
            echo -n "${chi2}" >> "${output_file}"
        else
            echo -n "\$\cdots\$" >> "${output_file}"
        fi
        echo "\\\\" >> "${output_file}"
    fi

    echo "\\cline{2-4}" >> "${output_file}"

    channel="AV"

    if [ "${rep}" == "F" ]
    then
        echo -n " & \$\\hat{m}^2_{\\rm AV}\$ & " >> "${output_file}"
    elif [ "${rep}" == "AS" ]
    then
        echo -n " & \$\\hat{m}^2_{\\rm av}\$ & " >> "${output_file}"
    else
        echo -n " & \$\\hat{m}^2_{\\mathcal{AV}}\$ & " >> "${output_file}"
    fi

    o=$(cat ${rep}_masses.txt | grep -w ${channel} | awk '{print $2}')
    oError=$(cat ${rep}_masses.txt | grep -w ${channel} | awk '{print $3}')

    if [ ${o} == "--" ]
    then
        echo "\$\cdots\$ & \$\cdots\$ \\\\" >> "${output_file}"
    else
        num=1
        tmp=${oError}
        while (( $(echo "$tmp < 1.0" | bc -l) ))
        do
            tmp=$(echo $tmp | awk '{print $1 * 10}')
            num=`expr $(($num +1 ))`
        done
        oErrRounded=$(printf "%.${num}f" ${oError})
        chi2=$(cat ${rep}_masses.txt | grep -w ${channel} | awk '{print $4}')
        if [[ "${chi2}" != "--" ]]
        then
            chi2=$(echo ${chi2} | xargs printf "%1.2f")
        fi

        echo -n $(printf "%.${num}f" $o) >> "${output_file}"
        echo -n "(" >> "${output_file}"
        echo -n $(printf "%s" "${oErrRounded}" | cut -c `expr $((${num}+1))`-) >> "${output_file}"
        echo -n ") & " >> "${output_file}"
        if [[ "${chi2}" != "--" ]] && (( $(echo "${chi2} > 3.0" |bc -l) ));
        then
            echo -n "\\textcolor{red}{${chi2}}" >> "${output_file}"
            if [[ "" == "${caption_chisquare}" ]]
            then
                caption_chisquare=${chi2}
            fi
        elif [[ "${chi2}" != "--" ]]
        then
            echo -n "${chi2}" >> "${output_file}"
        else
            echo -n "\$\cdots\$" >> "${output_file}"
        fi
        echo "\\\\" >> "${output_file}"
    fi

    echo "\\cline{2-4}" >> "${output_file}"

    channel="S"

    if [ "${rep}" == "F" ]
    then
        echo -n " & \$\\hat{m}^2_{\\rm S}\$ & " >> "${output_file}"
    elif [ "${rep}" == "AS" ]
    then
        echo -n " & \$\\hat{m}^2_{\\rm s}\$ & " >> "${output_file}"
    else
        echo -n " & \$\\hat{m}^2_{\\mathcal{S}}\$ & " >> "${output_file}"
    fi

    o=$(cat ${rep}_masses.txt | grep -w ${channel} | awk '{print $2}')
    oError=$(cat ${rep}_masses.txt | grep -w ${channel} | awk '{print $3}')

    if [ ${o} == "--" ]
    then
        echo "\$\cdots\$ & \$\cdots\$ \\\\" >> "${output_file}"
    else
        num=1
        tmp=${oError}
        while (( $(echo "$tmp < 1.0" | bc -l) ))
        do
            tmp=$(echo $tmp | awk '{print $1 * 10}')
            num=`expr $(($num +1 ))`
        done
        oErrRounded=$(printf "%.${num}f" ${oError})
        chi2=$(cat ${rep}_masses.txt | grep -w ${channel} | awk '{print $4}')
        if [[ "${chi2}" != "--" ]]
        then
            chi2=$(echo ${chi2} | xargs printf "%1.2f")
        fi

        echo -n $(printf "%.${num}f" $o) >> "${output_file}"
        echo -n "(" >> "${output_file}"
        echo -n $(printf "%s" "${oErrRounded}" | cut -c `expr $((${num}+1))`-) >> "${output_file}"
        echo -n ") & " >> "${output_file}"
        if [[ "${chi2}" != "--" ]] && (( $(echo "${chi2} > 3.0" |bc -l) ));
        then
            echo -n "\\textcolor{red}{${chi2}}" >> "${output_file}"
            if [[ "" == "${caption_chisquare}" ]]
            then
                caption_chisquare=${chi2}
            fi
        elif [[ "${chi2}" != "--" ]]
        then
            echo -n "${chi2}" >> "${output_file}"
        else
            echo -n "\$\cdots\$" >> "${output_file}"
        fi
        echo "\\\\" >> "${output_file}"
    fi

    echo "\\cline{2-4}" >> "${output_file}"

    channel="T"

    if [ "${rep}" == "F" ]
    then
        echo -n " & \$\\hat{m}^2_{\\rm T}\$ & " >> "${output_file}"
    elif [ "${rep}" == "AS" ]
    then
        echo -n " & \$\\hat{m}^2_{\\rm t}\$ & " >> "${output_file}"
    else
        echo -n " & \$\\hat{m}^2_{\\mathcal{T}}\$ & " >> "${output_file}"
    fi

    o=$(cat ${rep}_masses.txt | grep -w ${channel} | awk '{print $2}')
    oError=$(cat ${rep}_masses.txt | grep -w ${channel} | awk '{print $3}')

    if [ ${o} == "--" ]
    then
        echo "\$\cdots\$ & \$\cdots\$ \\\\" >> "${output_file}"
    else
        num=1
        tmp=${oError}
        while (( $(echo "$tmp < 1.0" | bc -l) ))
        do
            tmp=$(echo $tmp | awk '{print $1 * 10}')
            num=`expr $(($num +1 ))`
        done
        oErrRounded=$(printf "%.${num}f" ${oError})
        chi2=$(cat ${rep}_masses.txt | grep -w ${channel} | awk '{print $4}')
        if [[ "${chi2}" != "--" ]]
        then
            chi2=$(echo ${chi2} | xargs printf "%1.2f")
        fi

        echo -n $(printf "%.${num}f" $o) >> "${output_file}"
        echo -n "(" >> "${output_file}"
        echo -n $(printf "%s" "${oErrRounded}" | cut -c `expr $((${num}+1))`-) >> "${output_file}"
        echo -n ") & " >> "${output_file}"
        if [[ "${chi2}" != "--" ]] && (( $(echo "${chi2} > 3.0" |bc -l) ));
        then
            echo -n "\\textcolor{red}{${chi2}}" >> "${output_file}"
            if [[ "" == "${caption_chisquare}" ]]
            then
                caption_chisquare=${chi2}
            fi
        elif [[ "${chi2}" != "--" ]]
        then
            echo -n "${chi2}" >> "${output_file}"
        else
            echo -n "\$\cdots\$" >> "${output_file}"
        fi
        echo "\\\\" >> "${output_file}"
    fi

    echo "\\cline{2-4}" >> "${output_file}"

    channel="AT"

    if [ "${rep}" == "F" ]
    then
        echo -n " & \$\\hat{m}^2_{\\rm AT}\$ & " >> "${output_file}"
    elif [ "${rep}" == "AS" ]
    then
        echo -n " & \$\\hat{m}^2_{\\rm at}\$ & " >> "${output_file}"
    else
        echo -n " & \$\\hat{m}^2_{\\mathcal{AT}}\$ & " >> "${output_file}"
    fi

    o=$(cat ${rep}_masses.txt | grep -w ${channel} | awk '{print $2}')
    oError=$(cat ${rep}_masses.txt | grep -w ${channel} | awk '{print $3}')

    if [ ${o} == "--" ]
    then
        echo "\$\cdots\$ & \$\cdots\$ \\\\" >> "${output_file}"
    else
        num=1
        tmp=${oError}
        while (( $(echo "$tmp < 1.0" | bc -l) ))
        do
            tmp=$(echo $tmp | awk '{print $1 * 10}')
            num=`expr $(($num +1 ))`
        done
        oErrRounded=$(printf "%.${num}f" ${oError})
        chi2=$(cat ${rep}_masses.txt | grep -w ${channel} | awk '{print $4}')
        if [[ "${chi2}" != "--" ]]
        then
            chi2=$(echo ${chi2} | xargs printf "%1.2f")
        fi

        echo -n $(printf "%.${num}f" $o) >> "${output_file}"
        echo -n "(" >> "${output_file}"
        echo -n $(printf "%s" "${oErrRounded}" | cut -c `expr $((${num}+1))`-) >> "${output_file}"
        echo -n ") & " >> "${output_file}"
        if [[ "${chi2}" != "--" ]] && (( $(echo "${chi2} > 3.0" |bc -l) ));
        then
            echo -n "\\textcolor{red}{${chi2}}" >> "${output_file}"
            if [[ "" == "${caption_chisquare}" ]]
            then
                caption_chisquare=${chi2}
            fi
        elif [[ "${chi2}" != "--" ]]
        then
            echo -n "${chi2}" >> "${output_file}"
        else
            echo -n "\$\cdots\$" >> "${output_file}"
        fi
        echo "\\\\" >> "${output_file}"
    fi

    echo "\\hline" >> "${output_file}"
    echo "\\hline" >> "${output_file}"
done

sed -i '' -e '$ d' "${output_file}"

echo "\\end{tabular}" >> "${output_file}"
echo "captionchisquarelargeN ${caption_chisquare}" >> "${caption_file}"
