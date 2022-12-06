#!/bin/bash

echo "Enter dirctory for deleting trajectories in:"
read dir

function select_option {

    # little helpers for terminal print control and key input
    ESC=$( printf "\033")
    cursor_blink_on()  { printf "$ESC[?25h"; }
    cursor_blink_off() { printf "$ESC[?25l"; }
    cursor_to()        { printf "$ESC[$1;${2:-1}H"; }
    print_option()     { printf "   $1 "; }
    print_selected()   { printf "  $ESC[7m $1 $ESC[27m"; }
    get_cursor_row()   { IFS=';' read -sdR -p $'\E[6n' ROW COL; echo ${ROW#*[}; }
    key_input()        { read -s -n3 key 2>/dev/null >&2
                         if [[ $key = $ESC[A ]]; then echo up;    fi
                         if [[ $key = $ESC[B ]]; then echo down;  fi
                         if [[ $key = ""     ]]; then echo enter; fi; }

    # initially print empty new lines (scroll down if at bottom of screen)
    for opt; do printf "\n"; done

    # determine current screen position for overwriting the options
    local lastrow=`get_cursor_row`
    local startrow=$(($lastrow - $#))

    # ensure cursor and input echoing back on upon a ctrl+c during read -s
    trap "cursor_blink_on; stty echo; printf '\n'; exit" 2
    cursor_blink_off

    local selected=0
    while true; do
        # print options by overwriting the last lines
        local idx=0
        for opt; do
            cursor_to $(($startrow + $idx))
            if [ $idx -eq $selected ]; then
                print_selected "$opt"
            else
                print_option "$opt"
            fi
            ((idx++))
        done

        # user key control
        case `key_input` in
            enter) break;;
            up)    ((selected--));
                   if [ $selected -lt 0 ]; then selected=$(($# - 1)); fi;;
            down)  ((selected++));
                   if [ $selected -ge $# ]; then selected=0; fi;;
        esac
    done

    # cursor position back to normal
    cursor_to $lastrow
    printf "\n"
    cursor_blink_on

    return $selected
}


# pick the engine
echo "Select one option using up/down keys and enter to confirm:"
echo

options=("AMBER" "GROMACS" "SOMD" "ALL")

select_option "${options[@]}"
choice=$?

# pick if all trajectories
echo "Select one option using up/down keys and enter to confirm:"
echo

options2=("ALL" "KEEP LAMBDA 0 AND 1" "KEEP LAMBDA 0, 0.5 AND 1")

select_option "${options2[@]}"
choice2=$?

echo "Deleting trajectories for ${options[$choice]} and ${options2[$choice2]} trajectories..."

#removing all traj
# remove amber trajectories
if [ "${options[$choice]}" = "AMBER" ] || [ "${options[$choice]}" = "ALL" ] && [ "${options2[$choice2]}" = "ALL" ]; then
for i in $(find $dir -name 'amber.nc');
do
    rm -rf $i
done
fi

# remove gromacs trajectories
if [ "${options[$choice]}" = "GROMACS" ] || [ "${options[$choice]}" = "ALL" ] && [ "${options2[$choice2]}" = "ALL" ]; then
for i in $(find $dir -name 'gromacs.trr');
do
    rm -rf $i
done
fi

# remove somd trajectories
if [ "${options[$choice]}" = "SOMD" ] || [ "${options[$choice]}" = "ALL" ] && [ "${options2[$choice2]}" = "ALL" ]; then
for i in $(find $dir -name 'traj*.dcd');
do
    rm -rf $i
done
for i in $(find $dir -name '*.s3*');
do
    rm -rf $i
done
fi

#removing all except lambda 0 and 1
# remove amber trajectories
if [ "${options[$choice]}" = "AMBER" ] || [ "${options[$choice]}" = "ALL" ] && [ "${options2[$choice2]}" = "KEEP LAMBDA 0 AND 1" ]; then
for i in $(find $dir -name 'amber.nc');

do
if [[ "$i" != *"lambda_0.0000"* ]] && [[ "$i" != *"lambda_1.0000"* ]]; then
    rm -rf $i
fi
done

fi

# remove gromacs trajectories
if [ "${options[$choice]}" = "GROMACS" ] || [ "${options[$choice]}" = "ALL" ]&& [ "${options2[$choice2]}" = "KEEP LAMBDA 0 AND 1" ]; then
for i in $(find $dir -name 'gromacs.trr');

do
if [[ "$i" != *"lambda_0.0000"* ]] && [[ "$i" != *"lambda_1.0000"* ]]; then
    rm -rf $i
fi
done

fi

# remove somd trajectories
if [ "${options[$choice]}" = "SOMD" ] || [ "${options[$choice]}" = "ALL" ]&& [ "${options2[$choice2]}" = "KEEP LAMBDA 0 AND 1" ]; then
for i in $(find $dir -name 'traj*.dcd');

do
if [[ "$i" != *"lambda_0.0000"* ]] && [[ "$i" != *"lambda_1.0000"* ]]; then
    rm -rf $i
fi
done

for i in $(find $dir -name '*.s3*');

do
if [[ "$i" != *"lambda_0.0000"* ]] && [[ "$i" != *"lambda_1.0000"* ]]; then
    rm -rf $i
fi
done

fi

#removing all except lambda 0, 0.5 and 1
# remove amber trajectories
if [ "${options[$choice]}" = "AMBER" ] || [ "${options[$choice]}" = "ALL" ] && [ "${options2[$choice2]}" = "KEEP LAMBDA 0, 0.5 AND 1" ]; then
for i in $(find $dir -name 'amber.nc');

do
if [[ "$i" != *"lambda_0.0000"* ]] && [[ "$i" != *"lambda_0.5000"* ]] && [[ "$i" != *"lambda_1.0000"* ]]; then
    rm -rf $i
fi
done

fi

# remove gromacs trajectories
if [ "${options[$choice]}" = "GROMACS" ] || [ "${options[$choice]}" = "ALL" ]&& [ "${options2[$choice2]}" = "KEEP LAMBDA 0, 0.5 AND 1" ]; then
for i in $(find $dir -name 'gromacs.trr');

do
if [[ "$i" != *"lambda_0.0000"* ]] && [[ "$i" != *"lambda_0.5000"* ]] && [[ "$i" != *"lambda_1.0000"* ]]; then
    rm -rf $i
fi
done

fi

# remove somd trajectories
if [ "${options[$choice]}" = "SOMD" ] || [ "${options[$choice]}" = "ALL" ]&& [ "${options2[$choice2]}" = "KEEP LAMBDA 0, 0.5 AND 1" ]; then
for i in $(find $dir -name 'traj*.dcd');

do
if [[ "$i" != *"lambda_0.0000"* ]] && [[ "$i" != *"lambda_0.5000"* ]] && [[ "$i" != *"lambda_1.0000"* ]]; then
    rm -rf $i
fi
done

for i in $(find $dir -name '*.s3*');

do
if [[ "$i" != *"lambda_0.0000"* ]] && [[ "$i" != *"lambda_0.5000"* ]] && [[ "$i" != *"lambda_1.0000"* ]]; then
    rm -rf $i
fi
done

fi

echo "Done deleting."

# TODO fix so can pick tranjectories at only certain lambda windows ?