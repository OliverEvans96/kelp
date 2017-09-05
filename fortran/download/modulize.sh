#!/bin/bash
fname=$1
func_list=( $(awk '$1 == "function" {print $2}' < "$fname") )
cat $fname | awk -v funcs="${func_list[*]}" '
BEGIN {
    split(funcs, func_arr); 
    #https://stackoverflow.com/questions/26746361/check-if-array-contains-value
    for(i in func_arr) {
        values[func_arr[i]] = "";
    }
}

$1~/(function|subroutine)/ {current=$2}
$1~/(logical|integer|real)/ && ($NF in values) && ($NF != current) {printf "!123"}
{print}
END{}
'

