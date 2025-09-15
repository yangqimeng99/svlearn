#!/bin/bash

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 [INPUT_FILE]"
    echo ""
    echo "Replace numeric repeat_class labels in the SVLearn sv_feature.tsv"
    exit 0
fi

awk 'BEGIN{
    
    map["0"]="Low_Repeat";
    map["1"]="LINE";
    map["2"]="SINE";
    map["3"]="LTR";
    map["4"]="DNA";
    map["5"]="Mixed_TEs";
    map["6"]="Satellite";
    map["7"]="VNTR";
    map["8"]="STR";
    map["9"]="Mixed_Repeat";
    FS=OFS="\t"
}
NR==1{
    print $0;
    next;
}
{
    if ($4 in map) {
        $4 = map[$4];
    }
    print $0;
}' "${1:-/dev/stdin}"
