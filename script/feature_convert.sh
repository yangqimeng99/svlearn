#!/bin/bash

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 [INPUT_FILE]"
    echo ""
    echo "Replace numeric feature labels in the SVLearn feature importance"
    exit 0
fi

awk 'BEGIN{
    
    map["sv_type"]="Type";
    map["repeat_class"]="Repeat_class";
    map["SDs_class"]="SD_class";
    map["ref_GT"]="Ref_GT";
    map["alt_GT"]="Alt_GT";
    map["ref_FT"]="Ref_FT";
    map["alt_FT"]="Alt_FT";
    map["remainder__sv_length"]="Length";
    map["remainder__repeat_content"]="Repeat_content";
    map["remainder__SDs_content"]="SD_content";
    map["remainder__TR_length"]="TR_length";
    map["remainder__TR_content"]="TR_content";
    map["remainder__mappab_median"]="Mappability";
    map["remainder__sv_gc_content"]="GC_content";
    map["remainder__front_depth"]="RD_upstream";
    map["remainder__end_depth"]="RD_downstream";
    map["remainder__RD_lp_homref"]="RD_HomRef";
    map["remainder__RD_lp_het"]="RD_Het";
    map["remainder__RD_lp_homalt"]="RD_HomAlt";
    map["remainder__BP_lp_homref"]="BP_HomRef";
    map["remainder__BP_lp_het"]="BP_Het";
    map["remainder__BP_lp_homalt"]="BP_HomAlt";
    map["remainder__DP"]="DP";
    map["remainder__PL"]="PL";
    FS=OFS="\t"
}
NR==1{
    print $0;
    next;
}
{
    if ($2 in map) {
        $2 = map[$2];
    }
    print $0;
}' "${1:-/dev/stdin}"
