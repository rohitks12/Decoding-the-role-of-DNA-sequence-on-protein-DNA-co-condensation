#############################################################
#                                                           #
#  Functions                                               #
#                                                           # 
#############################################################
proc save_sim {cfile parinfo range } {
# write all available sim information to channel cfile
# in block format
 blockfile $cfile write variable all
 blockfile $cfile write tclvariable all
 blockfile $cfile write particles $parinfo $range
 blockfile $cfile write interactions
 blockfile $cfile write bonds $range
 blockfile $cfile write random
 blockfile $cfile write seed
 blockfile $cfile write bitrandom
 blockfile $cfile write bitseed
}

proc save_sim1 {cfile parinfo range } {
# write all available sim information to channel cfile
# in block format
blockfile $cfile write variable all
blockfile $cfile write tclvariable all
blockfile $cfile write particles $parinfo $range
}

proc readData {filename} {
    set result {}
    set f [open $filename r]
    foreach line [split [read $f] \n] {
        lappend result $line
    }
    return $result
}
