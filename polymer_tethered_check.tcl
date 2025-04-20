############################################################# #######################
#  A SAW linear polymer fixed at both ends. Proteins form droplets on top of it  ####
#             IMSc September 2022                       ######################## ####
#############################################################
  
source ./functions.tcl
puts " "
puts "======================================================="
puts "=               ESPResSo Simulation                   ="
puts "=                                                     ="
puts "======================================================="
puts " "
puts "[code_info]"

#############################################################
# 2 Preparing the System                                    #
#                                                           #
############################################################# 
set n_mono 500
set n_poly 1
set n_binders 3840 
set box_length_x 80.0
set box_length_y 80.0
set box_length_z 600.0
set n_part [expr $n_mono * $n_poly+ $n_binders]
setmd periodic 0 0 0
#############################################################
# 3 Initializing ESPResSo                                   #
#                                                           #
#############################################################
setmd box_l $box_length_x $box_length_y $box_length_z
setmd time_step 0.01
setmd skin 0.4
set kBT 1.0
set gamma 0.1
thermostat langevin $kBT $gamma

puts [setmd box_l]
puts [setmd time_step]
puts [setmd skin]
puts [integrate]
puts [thermostat]


############################################################################
##############################################################
set in [open "|gzip -cd - < ./check_point.gz" "r"] #loading the saved check_point file
while { [blockfile $in read auto] != "eof" } {}
close $in

#tethering the ends of the polymer
part 0 fix
part [expr $n_mono-1] fix
##############################################################
set warm_loop 1
set warm_step 0
set del_r [expr $rcap/$warm_loop]
for { set k 0 } { $k < $warm_loop} { incr k } {
    inter forcecap individual
    integrate $warm_step   
    
    set box_length_z [expr $box_length_z*1]
    
    change_volume $box_length_z z
    constraint delete
    
    constraint wall normal 1 0 0 dist 0 type 25
    constraint wall normal -1 0 0 dist -$box_length_x type 25
    constraint wall normal 0 1 0 dist 0 type 25
    constraint wall normal 0 -1 0 dist -$box_length_y type 25
    constraint wall normal 0 0 1 dist 0 type 25
    constraint wall normal 0 0 -1 dist -$box_length_z type 25
    
    set rcap [expr $rcap-$del_r]
    
    #mono-mono interaction
    for {set i 1} {$i < $type_max } {incr i} {
     for {set j $i} {$j < $type_max } {incr j} {
     inter $i $j lennard-jones 1.0 1.0 $rcut 0.25 0.0 $rcap 0.0
     }
    }

    inter 24 24 lennard-jones $epsilon_cc 1.0 $rcut_cc 0.25 0.0 $rcap 0.0

    #mono-protein interaction
    for {set i 1} {$i < $type_max} {incr i} {
     set j [split [lindex $epsilon_values $i] " "]
     inter $i 24 lennard-jones $j 1.0 $rcut_bind 0.25 0.0 $rcap 0.0
    }
    
    puts "Warm up integration $k"    
}
puts "AFTER WARM UP: [constraint]"
set rcap 0
inter forcecap 0



#mono-mono interaction
for {set i 1} {$i < $type_max } {incr i} {
 for {set j $i} {$j < $type_max } {incr j} {
 inter $i $j lennard-jones 1.0 1.0 $rcut 0.25 0.0 $rcap 0.0

 }
}


for {set i 1} {$i < $type_max} {incr i} {
 set j [split [lindex $epsilon_values $i] " "]
 
 inter $i 24 lennard-jones $j 1.0 $rcut_bind 0.25 0.0 $rcap 0.0
 
}
inter 24 24 lennard-jones $epsilon_cc 1.0 $rcut_cc 0.25 0.0 $rcap 0.0


###########################################################################################################
###########################################################################################################
set obs_en [open "./energy.dat" "w"]
set vsf [open "./config.vsf" "w"]
set vtf [open "./config.vtf" "w"]
###########################################################################################################
set f [open "./config/config_12001" "w"]
blockfile $f write particles {id pos type}
close $f

set tstart [setmd time]

writevsf $vsf
writevcf $vtf


############################################################
# 8 Integration                                            #
#                                                          #
############################################################


set n_steps 5000
set n_equi 3000

for {set i 0} { $i < $n_equi } { incr i} {
    puts $i
    integrate $n_steps
puts $obs_en "[expr $i] [analyze energy total] [analyze energy kinetic]"
writevcf $vtf
set j [expr $i+2+12000]
set f [open "./config/config_$j" "w"]
blockfile $f write particles {id pos type}
close $f

set check_point [open "|gzip -c - > ./check_point/check_point.gz" "w"]
save_sim $check_point "id pos v f type" "all"
close $check_point
}

close $obs_en
close $vsf
close $vtf
############################################################
