#########################################################################################
#  A semiflexible DNA polymer fixed at both ends. Proteins form droplets with of it  ####
#                         IMSc September 2024                                        ####
#########################################################################################
  
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
set n_binders 7680
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
# 4 Interactions                                            
#Type for specific interaction between different particles and constraints
#2 : for the monomers belonging to the homogeneous DNA 
#24 : proteins
#25 : for all constraints
############################################################################
#introducing bond parameters
inter 0 harmonic 100. 1
inter 1 angle_harmonic 15 [expr [PI]]


set rcap 1.12246
set rcut_rep 1.12246
set rcut_bind 2.5
set rcut_pp 2.5

#protein-protein interaction energy
set epsilon_pp 2.0

inter forcecap individual


#polymer with different types of diffrent beads
#2 for monomers of homogeneous DNA
#24 for proteins
#25 for walls
#1 to 21 binds to 24: attractive interactions, all other interactions are repulsive

#monomer-monomer interaction
set type_max 3
#type 0 starts for TypeID sequence, 

#mono-mono interaction
for {set i 1} {$i < $type_max } {incr i} {
 for {set j $i} {$j < $type_max } {incr j} {
 inter $i $j lennard-jones 1.0 1.0 $rcut_rep 0.25 0.0 $rcap 0.0
 }
}



set fp [open "./TypeID_homo.dat" r] ;#load the input file for sequence
set index {}
set TypeID {}

while {[gets $fp line]>=0} {
    if {[llength $line]>0} {
        lappend index [lindex [split $line ","] 0]
        lappend TypeID [lindex [split $line ","] 1]
    }
}
close $fp 


#monomer protein interaction
set epsilon_values {0.0 0.0 2.00 }


for {set i 1} {$i < $type_max} {incr i} {
 set j [split [lindex $epsilon_values $i] " "]
 inter $i 24 lennard-jones $j 1.0 $rcut_bind 0.25 0.0 $rcap 0.0
 puts $i
 puts $j
}

#protein-protein interaction
inter 24 24 lennard-jones $epsilon_pp 1.0 $rcut_pp 0.25 0.0 $rcap 0.0


puts "BEFORE WARM UP: [constraint]"

set cshift 0.6
set roff 0.0
set e1 10
set e2 4
set b1 0.4
set b2 1.0
#25 for walls, constraint between wall and all monomer types
#monomer-wall interaction
for {set i 1} {$i < $type_max } {incr i} {
 inter $i 25 lj-gen [expr 2*[PI]] 1.0 1.0 $cshift $roff $e1 $e2 $b1 $b2
 }


#protein wall interactions

inter 24 25 lj-gen [expr 2*[PI]] 1.0 1.0 $cshift $roff $e1 $e2 $b1 $b2
constraint wall normal 1 0 0 dist 0 type 25
constraint wall normal -1 0 0 dist -$box_length_x type 25
constraint wall normal 0 1 0 dist 0 type 25
constraint wall normal 0 -1 0 dist -$box_length_y type 25
constraint wall normal 0 0 1 dist 0 type 25
constraint wall normal 0 0 -1 dist -$box_length_z type 25

#############################################################
# Creating Particles                                        #
#############################################################
set bond_length 1.0
set bond_angle [expr [PI]]


set fp [open "./config_Re_0.6.dat" r] ;#load the polymer configuration file
set index {}
set pos_x {}
set pos_y {}
set pos_z {}
while {[gets $fp line]>=0} {
    if {[llength $line]>0} {
        lappend index [lindex [split $line ","] 0]
        lappend pos_x [lindex [split $line ","] 1]
        lappend pos_y [lindex [split $line ","] 2]
        lappend pos_z [lindex [split $line ","] 3]
    }
}
close $fp

#semi flexible polymer
for { set i 0 } { $i < $n_mono } { incr i } {	
 set j [split [lindex $TypeID $i] " "]
 set x [split [lindex $pos_x $i] " "]
 set y [split [lindex $pos_y $i] " "]
 set z [split [lindex $pos_z $i] " "]
 part $i pos $x $y $z type $j
 puts $i  
 puts $j
 puts $x
 puts "[part $i print pos type]"
}


#settng the harmonic bonds between the monomers
for { set i 1 } { $i < $n_mono } { incr i } {
part $i bond 0 [expr $i-1]
}

#setting harmonic angle potentials between the monomers
for { set i 1 } { $i < $n_mono-1 } { incr i } {
part $i bond 1 [expr $i-1] [expr $i+1]
}

#tethering the ends of the polymer
part 0 fix
part [expr $n_mono-1] fix
 



#Introducing the proteins
set n_start [expr $n_mono*$n_poly]
for { set i $n_start } { $i < $n_part} { incr i } {
set pos_x [expr 4+0.9*rand()*$box_length_x] 
set pos_y [expr 4+0.9*rand()*$box_length_y]
set pos_z [expr 3+0.99*rand()*$box_length_z]
part $i pos $pos_x $pos_y $pos_z type 24
}

#############################################################
# Warmup                                                    #
#                                                           #
#############################################################
set warm_loop 70
set warm_step 100
set del_r [expr $rcap/$warm_loop]
for { set k 0 } { $k < $warm_loop} { incr k } {
    inter forcecap individual
    integrate $warm_step   
    
    set box_length_z [expr $box_length_z*0.99]
    
    change_volume $box_length_z z
    constraint delete

    #reducing box length to desired box_length
    
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
     inter $i $j lennard-jones 1.0 1.0 $rcut_rep 0.25 0.0 $rcap 0.0
     }
    }
  
    inter 24 24 lennard-jones $epsilon_pp 1.0 $rcut_pp 0.25 0.0 $rcap 0.0

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
 inter $i $j lennard-jones 1.0 1.0 $rcut_rep 0.25 0.0 $rcap 0.0
 }
}

#monomer-protein interaction
for {set i 1} {$i < $type_max} {incr i} {
 set j [split [lindex $epsilon_values $i] " "]
 
 inter $i 24 lennard-jones $j 1.0 $rcut_bind 0.25 0.0 $rcap 0.0
 
}
#protein-protein interaction
inter 24 24 lennard-jones $epsilon_pp 1.0 $rcut_pp 0.25 0.0 $rcap 0.0

###########################################################################################################
set obs_en [open "./energy.dat" "w"] ;#saving the total energy and kinetic energy of the system
set vsf [open "./config.vsf" "w"] ;#saving the structure file
set vtf [open "./config.vtf" "w"] ;#saving the trajectory file
###########################################################################################################
set f [open "./config/config_0" "w"] ;#saving the configuration files
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
set n_equi 10000

for {set i 0} { $i < $n_equi } { incr i} {
    puts $i
    integrate $n_steps
puts $obs_en "[expr $i] [analyze energy total] [analyze energy kinetic]"
writevcf $vtf
set j [expr $i+1]
set f [open "./config/config_$j" "w"]
blockfile $f write particles {id pos type}
close $f

#saving the check point file
set check_point [open "|gzip -c - > ./check_point/check_point.gz" "w"]
save_sim $check_point "id pos v f type" "all"
close $check_point
}

close $obs_en
close $vsf
close $vtf
############################################################
