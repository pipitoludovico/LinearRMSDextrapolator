puts 'usage: vmd -dispdev text -e rmsd_selection.tcl -args "alignment selection" "rmsd selection"  "outputfilename.dat"'


puts "\n\nIterating through cwd"
set file_list [glob -nocomplain -directory "." *]

puts "\nChecking topology files..."
set sys_list [glob -nocomplain -directory "system" *]

set sys_File ""
set traj ""
set top_extension ""

# look for the xtc trajectory file
foreach file $file_list {
    if {[file isfile $file] && [string equal -nocase [file extension $file] ".xtc"]} {
        set traj $file
        break
    }
}
# look for the structural file
foreach file $sys_list {
    if {[file isfile $file]} {
        set extension [file extension $file]
        if {[string equal -nocase $extension ".prmtop"]} {
            set sys_File $file
            set top_extension "parm7"
            break
        } elseif {[string equal -nocase $extension ".psf"]} {
            set sys_File $file
            set top_extension "psf"
            break
        } elseif {[string equal -nocase $extension ".gro"]} {
            set sys_File $file
            set top_extension "gro"
            break
    }
  }
}
# looking for the coordinates
foreach file $sys_list {
    if {[file isfile $file]} {
        set extension [file extension $file]
        if {[string equal -nocase $extension ".gro"]} {
            set coordinates $file
            set coord_ext "gro"
            break
        } elseif {[string equal -nocase $extension ".pdb"]} {
            set coordinates $file
            set coord_ext "pdb"
            break
        }
    }
}

puts "Extension found: $top_extension"
puts "File path: $sys_File"
puts "Traj: $traj"
puts "Coordinates file: $coordinates"
puts "\n\n"

mol load $top_extension $sys_File $coord_ext $coordinates
mol addfile $traj type xtc first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

# aligning
proc align { rmolid smolid2 seltext } {
 set ref [atomselect 0 $seltext frame 0]
 set sel [atomselect $smolid2 $seltext]
 set all [atomselect $smolid2 all]
 set n [molinfo $smolid2 get numframes]

 for { set i 1 } { $i < $n } { incr i } {
 $sel frame $i
 $all frame $i
 $all move [measure fit $sel $ref]
 }
 return
}

set selection [lindex $::argv 0]
set rmsdof [lindex $::argv 1]

set ref [atomselect 0 $rmsdof frame 0]

align 0 0 "protein and name CA"
align 0 0 $selection

# RMSD 
if {[file exists $sys_File] && [file exists $traj]} {
    set compare [atomselect top "$rmsdof"]
    set num_steps [molinfo 0 get numframes]
    set file_name [lindex $::argv 2]
    set outfile [open $file_name w]
    for {set frame 0} {$frame < $num_steps} {incr frame} {
        $compare frame $frame
        set rmsd [measure rmsd $compare $ref]
        puts $outfile "$frame\t$rmsd"
    }
    close $outfile
    quit
} else {
    puts "Error: system or trajectories not found."
}

exit
