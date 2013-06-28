# export.tcl - CGNS export routines

proc export_default {w name cmd} {
  global ProgData Font Export tcl_platform

  if {$cmd != ""} {
    set exe [lindex $cmd 0]
    set exepath [get_executable $exe]
    if {$exepath != ""} {set exe $exepath}
    if {$tcl_platform(platform) == "windows"} {
      set Export(exefile) [join [split $exe /] \\]
    } else {
      set Export(exefile) $exe
    }
    set Export(options) [lrange $cmd 1 end]
  }
  set Export(cgnsfile) $ProgData(file,name)

  toplevel $w
  wm title $w $name
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Export(done) 0}

  set labw 12

  FrameCreate $w.top -text "Export Command" -font $Font(bold)
  pack $w.top -side top -padx 5 -pady 5 -fill x
  set top [FrameGet $w.top]

  set f [frame $top.exe]
  pack $f -side top -fill x -expand 1
  label $f.lab -text Executable -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Export(exefile) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 -command "export_browse $w exefile"
  pack $f.but -side right -fill y
  $f.ent xview [$f.ent index end]

  set f [frame $top.opts]
  pack $f -side top -fill x -expand 1
  label $f.lab -text Options -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Export(options) -width 30
  pack $f.ent -side left -fill x -expand 1

  set f [frame $top.input]
  pack $f -side top -fill x -expand 1
  label $f.lab -text "CGNS Input" -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Export(cgnsfile) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 -command "export_browse $w cgnsfile"
  pack $f.but -side right -fill y

  set f [frame $top.output]
  pack $f -side top -fill x -expand 1
  label $f.lab -text "Output File" -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Export(outputfile) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 \
    -command "export_browse $w outputfile {$name}"
  pack $f.but -side right -fill y

  if {[export_buttons $w export_check]} {
    set cmd [tools_unix_path $Export(exefile)]
    foreach opt $Export(options) {
      lappend cmd $opt
    }
    lappend cmd $Export(cgnsfile) $Export(outputfile)
    run_command $name $cmd
  }
}

proc export_input {w {base 0} {zone 0} {sol 0} {labw ""}} {
  global Export Font

  if {$labw == ""} {
    if {$sol} {
      set labw 12
    } else {
      set labw 10
    }
  }

  FrameCreate $w.input -text "CGNS Input" -font $Font(bold)
  pack $w.input -side top -padx 5 -pady 2 -fill x
  set input [FrameGet $w.input]

  set f [frame $input.file]
  pack $f -side top -fill x
  label $f.lab -text Filename -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Export(cgnsfile) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 -command "export_browse $w cgnsfile"
  pack $f.but -side right -fill y

  if {$base} {
    set f [frame $input.base]
    pack $f -side top -anchor w
    label $f.lab -text "Base Index" -width $labw -anchor w
    entry $f.ent -textvariable Export(basenum) -width 10
    checkbutton $f.but -text default \
      -variable Export(basenum,flag) -onvalue 1 -offvalue 0 \
      -command "export_state basenum $f.ent"
    pack $f.lab $f.ent $f.but -side left
    export_state basenum $f.ent
  }

  if {$zone} {
    set f [frame $input.zone]
    pack $f -side top -anchor w
    label $f.lab -text "Zone Index" -width $labw -anchor w
    entry $f.ent -textvariable Export(zonenum) -width 10
    checkbutton $f.but -text default \
      -variable Export(zonenum,flag) -onvalue 1 -offvalue 0 \
      -command "export_state zonenum $f.ent"
    pack $f.lab $f.ent $f.but -side left
    export_state zonenum $f.ent
  }

  if {$sol} {
    set f [frame $input.sol]
    pack $f -side top -anchor w
    label $f.lab -text "Solution Index" -width $labw -anchor w
    entry $f.ent -textvariable Export(solnum) -width 10
    checkbutton $f.but -text default \
      -variable Export(solnum,flag) -onvalue 1 -offvalue 0 \
      -command "export_state solnum $f.ent"
    pack $f.lab $f.ent $f.but -side left
    export_state solnum $f.ent
  }
}

proc export_output {w ftype name exts {labw 10}} {
  global Export Font

  FrameCreate $w.output -text "$name Output" -font $Font(bold)
  pack $w.output -side top -padx 5 -pady 2 -fill x
  set f [FrameGet $w.output]

  label $f.lab -text Filename -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Export($ftype) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 \
    -command "export_browse $w $ftype {$name} {$exts}"
  pack $f.but -side right -fill y
}

proc export_options {w {labw ""}} {
  global Export Font

  FrameCreate $w.opts -text "Command Line Options" -font $Font(bold)
  pack $w.opts -side top -padx 5 -pady 2 -fill x
  set f [FrameGet $w.opts]

  label $f.lab -text Options
  pack $f.lab -side left
  if {$labw != ""} {
    $f.lab configure -width $labw -anchor w
  }
  entry $f.ent -textvariable Export(options) -width 30
  pack $f.ent -side left -fill x -expand 1
}

proc export_browse {w ftype {name ""} {exts {}}} {
  global Export tcl_platform
  if {$ftype == "cgnsfile"} {
    set fname [FileOpen "CGNS Input File" $Export(cgnsfile) $w \
      {{{CGNS Files} {.cgns .cgn .adf}} {{All Files} {*}}}]
    if {$fname != ""} {
      if {$tcl_platform(platform) == "windows"} {
        set Export(cgnsfile) [join [split $fname /] \\]
      } else {
        set Export(cgnsfile) $fname
      }
    }
  } elseif {$ftype == "exefile"} {
    if {$tcl_platform(platform) == "windows"} {
      set types {{{Executable Files} {.exe .cmd .bat .com}} {{All Files} {*}}}
    } else {
      set types {{{All Files} {*}}}
    }
    set fname [FileOpen "Export Executable" $Export(exefile) $w $types]
    if {$fname != ""} {
      if {$tcl_platform(platform) == "windows"} {
        set Export(exefile) [join [split $fname /] \\]
      } else {
        set Export(exefile) $fname
      }
    }
  } else {
    if {$name == ""} {set name Output}
    if {$exts == {}} {
      set types {{{All Files} {*}}}
    } else {
      set types [list [list "$name Files" $exts] [list "All Files" *]]
    }
    set fname [FileSave "$name File" $Export($ftype) $w $types]
    if {$fname != ""} {
      if {$tcl_platform(platform) == "windows"} {
        set Export($ftype) [join [split $fname /] \\]
      } else {
        set Export($ftype) $fname
      }
    }
  }
}

proc export_state {what ent {def 1}} {
  global Export
  if {$Export($what,flag)} {
    set Export($what) ""
    $ent configure -state disabled -cursor {}
  } else {
    if {$Export($what) == ""} {
      set Export($what) $def
    }
    $ent configure -state normal -cursor xterm
  }
}

proc export_buttons {w {check ""}} {
  global Export
  set f [frame $w.but]
  pack $f -side top -pady 5
  button $f.accept -text Accept -width 6 -default active
  button $f.cancel -text Cancel -width 6 -command {set Export(done) 0}
  pack $f.accept $f.cancel -side left -padx 5

  if {$check == ""} {
    $f.accept configure -command {set Export(done) 1}
  } else {
    $f.accept configure -command "$check $w"
    bind $w <Return> "$check $w"
  }

  center_window $w .
  set oldFocus [focus]
  set oldGrab [grab current $w]
  if {$oldGrab != ""} {
    set grabStatus [grab status $oldGrab]
  }
  catch {grab $w}
  tkwait visibility $w
  focus $w
  tkwait variable Export(done)
  catch {focus $oldFocus}
  destroy $w
  if {$oldGrab != ""} {
    if {$grabStatus == "global"} {
      grab -global $oldGrab
    } else {
      grab $oldGrab
    }
  }
  return $Export(done)
}

proc export_check {w} {
  global Export
  if {$Export(exefile) == "" ||
      $Export(outputfile) == "" ||
      $Export(cgnsfile) == ""} {
    errormsg "must specify an executable, CGNS file and an output file" $w
    return
  }
  if {![file exists $Export(exefile)] ||
      ![file executable $Export(exefile)]} {
    errormsg "the export executable does not exist or is not executable" $w
    return
  }
  if {![file exists $Export(cgnsfile)]} {
    errormsg "the CGNS input file doesn't exist" $w
    return
  }
  set Export(done) 1
}

