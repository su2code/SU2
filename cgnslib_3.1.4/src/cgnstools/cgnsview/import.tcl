# import.tcl - CGNS import routines

proc import_default {w name cmd} {
  global ProgData Font Import tcl_platform

  if {$cmd != ""} {
    set exe [lindex $cmd 0]
    set exepath [get_executable $exe]
    if {$exepath != ""} {set exe $exepath}
    if {$tcl_platform(platform) == "windows"} {
      set Import(exefile) [join [split $exe /] \\]
    } else {
      set Import(exefile) $exe
    }
    set Import(options) [lrange $cmd 1 end]
  }
  set Import(cgnsfile) $ProgData(file,name)

  toplevel $w
  wm title $w $name
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Import(done) 0}

  set labw 12

  FrameCreate $w.top -text "Import Command" -font $Font(bold)
  pack $w.top -side top -padx 5 -pady 5 -fill x
  set top [FrameGet $w.top]

  set f [frame $top.exe]
  pack $f -side top -fill x -expand 1
  label $f.lab -text Executable -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Import(exefile) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 -command "import_browse $w exefile"
  pack $f.but -side right -fill y
  $f.ent xview [$f.ent index end]

  set f [frame $top.opts]
  pack $f -side top -fill x -expand 1
  label $f.lab -text Options -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Import(options) -width 30
  pack $f.ent -side left -fill x -expand 1

  set f [frame $top.input]
  pack $f -side top -fill x -expand 1
  label $f.lab -text "Input File" -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Import(inputfile) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 \
    -command "import_browse $w inputfile {$name}"
  pack $f.but -side right -fill y

  set f [frame $top.output]
  pack $f -side top -fill x -expand 1
  label $f.lab -text "CGNS Output" -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Import(cgnsfile) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 -command "import_browse $w cgnsfile"
  pack $f.but -side right -fill y

  if {[import_buttons $w import_check]} {
    set cmd [tools_unix_path $Import(exefile)]
    foreach opt $Import(options) {
      lappend cmd $opt
    }
    lappend cmd $Import(inputfile) $Import(cgnsfile)
    import_run $name $cmd $Import(cgnsfile)
  }
}

proc import_input {w ftype name exts {labw 10}} {
  global Import Font

  FrameCreate $w.input -text "$name Input" -font $Font(bold)
  pack $w.input -side top -padx 5 -pady 2 -fill x
  set f [FrameGet $w.input]

  label $f.lab -text Filename -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Import($ftype) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 \
    -command "import_browse $w $ftype {$name} {$exts}"
  pack $f.but -side right -fill y
}

proc import_output {w {use_base 0} {labw 10}} {
  global Import Font

  FrameCreate $w.output -text "CGNS Output" -font $Font(bold)
  pack $w.output -side top -padx 5 -pady 2 -fill x
  set output [FrameGet $w.output]

  set f [frame $output.file]
  pack $f -side top -fill x
  label $f.lab -text Filename -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Import(cgnsfile) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 -command "import_browse $w cgnsfile"
  pack $f.but -side right -fill y

  if {$use_base} {
    set f [frame $output.basename]
    pack $f -side top -anchor w
    label $f.lab -text "Base Name" -width $labw -anchor w
    entry $f.ent -textvariable Import(basename) -width 30
    checkbutton $f.but -text default \
      -variable Import(basename,flag) -onvalue 1 -offvalue 0 \
      -command "import_set_state basename $f.ent"
    pack $f.lab $f.ent $f.but -side left
    import_set_state basename $f.ent
  }
}

proc import_options {w {labw ""}} {
  global Import Font

  FrameCreate $w.opts -text "Command Line Options" -font $Font(bold)
  pack $w.opts -side top -padx 5 -pady 2 -fill x
  set f [FrameGet $w.opts]

  label $f.lab -text Options
  if {$labw != ""} {
    $f.lab configure -width $labw -anchor w
  }
  pack $f.lab -side left
  entry $f.ent -textvariable Import(options) -width 30
  pack $f.ent -side left -fill x -expand 1
}

proc import_browse {w ftype {name ""} {exts {}}} {
  global Import tcl_platform
  if {$ftype == "cgnsfile"} {
    set fname [FileSave "CGNS Output File" $Import(cgnsfile) $w \
      {{{CGNS Files} {.cgns .cga .cgh .cgx}} {{All Files} {*}}} cgns]
    if {$fname != ""} {
      if {$tcl_platform(platform) == "windows"} {
        set Import(cgnsfile) [join [split $fname /] \\]
      } else {
        set Import(cgnsfile) $fname
      }
    }
  } elseif {$ftype == "exefile"} {
    if {$tcl_platform(platform) == "windows"} {
      set types {{{Executable Files} {.exe .cmd .bat .com}} {{All Files} {*}}}
    } else {
      set types {{{All Files} {*}}}
    }
    set fname [FileOpen "Import Executable" $Import(exefile) $w $types]
    if {$fname != ""} {
      if {$tcl_platform(platform) == "windows"} {
        set Import(exefile) [join [split $fname /] \\]
      } else {
        set Import(exefile) $fname
      }
    }
  } else {
    if {$name == ""} {set name Input}
    if {$exts == {}} {
      set types {{{All Files} {*}}}
    } else {
      set types [list [list "$name Files" $exts] [list "All Files" *]]
    }
    set fname [FileOpen "$name File" $Import($ftype) $w $types]
    if {$fname != ""} {
      if {$tcl_platform(platform) == "windows"} {
        set Import($ftype) [join [split $fname /] \\]
      } else {
        set Import($ftype) $fname
      }
      set Import(cgnsfile) [file rootname $Import($ftype)].cgns
    }
  }
}

proc import_check {w} {
  global Import
  if {$Import(exefile) == "" ||
      $Import(inputfile) == "" ||
      $Import(cgnsfile) == ""} {
    errormsg "must specify an executable, input file and a CGNS file" $w
    return
  }
  if {![file exists $Import(exefile)] ||
      ![file executable $Import(exefile)]} {
    errormsg "the import executable does not exist or is not executable" $w
    return
  }
  if {![file exists $Import(inputfile)]} {
    errormsg "the input file doesn't exist" $w
    return
  }
  set Import(done) 1
}

proc import_buttons {w {check ""}} {
  global Import
  set f [frame $w.but]
  pack $f -side top -pady 5
  button $f.accept -text Accept -width 6 -default active
  button $f.cancel -text Cancel -width 6 -command {set Import(done) 0}
  pack $f.accept $f.cancel -side left -padx 5

  if {$check == ""} {
    $f.accept configure -command {set Import(done) 1}
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
  tkwait variable Import(done)
  catch {focus $oldFocus}
  destroy $w
  if {$oldGrab != ""} {
    if {$grabStatus == "global"} {
      grab -global $oldGrab
    } else {
      grab $oldGrab
    }
  }
  return $Import(done)
}

proc import_run {title cmd {outfile ""}} {
  global ProgData
  if {$ProgData(file,name) == ""} {
    set same 0
  } elseif {$outfile == ""} {
    set same 1
  } elseif {[file dirname $outfile] == "."} {
    set same [same_file $outfile $ProgData(file,name)]
  } else {
    set same [same_file $outfile $ProgData(file)]
  }
  if {$same} {
    catch CGIOclose
    do_backup
  }
  update
  if {[run_command $title $cmd] && !$same} {
    file_load $outfile
  } else {
    file_reload
  }
}

proc import_dup_check {f} {
  global Import

  checkbutton $f.dup -text "Duplicate Node Checking" \
    -variable Import(dupcheck,flag) -onvalue 1 -offvalue 0
  pack $f.dup -side top -anchor w

  frame $f.tol
  pack $f.tol -side left
  label $f.tol.lab -text Tolerance
  entry $f.tol.ent -textvariable Import(duptol) -width 10
  pack $f.tol.lab $f.tol.ent -side top

  frame $f.type
  pack $f.type -side left
  radiobutton $f.type.rel -text Relative -variable Import(dupcheck) -value d
  radiobutton $f.type.abs -text Absolute -variable Import(dupcheck) -value D
  pack $f.type.rel $f.type.abs -side top -anchor w
}

proc import_set_state {what ent} {
  global Import
  if {$Import($what,flag)} {
    set Import($what) ""
    $ent configure -state disabled -cursor {}
  } else {
    $ent configure -state normal -cursor xterm
  }
}

