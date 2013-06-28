# convert.tcl - CGNS conversion routines

array set Tools {
  cnvfile ""
  location v
  variable p
  varcnv ""
  varcnv,flag 0
  rmini ""
  rmaxi ""
  rminj ""
  rmaxj ""
  rmink ""
  rmaxk ""
  gamma 1.4
  dataclass d
  datacnv ""
  datacnv,flag 0
  datacnv,errs 0
}

array set ConvFiles {
  p primitive.cnv
  c conserved.cnv
  d dimensional.cnv
  n nondimensional.cnv
  u undimensional.cnv
}

array set ConvLabels {
  p Primitive
  c Conserved
  d Dimensional
  n NormalizedByDimensional
  u NormalizedByUnknownDimensional
}

proc convert_location {w name exe} {
  global ProgData Tools Font

  if {$ProgData(file,name) == ""} return
  set cmd [get_executable $exe 1]
  if {$cmd == ""} return

  set Tools(cgnsinput) $ProgData(file,name)
  if {$Tools(cgnsoutput) == ""} {
    set Tools(cgnsoutput) $ProgData(file,name)
  }

  toplevel $w
  wm title $w "Solution Location"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Tools(done) 0}

  FrameCreate $w.loc -text Location -font $Font(bold)
  pack $w.loc -side top -pady 2 -fill x
  set f [FrameGet $w.loc]

  radiobutton $f.v -text "Vertex" \
    -variable Tools(location) -value v \
    -command "catch {pack forget $w.rind}"
  radiobutton $f.c -text "Cell Center" \
    -variable Tools(location) -value c \
    -command "catch {pack $w.rind -side top -pady 2 -fill x}"
  pack $f.v $f.c -side left -expand 1

  tools_cgnsinput $w 2
  tools_cgnsoutput $w solname
  tools_averaging $w

  FrameCreate $w.rind -text "Rind Cells" -font $Font(bold)
  if {$Tools(location) == "c"} {
    pack $w.rind -side top -pady 2 -fill x
  }
  set rind [FrameGet $w.rind]

  foreach i {i j k} {
    set f [frame $rind.$i]
    pack $f -side left -expand 1
    checkbutton $f.min -text "min $i" \
      -variable Tools(rmin$i) -onvalue $i -offvalue ""
    checkbutton $f.max -text "max $i" \
      -variable Tools(rmax$i) -onvalue [string toupper $i] -offvalue ""
    pack $f.min $f.max -side top -anchor w
  }

  if {![tools_interact $w]} return

  set opts -$Tools(location)$Tools(weight)
  if {$Tools(location) == "c"} {
    foreach i {i j k} {
      foreach j {rmin rmax} {
        if {$Tools($j$i) != ""} {
          append opts $Tools($j$i)
        }
      }
    }
  }
  lappend cmd $opts
  foreach i {basenum zonenum solnum} {
    if {$Tools($i) != ""} {
      lappend cmd "-[string index $i 0]$Tools($i)"
    }
  }
  if {$Tools(solname) != ""} {
    lappend cmd "-S$Tools(solname)"
  }
  lappend cmd $Tools(cgnsinput)
  if {$Tools(cgnsoutput) != ""} {
    lappend cmd $Tools(cgnsoutput)
  }

  if {$Tools(location) == "c"} {
    tools_run "Vertex to CellCenter" $cmd $Tools(cgnsoutput)
  } else {
    tools_run "CellCenter to Vertex" $cmd $Tools(cgnsoutput)
  }
}

proc convert_variables {w name exe} {
  global ProgData Tools Font ConvLabels

  if {$ProgData(file,name) == ""} return
  set cmd [get_executable $exe 1]
  if {$cmd == ""} return

  set Tools(cgnsinput) $ProgData(file,name)
  if {$Tools(cgnsoutput) == ""} {
    set Tools(cgnsoutput) $ProgData(file,name)
  }

  toplevel $w
  wm title $w "Solution Variables"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Tools(done) 0}

  FrameCreate $w.var -text Variables -font $Font(bold)
  pack $w.var -side top -pady 2 -fill x
  set var [FrameGet $w.var]

  set f [frame $var.f1]
  pack $f -side top -fill x
  radiobutton $f.p -text $ConvLabels(p) -variable Tools(variable) \
    -value p -command {convert_file p varcnv}
  radiobutton $f.c -text $ConvLabels(c) -variable Tools(variable) \
    -value c -command {convert_file c varcnv}
  frame $f.g
  label $f.g.lab -text gamma
  entry $f.g.ent -textvariable Tools(gamma) -width 10
  pack $f.g.lab $f.g.ent -side left
  pack $f.p $f.c $f.g -side left -expand 1

  frame $var.sep -bd 1 -height 2 -relief sunken
  pack $var.sep -side top -fill x -pady 3

  if {$Tools(varcnv) == ""} {
    convert_file $Tools(variable) varcnv
  }

  set f [frame $var.f2]
  pack $f -side top -fill x
  checkbutton $f.use -text "Use Conversion File" \
    -variable Tools(varcnv,flag) -onvalue 1 -offvalue 0
  button $f.edit -text "Edit Conversion File" \
    -command "convert_edit varcnv $w"
  pack $f.use $f.edit -side left -expand 1

  set f [frame $var.f3]
  pack $f -side top -fill x
  label $f.lab -text Filename -width 13 -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Tools(varcnv) -width 20
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -command "convert_browse varcnv $f.but"
  pack $f.but -side right

  tools_cgnsinput $w 2
  tools_cgnsoutput $w solname

  if {![tools_interact $w check_variables]} return

  if {$Tools(varcnv,flag)} {
    lappend cmd -f$Tools(cnvfile)
  } else {
    lappend cmd -$Tools(variable)
  }
  foreach i {basenum zonenum solnum} {
    if {$Tools($i) != ""} {
      lappend cmd "-[string index $i 0]$Tools($i)"
    }
  }
  if {$Tools(solname) != ""} {
    lappend cmd "-S$Tools(solname)"
  }
  lappend cmd $Tools(cgnsinput)
  if {$Tools(cgnsoutput) != ""} {
    lappend cmd $Tools(cgnsoutput)
  }

  set n $Tools(variable)
  tools_run "$ConvLabels($n) Variables" $cmd $Tools(cgnsoutput)
}

proc convert_file {n ftype} {
  global Tools ConvFiles tcl_platform
  set name $ConvFiles($n)
  set cnvpath [get_file $name]
  if {$cnvpath == ""} {
    set Tools($ftype) $name
  } else {
    if {$tcl_platform(platform) == "windows"} {
      set Tools($ftype) [join [split $cnvpath /] \\]
    } else {
      set Tools($ftype) $cnvpath
    }
  }
}

proc convert_edit {name {parent .}} {
  global Tools _EditFile
  if {$Tools($name) == ""} return
  edit_file .editcnv $Tools($name) $parent
  tkwait window .editcnv
  set Tools($name) $_EditFile(.editcnv,name)
}

proc convert_browse {name w} {
  global Tools tcl_platform
  set fname [FileOpen "Conversion File" $Tools($name) $w \
    {{{Conversion Scripts} {.cnv .clc .txt}} {{All Files} {*}}}]
  if {$fname != ""} {
    if {$tcl_platform(platform) == "windows"} {
      set Tools($name) [join [split $fname /] \\]
    } else {
      set Tools($name) $fname
    }
  }
}

proc check_variables {w} {
  global Tools tcl_platform
  if {$Tools(varcnv,flag)} {
    if {$Tools(varcnv) == ""} {
      errormsg "variable conversion file not specified" $w
      return 0
    }
    set cnvpath [get_file $Tools(varcnv)]
    if {$cnvpath == ""} {
      errormsg "can't find conversion file $Tools(varcnv)" $w
      return 0
    }
    if {$tcl_platform(platform) == "windows"} {
      set Tools(cnvfile) [join [split $cnvpath /] \\]
    } else {
      set Tools(cnvfile) $cnvpath
    }
  }
  if {$Tools(gamma) != ""} {
    if {[catch {expr $Tools(gamma) <= 1.0} nogood] || $nogood} {
      errormsg "invalid value for gamma" $w
      return 0
    }
  }
  return [tools_check_input $w]
}

proc convert_dimensional {w name exe} {
  global ProgData Tools Font ConvLabels

  if {$ProgData(file,name) == ""} return
  set cmd [get_executable $exe 1]
  if {$cmd == ""} return

  set Tools(cgnsinput) $ProgData(file,name)
  if {$Tools(cgnsoutput) == ""} {
    set Tools(cgnsoutput) $ProgData(file,name)
  }

  toplevel $w
  wm title $w "Dimensionalization"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Tools(done) 0}

  FrameCreate $w.var -text "Data Class" -font $Font(bold)
  pack $w.var -side top -pady 2 -fill x
  set var [FrameGet $w.var]

  set f [frame $var.f1]
  pack $f -side top -fill x
  foreach i {d n u} {
    radiobutton $f.$i -text $ConvLabels($i) -variable Tools(dataclass) \
      -value $i -command "convert_file $i datacnv"
    pack $f.$i -side left -expand 1
  }

  frame $var.sep -bd 1 -height 2 -relief sunken
  pack $var.sep -side top -fill x -pady 3

  if {$Tools(datacnv) == ""} {
    convert_file $Tools(dataclass) datacnv
  }

  set f [frame $var.f2]
  pack $f -side top -fill x
  checkbutton $f.use -text "Use Conversion File" \
    -variable Tools(datacnv,flag) -onvalue 1 -offvalue 0
  checkbutton $f.errs -text "Show Parsing Errors" \
    -variable Tools(datacnv,errs) -onvalue 1 -offvalue 0
  button $f.edit -text "Edit Conversion File" \
    -command "convert_edit datacnv $w"
  pack $f.use $f.errs $f.edit -side left -expand 1

  set f [frame $var.f3]
  pack $f -side top -fill x
  label $f.lab -text "Filename" -width 13 -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Tools(datacnv) -width 20
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -command "convert_browse datacnv $f.but"
  pack $f.but -side right

  tools_cgnsinput $w 2
  tools_cgnsoutput $w

  if {![tools_interact $w check_dimensional]} return

  lappend cmd -$Tools(dataclass)
  if {$Tools(datacnv,errs)} {
    lappend cmd -v
  }
  if {$Tools(datacnv,flag)} {
    lappend cmd -f$Tools(cnvfile)
  }
  foreach i {basenum zonenum solnum} {
    if {$Tools($i) != ""} {
      lappend cmd "-[string index $i 0]$Tools($i)"
    }
  }
  lappend cmd $Tools(cgnsinput)
  if {$Tools(cgnsoutput) != ""} {
    lappend cmd $Tools(cgnsoutput)
  }

  set n $Tools(dataclass)
  tools_run $ConvLabels($n) $cmd $Tools(cgnsoutput)
}

proc check_dimensional {w} {
  global Tools tcl_platform
  if {$Tools(datacnv,flag)} {
    if {$Tools(datacnv) == ""} {
      errormsg "data class conversion file not specified" $w
      return 0
    }
    set cnvpath [get_file $Tools(datacnv)]
    if {$cnvpath == ""} {
      errormsg "can't find conversion file $Tools(datacnv)" $w
      return 0
    }
    if {$tcl_platform(platform) == "windows"} {
      set Tools(cnvfile) [join [split $cnvpath /] \\]
    } else {
      set Tools(cnvfile) $cnvpath
    }
  }
  return [tools_check_input $w]
}

