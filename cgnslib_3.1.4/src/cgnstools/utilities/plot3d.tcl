# plot3d.tcl - PLOT3D import/export

array set Plot3d {
  mode ""
  cgns ""
  xyz ""
  q ""
  basenum ""
  basename ""
  solnum ""
  basenum,flag 1
  solnum,flag 1
  block ""
  format ""
  iblank ""
  type u
  mach ""
  ignore ""
  double ""
  base ""
  convert ""
  weight ""
  gamma 1.4
}

proc plot3d_import {w name exe} {
  global ProgData Plot3d Import
  set cmd [get_executable $exe 1]
  if {$cmd == ""} return

  set Plot3d(mode) import
  set Plot3d(cgns) $ProgData(file,name)
  if ![plot3d_panel $w] return
  update idletasks

  set opts $Plot3d(block)$Plot3d(format)$Plot3d(type)
  if {$Plot3d(ignore) != ""} {
    append opts $Plot3d(ignore)
  } else {
    append opts $Plot3d(iblank)
  }
  append opts $Plot3d(double)$Plot3d(convert)
  if {$Plot3d(mach) != ""} {
    append opts M$Plot3d(mach)
  }
  if {$opts != ""} {
    lappend cmd -$opts
  }
  if {$Plot3d(basenum) != ""} {
    lappend cmd -b$Plot3d(basenum)
  }
  if {$Plot3d(basename) != ""} {
    lappend cmd -B$Plot3d(basename)
  }
  if {$Plot3d(gamma) != ""} {
    lappend cmd -g$Plot3d(gamma)
  }

  lappend cmd $Plot3d(xyz)
  if {$Plot3d(q) != ""} {
    lappend cmd $Plot3d(q)
  }
  lappend cmd $Plot3d(cgns)

  import_run "PLOT3D Import" $cmd $Plot3d(cgns)
}

proc plot3d_export {w name exe} {
  global ProgData Plot3d
  set cmd [get_executable $exe 1]
  if {$cmd == ""} return

  set Plot3d(mode) export
  set Plot3d(cgns) $ProgData(file,name)
  set Plot3d(xyz) [file rootname $ProgData(file,name)].xyz
  set Plot3d(q) [file rootname $ProgData(file,name)].q
  if ![plot3d_panel $w] return
  update idletasks

  set opts $Plot3d(block)$Plot3d(format)$Plot3d(type)
  append opts $Plot3d(ignore)$Plot3d(double)$Plot3d(weight)
  if {$opts != ""} {
    lappend cmd -$opts
  }
  if {$Plot3d(basenum) != ""} {
    lappend cmd -b$Plot3d(basenum)
  }
  if {$Plot3d(solnum) != ""} {
    lappend cmd -S$Plot3d(solnum)
  }
  if {$Plot3d(gamma) != ""} {
    lappend cmd -g$Plot3d(gamma)
  }

  lappend cmd $Plot3d(cgns) $Plot3d(xyz)
  if {$Plot3d(q) != ""} {
    lappend cmd $Plot3d(q)
  }

  update
  run_command "PLOT3D Export" $cmd
}

proc plot3d_panel {w} {
  global ProgData Font Plot3d

  toplevel $w
  wm title $w "Plot3d $Plot3d(mode)"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Plot3d(done) 0}

  FrameCreate $w.file -text "Files" -font $Font(bold)
  pack $w.file -side top -pady 2 -fill x
  set file [FrameGet $w.file]

  foreach i {cgns xyz q} {
    set f [frame $file.$i]
    pack $f -side top -fill x
    label $f.lab -text [string toupper $i] -width 6 -anchor w
    pack $f.lab -side left
    entry $f.ent -textvariable Plot3d($i) -width 30
    pack $f.ent -side left -fill x -expand 1
    button $f.but -text Browse -pady 0 -command "plot3d_browse_$i $w"
    pack $f.but -side right -fill y
  }

  FrameCreate $w.base -text "CGNS Location" -font $Font(bold)
  pack $w.base -side top -pady 2 -fill x
  set base [FrameGet $w.base]
  set labw 12

  set f [frame $base.num]
  pack $f -side top -anchor w
  label $f.lab -text "Base Index" -width $labw -anchor w
  entry $f.ent -textvariable Plot3d(basenum) -width 10
  checkbutton $f.but -text default \
    -variable Plot3d(basenum,flag) -onvalue 1 -offvalue 0 \
    -command "plot3d_state basenum $f.ent"
  pack $f.lab $f.ent $f.but -side left
  plot3d_state basenum $f.ent

  set f [frame $base.name]
  if {$Plot3d(mode) == "import"} {
    pack $f -side top -anchor w
  }
  label $f.lab -text "Base Name" -width $labw -anchor w
  entry $f.ent -textvariable Plot3d(basename) -width 30
  pack $f.lab $f.ent -side left

  set f [frame $base.sol]
  if {$Plot3d(mode) == "export"} {
    pack $f -side top -anchor w
  }
  label $f.lab -text "Solution Index" -width $labw -anchor w
  entry $f.ent -textvariable Plot3d(solnum) -width 10
  checkbutton $f.but -text default \
    -variable Plot3d(solnum,flag) -onvalue 1 -offvalue 0 \
    -command "plot3d_state solnum $f.ent"
  pack $f.lab $f.ent $f.but -side left
  plot3d_state solnum $f.ent

  FrameCreate $w.format -text "Plot3d File Format" -font $Font(bold)
  pack $w.format -side top -pady 2 -fill x
  set fmt [FrameGet $w.format]

  set f [frame $fmt.block]
  pack $f -side left -expand 1
  radiobutton $f.m -text "Multi-Block" \
    -variable Plot3d(block) -value ""
  radiobutton $f.s -text "Single-Block" \
    -variable Plot3d(block) -value s
  pack $f.m $f.s -side top -anchor w

  set f [frame $fmt.format]
  pack $f -side left -expand 1
  radiobutton $f.w -text "Whole Format" \
    -variable Plot3d(format) -value ""
  radiobutton $f.p -text "Planar Format" \
    -variable Plot3d(format) -value p
  pack $f.w $f.p -side top -anchor w

  set f [frame $fmt.iblank]
  if {$Plot3d(mode) == "import"} {
    pack $f -side left -expand 1
  }
  radiobutton $f.n -text "No iblank" \
    -variable Plot3d(iblank) -value ""
  radiobutton $f.y -text "Has iblank" \
    -variable Plot3d(iblank) -value i
  pack $f.n $f.y -side top -anchor w

  FrameCreate $w.type -text "Plot3d File Type" -font $Font(bold)
  pack $w.type -side top -pady 2 -fill x
  set type [FrameGet $w.type]

  radiobutton $type.b -text binary \
    -variable Plot3d(type) -value ""
  radiobutton $type.u -text unformatted \
    -variable Plot3d(type) -value u
  radiobutton $type.f -text formatted \
    -variable Plot3d(type) -value f
  checkbutton $type.d -text "64-bit" \
    -variable Plot3d(double) -onvalue d -offvalue ""
  pack $type.b $type.u $type.f $type.d -side left -expand 1

  FrameCreate $w.mach -text "Plot3d Machine Type" -font $Font(bold)
  if {$Plot3d(mode) == "import"} {
    pack $w.mach -side top -pady 2 -fill x
  }
  set mach [FrameGet $w.mach]

  set n 1
  foreach type {\
    {default IEEE BSIEEE} \
    {Iris Sun HP} \
    {Alpha DEC IBM} \
    {Cray Convex NT}} {
    set f [frame $mach.f$n]
    pack $f -side left -expand 1
    foreach i $type {
      set j [string tolower $i]
      radiobutton $f.$j -text $i \
        -variable Plot3d(mach) -value $j
      if {$i == "default"} {
        $f.$j configure -value ""
      }
      pack $f.$j -side top -anchor w
    }
    incr n
  }

  FrameCreate $w.opts -text "Options" -font $Font(bold)
  pack $w.opts -side top -pady 2 -fill x
  set opts [FrameGet $w.opts]

  set f [frame $opts.f1]
  pack $f -side left -expand 1
  if {$Plot3d(mode) == "import"} {
    checkbutton $f.i -text "Don't use iblank" \
      -variable Plot3d(ignore) -onvalue n -offvalue ""
    checkbutton $f.c -text "Write primitive variables" \
      -variable Plot3d(convert) -onvalue c -offvalue ""
  } else {
    checkbutton $f.i -text "Don't write iblank" \
      -variable Plot3d(ignore) -onvalue n -offvalue ""
    checkbutton $f.c -text "Use volume weighting" \
      -variable Plot3d(weight) -onvalue w -offvalue ""
  }
  pack $f.i $f.c -side top -anchor w

  set f [frame $opts.gamma]
  pack $f -side left -expand 1
  label $f.lab -text gamma
  entry $f.ent -textvariable Plot3d(gamma) -width 10
  pack $f.lab $f.ent -side left

  set f [frame $w.but]
  pack $f -side top -pady 5
  button $f.accept -text Accept -width 6 -default active \
    -command "plot3d_check $w"
  button $f.cancel -text Cancel -width 6 -command {set Plot3d(done) 0}
  pack $f.accept $f.cancel -side left -padx 5

  bind $w <Return> "plot3d_check $w"

  center_window $w .
  set oldFocus [focus]
  set oldGrab [grab current $w]
  if {$oldGrab != ""} {
    set grabStatus [grab status $oldGrab]
  }
  catch {grab $w}
  tkwait visibility $w
  focus $w
  tkwait variable Plot3d(done)
  catch {focus $oldFocus}
  destroy $w
  if {$oldGrab != ""} {
    if {$grabStatus == "global"} {
      grab -global $oldGrab
    } else {
      grab $oldGrab
    }
  }

  return $Plot3d(done)
}

proc plot3d_state {what ent {def 1}} {
  global Plot3d
  if {$Plot3d($what,flag)} {
    set Plot3d($what) ""
    $ent configure -state disabled -cursor {}
  } else {
    if {$Plot3d($what) == ""} {
      set Plot3d($what) $def
    }
    $ent configure -state normal -cursor xterm
  }
}

proc plot3d_check {w} {
  global Plot3d
  if {[string trim $Plot3d(cgns)] == "" ||
      [string trim $Plot3d(xyz)] == ""} {
    errormsg "must specify a CGNS and a Plot3d XYZ file" $w
    return
  }
  if {$Plot3d(mode) == "import"} {
    if {![file exists $Plot3d(xyz)]} {
      errormsg "Plot3d XYZ file doesn't exist" $w
      return
    }
    if {$Plot3d(q) != "" && ![file exists $Plot3d(q)]} {
      errormsg "Plot3d Q file doesn't exist" $w
      return
    }
  } else {
    if {![file exists $Plot3d(cgns)]} {
      errormsg "CGNS file doesn't exist" $w
      return
    }
  }
  if {$Plot3d(gamma) != ""} {
    if {[catch {expr $Plot3d(gamma) <= 1.0} nogood] || $nogood} {
      errormsg "invalid value for gamma" $w
      return
    }
  }
  set Plot3d(done) 1
}

proc plot3d_browse_cgns {w} {
  global Plot3d tcl_platform
  if {$Plot3d(mode) == "import"} {
    set fname [FileSave "CGNS Output File" $Plot3d(cgns) $w \
      {{{CGNS Files} {.cgns .cgn .adf}} {{All Files} {*}}} cgns]
  } else {
    set fname [FileOpen "CGNS Input File" $Plot3d(cgns) $w \
      {{{CGNS Files} {.cgns .cgn .adf}} {{All Files} {*}}}]
  }
  if {$fname != ""} {
    if {$tcl_platform(platform) == "windows"} {
      set Plot3d(cgns) [join [split $fname /] \\]
    } else {
      set Plot3d(cgns) $fname
    }
  }
}

proc plot3d_browse_xyz {w} {
  global Plot3d tcl_platform
  if {$Plot3d(mode) == "import"} {
    set fname [FileOpen "Plot3d XYZ Input File" $Plot3d(xyz) $w \
      {{{Plot3d Files} {.bin .dat .fmt .xyz}} {{All Files} {*}}}]
  } else {
    set fname [FileSave "Plot3d XYZ Output File" $Plot3d(xyz) $w \
      {{{Plot3d Files} {.bin .dat .fmt .xyz}} {{All Files} {*}}}]
  }
  if {$fname != ""} {
    if {$tcl_platform(platform) == "windows"} {
      set Plot3d(xyz) [join [split $fname /] \\]
    } else {
      set Plot3d(xyz) $fname
    }
  }
}

proc plot3d_browse_q {w} {
  global Plot3d tcl_platform
  if {$Plot3d(mode) == "import"} {
    set fname [FileOpen "Plot3d Q Input File" $Plot3d(q) $w \
      {{{Plot3d Files} {.bin .dat .fmt .q}} {{All Files} {*}}}]
  } else {
    set fname [FileSave "Plot3d Q Output File" $Plot3d(q) $w \
      {{{Plot3d Files} {.bin .dat .fmt .q}} {{All Files} {*}}}]
  }
  if {$fname != ""} {
    if {$tcl_platform(platform) == "windows"} {
      set Plot3d(q) [join [split $fname /] \\]
    } else {
      set Plot3d(q) $fname
    }
  }
}

