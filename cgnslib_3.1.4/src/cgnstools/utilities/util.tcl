# util.tcl - CGNS utility routines

array set Tools {
  subseti i
  subsetj j
  subsetk k
  cgnssol ""
  solbase ""
  solbase,flag 1
  bboxpad ""
  bboxpad,flag 1
  maxdepth ""
  maxdepth,flag 1
  maxelem ""
  maxelem,flag 1
  extrapolate 0
}

proc interpolate_cgns {w name exe} {
  global ProgData Tools Font

  if {$ProgData(file,name) == ""} return
  set cmd [get_executable $exe 1]
  if {$cmd == ""} return

  set Tools(cgnsinput) $ProgData(file,name)
  if {$Tools(cgnsoutput) == ""} {
    set Tools(cgnsoutput) $ProgData(file,name)
  }

  toplevel $w
  wm title $w "Interpolate Solution"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Tools(done) 0}

  set labw 13

  FrameCreate $w.cgnssol -text "CGNS Solution" -font $Font(bold)
  pack $w.cgnssol -side top -pady 2 -fill x
  set wf [FrameGet $w.cgnssol]

  set f [frame $wf.file]
  pack $f -side top -fill x
  label $f.lab -text Filename -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Tools(cgnssol) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 -command "tools_browse $f.but cgnssol"
  pack $f.but -side right -fill y

  set f [frame $wf.basenum]
  pack $f -side top -fill x
  label $f.lab -text "Base Index" -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Tools(solbase) -width 10
  pack $f.ent -side left
  checkbutton $f.but -text default \
    -variable Tools(solbase,flag) -onvalue 1 -offvalue 0 \
    -command "tools_state solbase $f.ent"
  pack $f.but -side left

  tools_state solbase $f.ent

  FrameCreate $w.cgnsinput -text "CGNS Grid" -font $Font(bold)
  pack $w.cgnsinput -side top -pady 2 -fill x
  set wf [FrameGet $w.cgnsinput]

  set f [frame $wf.file]
  pack $f -side top -fill x
  label $f.lab -text Filename -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Tools(cgnsinput) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 -command "tools_browse $f.but cgnsinput"
  pack $f.but -side right -fill y

  set f [frame $wf.basenum]
  pack $f -side top -fill x
  label $f.lab -text "Base Index" -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Tools(basenum) -width 10
  pack $f.ent -side left
  checkbutton $f.but -text default \
    -variable Tools(basenum,flag) -onvalue 1 -offvalue 0 \
    -command "tools_state basenum $f.ent"
  pack $f.but -side left

  tools_state basenum $f.ent

  FrameCreate $w.cgnsoutput -text "CGNS Output" -font $Font(bold)
  pack $w.cgnsoutput -side top -pady 2 -fill x
  set wf [FrameGet $w.cgnsoutput]

  set f [frame $wf.file]
  pack $f -side top -fill x
  label $f.lab -text Filename -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Tools(cgnsoutput) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 -command "tools_browse $f.but cgnsoutput"
  pack $f.but -side right -fill y

  set f [frame $wf.basename]
  pack $f -side top -fill x
  label $f.lab -text "Base Name" -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Tools(basename) -width 30
  pack $f.ent -side left -fill x -expand 1
  checkbutton $f.but -text default \
    -variable Tools(basename,flag) -onvalue 1 -offvalue 0 \
    -command "tools_state basename $f.ent {}"
  pack $f.but -side left

  tools_state basename $f.ent ""

  set f [frame $wf.solname]
  pack $f -side top -fill x
  label $f.lab -text "Solution Name" -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Tools(solname) -width 30
  pack $f.ent -side left -fill x -expand 1
  checkbutton $f.but -text default \
    -variable Tools(solname,flag) -onvalue 1 -offvalue 0 \
    -command "tools_state solname $f.ent {}"
  pack $f.but -side left

  tools_state solname $f.ent ""

  tools_averaging $w

  FrameCreate $w.treeopts -text "Octtree Options" -font $Font(bold)
  pack $w.treeopts -side top -pady 2 -fill x
  set treeopts [FrameGet $w.treeopts]

  foreach opt {\
    {bboxpad "Bounding Box Padding" 0.01} \
    {maxdepth "Max Octtree Depth" 16} \
    {maxelem "Max Octtree Elements" 256}} {
    set tag [lindex $opt 0]
    set f [frame $treeopts.$tag]
    pack $f -side top -anchor w
    label $f.lab -text [lindex $opt 1] -width 20 -anchor w
    entry $f.ent -textvariable Tools($tag) -width 15
    checkbutton $f.but -text default \
      -variable Tools($tag,flag) -onvalue 1 -offvalue 0 \
      -command "tools_state $tag $f.ent [lindex $opt 2]"
    pack $f.lab $f.ent $f.but -side left
    tools_state $tag $f.ent [lindex $opt 2]
  }

  checkbutton $treeopts.extrap -text "Allow Extrapolation" \
    -variable Tools(extrapolate) -onvalue 1 -offvalue 0
  pack $treeopts.extrap -side top -anchor w

  if {![tools_interact $w check_interpolation]} return

  foreach i {{c solbase} {b basenum} {B basename} {S solname} \
    {p bboxpad} {d maxdepth} {e maxelem}} {
    set opt $Tools([lindex $i 1])
    if {$opt != ""} {
      lappend cmd "-[lindex $i 0]$opt"
    }
  }
  if {$Tools(weight) != ""} {
    lappend cmd -w
  }
  if {$Tools(extrapolate)} {
    lappend cmd -a
  }
  lappend cmd $Tools(cgnssol) $Tools(cgnsinput)
  if {$Tools(cgnsoutput) != ""} {
    lappend cmd $Tools(cgnsoutput)
  }

  tools_run "Solution Interpolation" $cmd $Tools(cgnsoutput)
}

proc check_interpolation {w} {
  global Tools
  if {$Tools(cgnssol) == "" || ![file exists $Tools(cgnssol)]} {
    errormsg "CGNS solution file not specified or does not exist" $w
    return 0
  }
  if {$Tools(cgnsinput) == "" || ![file exists $Tools(cgnsinput)]} {
    errormsg "CGNS grid file not specified or does not exist" $w
    return 0
  }
  if {$Tools(bboxpad) != "" && [catch {expr double($Tools(bboxpad))}]} {
    errormsg "invalid value for bounding box padding" $w
    return 0
  }
  if {$Tools(maxdepth) != "" && [catch {expr int($Tools(maxdepth))}]} {
    errormsg "invalid value for max octtree depth" $w
    return 0
  }
  if {$Tools(maxelem) != "" && [catch {expr int($Tools(maxelem))}]} {
    errormsg "invalid value for max octtree elements" $w
    return 0
  }
  return 1
}

proc extract_subset {w name exe} {
  global ProgData Tools Font

  if {$ProgData(file,name) == ""} return
  set cmd [get_executable $exe 1]
  if {$cmd == ""} return

  set Tools(cgnsinput) $ProgData(file,name)
  if {$Tools(cgnsoutput) == ""} {
    set Tools(cgnsoutput) $ProgData(file,name)
  }

  toplevel $w
  wm title $w "Extract Subset"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Tools(done) 0}

  FrameCreate $w.loc -text Directions -font $Font(bold)
  pack $w.loc -side top -pady 2 -fill x
  set f [FrameGet $w.loc]

  foreach n {i j k} {
    checkbutton $f.$n -text "[string toupper $n]-subset" \
    -variable Tools(subset$n) -onvalue $n -offvalue ""
    pack $f.$n -side left -expand 1
  }

  tools_cgnsinput $w 1
  tools_cgnsoutput $w basename
  tools_averaging $w

  if {![tools_interact $w check_subset]} return

  set opts $Tools(weight)
  foreach i {i j k} {
    if {$Tools(subset$i) != ""} {
      append opts $Tools(subset$i)
    }
  }
  lappend cmd -$opts
  foreach i {basenum zonenum} {
    if {$Tools($i) != ""} {
      lappend cmd "-[string index $i 0]$Tools($i)"
    }
  }
  if {$Tools(basename) != ""} {
    lappend cmd "-B$Tools(basename)"
  }
  lappend cmd $Tools(cgnsinput)
  if {$Tools(cgnsoutput) != ""} {
    lappend cmd $Tools(cgnsoutput)
  }

  tools_run "Subset Extraction" $cmd $Tools(cgnsoutput)
}

proc check_subset {w} {
  global Tools
  foreach i {i j k} {
    if {$Tools(subset$i) != ""} {
      return [tools_check_input $w]
    }
  }
  errormsg "no subset directions selected" $w
  return 0
}

