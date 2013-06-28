# tools.tcl - CGNS tool routines

proc update_cgns {} {
  global ProgData Tools

  if {$ProgData(cgnsversion) == "" || $ProgData(file,name) == ""} return
  set Tools(cgnsinput) $ProgData(file,name)
  if {$Tools(cgnsoutput) == ""} {
    set Tools(cgnsoutput) $ProgData(file,name)
  }
  set Tools(version) ""

  set w .cgnsvers
  toplevel $w
  wm title $w "Change CGNS Version"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Tools(done) 0}

  FrameCreate $w.top
  pack $w.top -side top -padx 5 -pady 5 -fill x
  set top [FrameGet $w.top]

  label $top.ver -text "Change $ProgData(file,name) from CGNS Version\
$ProgData(filevers) to"
  pack $top.ver -side top -anchor w

  set f [frame $top.newvers]
  pack $f -side top -fill x -expand 1
  set n 1
  foreach v {1.2 2.0 2.1 2.2 2.3 2.4 2.5 3.0 3.1} {
    if {[expr $v - $ProgData(libvers) < 0.01]} {
      radiobutton $f.f$n -text $v -variable Tools(version) -value $v
      pack $f.f$n -side left -fill x -expand 1
      if {[expr abs($ProgData(filevers) - $v) < 0.01]} {
        $f.f$n configure -state disabled
      }
      incr n
    }
  }

  set f [frame $top.output]
  pack $f -side top -fill x -expand 1
  label $f.lab -text "CGNS Output"
  pack $f.lab -side left
  entry $f.ent -textvariable Tools(cgnsoutput) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 \
    -command "tools_browse $f.but cgnsoutput"
  pack $f.but -side right -fill y

  if {[tools_interact $w] && $Tools(version) != ""} {
    set cmd [tools_unix_path $ProgData(cgnsversion)]
    lappend cmd $Tools(version) $Tools(cgnsinput)
    if {$Tools(cgnsoutput) != ""} {
      lappend cmd $Tools(cgnsoutput)
    }
    tools_run "Change CGNS Version" $cmd $Tools(cgnsoutput)
    file_stats
  }
}

proc check_cgns {} {
  global ProgData Tools

  if {$ProgData(cgnscheck) == "" || $ProgData(file,name) == ""} return

  set w .cgnscheck
  toplevel $w
  wm title $w "Check CGNS File"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Tools(done) 0}

  FrameCreate $w.top
  pack $w.top -side top -padx 5 -pady 5 -fill x
  set top [FrameGet $w.top]

  set f [frame $top.opts]
  pack $f -side top -fill x -expand 1
  checkbutton $f.verbose -text "Verbose Output" -variable Tools(verbose) \
    -onvalue "-v" -offvalue ""
  checkbutton $f.errors -text "Show Errors" -variable Tools(errors) \
    -onvalue "" -offvalue "-e"
  pack $f.verbose $f.errors -side left -fill x -expand 1

  set f [frame $top.warn]
  pack $f -side top -fill x -expand 1
  label $f.lab -text "Warning Level"
  pack $f.lab -side left
  foreach i {0 1 2 3} {
    radiobutton $f.f$i -text $i -variable Tools(warnings) -value $i
    pack $f.f$i -side left -fill x -expand 1
  }

  if {[tools_interact $w ""]} {
    set cmd [tools_unix_path $ProgData(cgnscheck)]
    foreach i {verbose errors} {
      if {$Tools($i) != ""} {lappend cmd $Tools($i)}
    }
    lappend cmd -w$Tools(warnings) $ProgData(file,name)
    run_command "Check CGNS" $cmd 25
  }
}

proc plot_cgns {} {
  global ProgData tcl_platform
  if {$ProgData(cgnsplot) == "" || $ProgData(file,name) == ""} return
  set cmd [tools_unix_path $ProgData(cgnsplot)]
  lappend cmd $ProgData(file,name)
  if {[catch {eval exec $cmd &} msg]} {
    errormsg $msg
    return
  }
}

proc calc_cgns {} {
  global ProgData
  if {$ProgData(cgnscalc) == "" || $ProgData(file,name) == ""} return
  set cmd [tools_unix_path $ProgData(cgnscalc)]
  lappend cmd $ProgData(file,name)
  if {[catch {eval exec $cmd &} msg]} {
    errormsg $msg
    return
  }
}

proc tools_default {w name exe} {
  global ProgData Tools

  tools_window $w $name $exe run

  if {[tools_interact $w tools_check_exe]} {
    set cmd [tools_unix_path $Tools(exefile)]
    foreach opt $Tools(options) {
      lappend cmd $opt
    }
    if {$Tools(background)} {
      if {[catch {eval exec $cmd &} msg]} {
        errormsg $msg
      }
    } else {
      run_command $name $cmd
    }
  }
}

proc tools_utility {w name exe} {
  global ProgData Tools

  set Tools(cgnsinput) $ProgData(file,name)

  tools_window $w $name $exe utility

  if {[tools_interact $w tools_check]} {
    set cmd [tools_unix_path $Tools(exefile)]
    foreach opt $Tools(options) {
      lappend cmd $opt
    }
    lappend cmd $Tools(cgnsinput)
    if {$Tools(background)} {
      if {[catch {eval exec $cmd &} msg]} {
        errormsg $msg
      }
    } else {
        $name $cmd
    }
  }
}

proc tools_convert {w name exe} {
  global ProgData Tools

  set Tools(cgnsinput) $ProgData(file,name)
  if {$Tools(cgnsoutput) == ""} {
    set Tools(cgnsoutput) $ProgData(file,name)
  }

  tools_window $w $name $exe convert

  if {[tools_interact $w tools_check]} {
    set cmd [tools_unix_path $Tools(exefile)]
    foreach opt $Tools(options) {
      lappend cmd $opt
    }
    lappend cmd $Tools(cgnsinput)
    if {$Tools(cgnsoutput) != ""} {
      lappend cmd $Tools(cgnsoutput)
    }
    tools_run $name $cmd $Tools(cgnsoutput)
  }
}

proc tools_window {w name exe type} {
  global ProgData Tools Font tcl_platform

  if {$exe != ""} {
    set exename [lindex $exe 0]
    set exepath [get_executable $exename]
    if {$exepath != ""} {set exename $exepath}
    if {$tcl_platform(platform) == "windows"} {
      set Tools(exefile) [join [split $exename /] \\]
    } else {
      set Tools(exefile) $exename
    }
    set Tools(options) [lrange $exe 1 end]
  }

  toplevel $w
  wm title $w $name
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Tools(done) 0}

  set labw 12

  set lab "[string toupper [string index $type 0]][string range $type 1 end]"
  FrameCreate $w.top -text "$lab Command" -font $Font(bold)
  pack $w.top -side top -padx 5 -pady 5 -fill x
  set top [FrameGet $w.top]

  set f [frame $top.exe]
  pack $f -side top -fill x -expand 1
  label $f.lab -text Executable -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Tools(exefile) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 -command "tools_browse $w exefile"
  pack $f.but -side right -fill y
  $f.ent xview [$f.ent index end]

  set f [frame $top.opts]
  pack $f -side top -fill x -expand 1
  if {$type == "run"} {
    set lab Arguments
  } else {
    set lab Options
  }
  label $f.lab -text $lab -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Tools(options) -width 30
  pack $f.ent -side left -fill x -expand 1

  if {$type != "run"} {
    set f [frame $top.input]
    pack $f -side top -fill x
    label $f.lab -text "CGNS Input" -width $labw -anchor w
    pack $f.lab -side left
    entry $f.ent -textvariable Tools(cgnsinput) -width 30
    pack $f.ent -side left -fill x -expand 1
    button $f.but -text Browse -pady 0 -command "tools_browse $f.but cgnsinput"
    pack $f.but -side right -fill y
    if {$type == "convert"} {
      set f [frame $top.output]
      pack $f -side top -fill x
      label $f.lab -text "CGNS Output" -width $labw -anchor w
      pack $f.lab -side left
      entry $f.ent -textvariable Tools(cgnsoutput) -width 30
      pack $f.ent -side left -fill x -expand 1
      button $f.but -text Browse -pady 0 \
        -command "tools_browse $f.but cgnsoutput"
      pack $f.but -side right -fill y
    }
  }

  if {$type != "convert"} {
    checkbutton $top.bg -text "Run Command in Background" \
      -variable Tools(background) -onvalue 1 -offvalue 0
    pack $top.bg -side top -anchor w
  }
}

proc tools_check_input {w} {
  global Tools
  if {$Tools(cgnsinput) == "" || ![file exists $Tools(cgnsinput)]} {
    errormsg "CGNS input file not specified or does not exist" $w
    return 0
  }
  return 1
}

proc tools_check_exe {w} {
  global Tools
  if {$Tools(exefile) == ""} {
    errormsg "must specify an executable" $w
    return 0
  }
  if {![file exists $Tools(exefile)] ||
      ![file executable $Tools(exefile)]} {
    errormsg "the utility executable does not exist or is not executable" $w
    return 0
  }
  return 1
}

proc tools_check {w} {
  if {[tools_check_exe $w] &&
      [tools_check_input $w]} {return 1}
  return 0
}

proc tools_unix_path {exepath} {
  global tcl_platform
  if {$tcl_platform(platform) == "windows"} {
    if {![catch {file attributes $exepath -shortname} name]} {
      return [join [split $name \\] /]
    }
    return [join [split $exepath \\] /]
  }
  return $exepath
}

proc tools_cgnsinput {w {loc 0} {labw 13}} {
  global ProgData Font Tools

  FrameCreate $w.cgnsinput -text "CGNS Input" -font $Font(bold)
  pack $w.cgnsinput -side top -padx 5 -pady 2 -fill x
  set wf [FrameGet $w.cgnsinput]

  set f [frame $wf.file]
  pack $f -side top -fill x
  label $f.lab -text Filename -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Tools(cgnsinput) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 -command "tools_browse $f.but cgnsinput"
  pack $f.but -side right -fill y

  if {!$loc} return

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

  set f [frame $wf.zonenum]
  pack $f -side top -fill x
  label $f.lab -text "Zone Index" -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Tools(zonenum) -width 10
  pack $f.ent -side left
  checkbutton $f.but -text "all zones" \
    -variable Tools(zonenum,flag) -onvalue 1 -offvalue 0 \
    -command "tools_state zonenum $f.ent"
  pack $f.but -side left

  tools_state zonenum $f.ent

  if {$loc < 2} return

  set f [frame $wf.solnum]
  pack $f -side top -fill x
  label $f.lab -text "Solution Index" -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Tools(solnum) -width 10
  pack $f.ent -side left
  checkbutton $f.but -text "all solutions" \
    -variable Tools(solnum,flag) -onvalue 1 -offvalue 0 \
    -command "tools_state solnum $f.ent"
  pack $f.but -side left

  tools_state solnum $f.ent
}

proc tools_cgnsoutput {w {newloc ""} {labw 13}} {
  global Font Tools
  FrameCreate $w.cgnsoutput -text "CGNS Output" -font $Font(bold)
  pack $w.cgnsoutput -side top -padx 5 -pady 2 -fill x
  set wf [FrameGet $w.cgnsoutput]

  set f [frame $wf.file]
  pack $f -side top -fill x
  label $f.lab -text Filename -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Tools(cgnsoutput) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 -command "tools_browse $f.but cgnsoutput"
  pack $f.but -side right -fill y

  if {$newloc == "solname"} {
    set lab "Solution Name"
  } elseif {$newloc == "basename"} {
    set lab "Base Name"
  } else {
    return
  }
  set f [frame $wf.newloc]
  pack $f -side top -fill x
  label $f.lab -text $lab -width $labw -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Tools($newloc) -width 30
  pack $f.ent -side left -fill x -expand 1
  checkbutton $f.but -text default \
    -variable Tools($newloc,flag) -onvalue 1 -offvalue 0 \
    -command "tools_state $newloc $f.ent {}"
  pack $f.but -side left

  tools_state $newloc $f.ent ""
}

proc tools_averaging {w} {
  global Font Tools
  FrameCreate $w.avg -text Averaging -font $Font(bold)
  pack $w.avg -side top -padx 5 -pady 2 -fill x
  set avg [FrameGet $w.avg]

  set f [frame $avg.f]
  pack $f -side top -fill x -expand 1
  radiobutton $f.a -text "Simple Averaging" \
    -variable Tools(weight) -value ""
  radiobutton $f.w -text "Volume Weighting" \
    -variable Tools(weight) -value w
  pack $f.a $f.w -side left -expand 1
  return $avg
}

proc tools_options {w {bgopt 0} {labw ""}} {
  global Tools Font

  FrameCreate $w.opts -text "Command Line Options" -font $Font(bold)
  pack $w.opts -side top -padx 5 -pady 2 -fill x
  set opts [FrameGet $w.opts]

  set f [frame $opts.f]
  pack $f -side top -fill x -expand 1
  label $f.lab -text Options
  if {$labw != ""} {
    $f.lab configure -width $labw -anchor w
  }
  pack $f.lab -side left
  entry $f.ent -textvariable Tools(options) -width 30
  pack $f.ent -side left -fill x -expand 1

  if {$bgopt} {
    checkbutton $opts.bg -text "Run Command in Background" \
      -variable Tools(background) -onvalue 1 -offvalue 0
    pack $opts.bg -side top -anchor w
  }
  return $opts
}

proc tools_browse {w ftype {name ""}} {
  global Tools tcl_platform
  if {$ftype == "exefile"} {
    if {$tcl_platform(platform) == "windows"} {
      set types {{{Executable Files} {.exe .cmd .bat .com}} {{All Files} {*}}}
    } else {
      set types {{{All Files} {*}}}
    }
    set fname [FileOpen Executable $Tools(exefile) $w $types]
    if {$fname != ""} {
      if {$tcl_platform(platform) == "windows"} {
        set Tools(exefile) [join [split $fname /] \\]
      } else {
        set Tools(exefile) $fname
      }
    }
  } elseif {$ftype == "cgnsoutput"} {
    set fname [FileSave "CGNS Output File" $Tools(cgnsoutput) . \
      {{{CGNS Files} {.cgns .cgn}} {{All Files} {*}}} cgns]
    if {$fname != ""} {
      if {$tcl_platform(platform) == "windows"} {
        set Tools(cgnsoutput) [join [split $fname /] \\]
      } else {
        set Tools(cgnsoutput) $fname
      }
    }
  } else {
    if {$name == ""} {set name Input}
    set fname [FileOpen "CGNS $name File" $Tools($ftype) . \
      {{{CGNS Files} {.cgns .cga .cgh .cgx}} {{All Files} {*}}}]
    if {$fname != ""} {
      if {$tcl_platform(platform) == "windows"} {
        set Tools($ftype) [join [split $fname /] \\]
      } else {
        set Tools($ftype) $fname
      }
    }
  }
}

proc tools_state {what ent {def 1}} {
  global Tools
  if {$Tools($what,flag)} {
    set Tools($what) ""
    $ent configure -state disabled -cursor {}
  } else {
    if {$Tools($what) == ""} {
      set Tools($what) $def
    }
    $ent configure -state normal -cursor xterm
  }
}

proc tools_interact {w {check_proc tools_check_input}} {
  global Tools

  set f [frame $w.but]
  pack $f -side bottom -pady 5
  button $f.accept -text Accept -width 6 -default active
  button $f.cancel -text Cancel -width 6 -command {set Tools(done) 0}
  pack $f.accept $f.cancel -side left -padx 5

  if {$check_proc == ""} {
    $f.accept configure -command {set Tools(done) 1}
    bind $w <Return> {set Tools(done) 1}
  } else {
    $f.accept configure -command "
      if {!\[catch {$check_proc $w} result\] && \$result} {
        set Tools(done) 1
      }
    "
    bind $w <Return> "
      if {!\[catch {$check_proc $w} result\] && \$result} {
        set Tools(done) 1
      }
    "
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
  tkwait variable Tools(done)
  catch {focus $oldFocus}
  destroy $w
  if {$oldGrab != ""} {
    if {$grabStatus == "global"} {
      grab -global $oldGrab
    } else {
      grab $oldGrab
    }
  }
  return $Tools(done)
}

proc tools_run {title cmd {outfile ""}} {
  global ProgData
  if {$outfile == ""} {
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

