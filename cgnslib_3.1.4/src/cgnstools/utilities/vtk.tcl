# VTK.tcl - VTK import/export

array set Export {
  vtkdir ""
  elemset 0
}

proc vtk_outdir {w} {
  global Export
  set dir [GetDirectory "VTK Output Directory" $Export(vtkdir) $w]
  if {$dir != ""} {
    set Export(vtkdir) $dir
  }
}

proc vtk_export {w name exe} {
  global ProgData Font Export
  set cmd [get_executable $exe 1]
  if {$cmd == ""} return

  toplevel $w
  wm title $w "VTK Export"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Export(done) 0}

  export_input $w 1 1 1

  FrameCreate $w.output -text "VTK Output Directory" -font $Font(bold)
  pack $w.output -side top -padx 5 -pady 2 -fill x
  set f [FrameGet $w.output]

  label $f.lab -text Directory -anchor w
  pack $f.lab -side left
  entry $f.ent -textvariable Export(vtkdir) -width 30
  pack $f.ent -side left -fill x -expand 1
  button $f.but -text Browse -pady 0 -command "vtk_outdir $w"
  pack $f.but -side right -fill y

  FrameCreate $w.opts -text "Options" -font $Font(bold)
  pack $w.opts -side top -pady 2 -fill x
  set opts [FrameGet $w.opts]

  checkbutton $opts.fmt -text "ASCII VTK Output" \
    -variable Export(ascii) -onvalue 1 -offvalue 0 \
    -command tecplot_extension
  pack $opts.fmt -side top -anchor w

  checkbutton $opts.sol -text "Write Element Sets" \
    -variable Export(elemset) -onvalue 1 -offvalue 0
  pack $opts.sol -side top -anchor w

  set Export(cgnsfile) $ProgData(file,name)
  set dir $ProgData(file,dir)
  if {$dir == "" || $dir == "."} {
    set dir [pwd]
  }
  set prefix [file rootname $ProgData(file,name)]
  set Export(vtkdir) "${prefix}_vtk"

  if {[export_buttons $w vtk_export_check]} {
  update idletasks
    if {$Export(basenum) != ""} {
      lappend cmd -b$Export(basenum)
    }
    if {$Export(zonenum) != ""} {
      lappend cmd -z$Export(basenum)
    }
    if {$Export(solution)} {
      if {$Export(solnum) != ""} {
        lappend cmd -s$Export(solnum)
      }
    } else {
      lappend cmd -s0
    }
    if {$Export(elemset)} {
      lappend cmd -e
    }
    lappend cmd $Export(cgnsfile)
    if {$Export(vtkdir) != ""} {
      lappend cmd $Export(vtkdir)
    }
    update
    run_command "VTK Export" $cmd
  }
}

proc vtk_export_check {w} {
  global Export
  if {[string trim $Export(cgnsfile)] == ""} {
    errormsg "must specify a CGNS file" $w
    return
  }
  if {![file exists $Export(cgnsfile)]} {
    errormsg "CGNS input file doesn't exist" $w
    return
  }
  set dir [string trim $Export(vtkdir)]
  if {$dir != "" && [file exists $dir] && ![file isdirectory $dir]} {
    errormsg "Output directory $dir exists and is not a directory" $w
    return
  }
  set Export(done) 1
}

