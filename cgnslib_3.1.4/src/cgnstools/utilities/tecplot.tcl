# tecplot.tcl - Tecplot import/export

array set Import {
  tecfile ""
  fixbricks 0
  singlezone 0
}

set Export(tecfile) ""

proc tecplot_import {w name exe} {
  global ProgData Font Import

  set cmd [get_executable $exe 1]
  if {$cmd == ""} return
  set Import(cgnsfile) $ProgData(file,name)

  toplevel $w
  wm title $w "Tecplot File Import"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Import(done) 0}

  import_input $w tecfile Tecplot {.dat .plt .tec}
  import_output $w 1

  FrameCreate $w.opts -text "Options" -font $Font(bold)
  pack $w.opts -side top -padx 5 -pady 2 -fill x
  set opts [FrameGet $w.opts]

  set f [frame $opts.f]
  pack $f -side left

  checkbutton $f.fb -text "Fix Degenerate Bricks" \
    -variable Import(fixbricks) -onvalue 1 -offvalue 0
  checkbutton $f.sz -text "Import as Single Zone" \
    -variable Import(singlezone) -onvalue 1 -offvalue 0 -state disabled
  pack $f.fb $f.sz -side top -anchor w

  frame $opts.dup
  pack $opts.dup -side right
  import_dup_check $opts.dup

  if {[import_buttons $w tecplot_import_check]} {
    if {$Import(fixbricks)} {
      lappend cmd -f
    }
    if {$Import(dupcheck,flag)} {
      lappend cmd -$Import(dupcheck)
      if {$Import(duptol) != ""} {
        lappend cmd -t$Import(duptol)
      }
    }
    if {$Import(basename) != ""} {
      lappend cmd -B$Import(basename)
    }
    lappend cmd $Import(tecfile) $Import(cgnsfile)
    import_run "Tecplot Import" $cmd $Import(cgnsfile)
  }
}

proc tecplot_import_check {w} {
  global Import
  if {[string trim $Import(tecfile)] == "" ||
      [string trim $Import(cgnsfile)] == ""} {
    errormsg "must specify a Tecplot and a CGNS file" $w
    return
  }
  if {![file exists $Import(tecfile)]} {
    errormsg "Tecplot input file doesn't exist" $w
    return
  }
  set Import(done) 1
}

proc tecplot_export {w name exe} {
  global ProgData Font Export
  set cmd [get_executable $exe 1]
  if {$cmd == ""} return

  toplevel $w
  wm title $w "Tecplot Export"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Export(done) 0}

  export_input $w 1 0 1
  export_output $w tecfile Tecplot {.dat .plt .tec}

  FrameCreate $w.opts -text "Options" -font $Font(bold)
  pack $w.opts -side top -pady 2 -fill x
  set opts [FrameGet $w.opts]

  checkbutton $opts.fmt -text "ASCII Tecplot Output" \
    -variable Export(ascii) -onvalue 1 -offvalue 0 \
    -command tecplot_extension
  pack $opts.fmt -side top -anchor w

  checkbutton $opts.sol -text "Grid Ouput Only (no solution)" \
    -variable Export(solution) -onvalue 0 -offvalue 1
  pack $opts.sol -side top -anchor w

  checkbutton $opts.weight -text "Use Volume Weighting" \
    -variable Export(weight) -onvalue 1 -offvalue 0
  pack $opts.weight -side top -anchor w

  set Export(cgnsfile) $ProgData(file,name)
  set Export(tecfile) [file rootname $ProgData(file,name)]
  tecplot_extension

  if {[export_buttons $w tecplot_export_check]} {
  update idletasks
    if {$Export(basenum) != ""} {
      lappend cmd -b$Export(basenum)
    }
    if {$Export(solution)} {
      if {$Export(solnum) != ""} {
        lappend cmd -S$Export(solnum)
      }
    } else {
      lappend cmd -S0
    }
    if {$Export(ascii)} {
      lappend cmd -a
    }
    if {$Export(weight)} {
      lappend cmd -w
    }
    lappend cmd $Export(cgnsfile) $Export(tecfile)
    update
    run_command "Tecplot Export" $cmd
  }
}

proc tecplot_extension {} {
  global Export
  if {$Export(tecfile) == ""} return
  if {$Export(ascii)} {
    set ext dat
  } else {
    set ext plt
  }
  set Export(tecfile) [file rootname $Export(tecfile)].$ext
}

proc tecplot_export_check {w} {
  global Export
  if {[string trim $Export(cgnsfile)] == "" ||
      [string trim $Export(tecfile)] == ""} {
    errormsg "must specify a CGNS and a Tecplot file" $w
    return
  }
  if {![file exists $Export(cgnsfile)]} {
    errormsg "CGNS input file doesn't exist" $w
    return
  }
  set Export(done) 1
}

