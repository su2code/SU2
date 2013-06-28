
array set Import {
  patfile ""
  distloads 0
}

proc patran_import {w name exe} {
  global ProgData Font Import

  set cmd [get_executable $exe 1]
  if {$cmd == ""} return
  set Import(cgnsfile) $ProgData(file,name)

  toplevel $w
  wm title $w "Patran Neutral File Import"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set Import(done) 0}

  import_input $w patfile Patran {.pat .neu .txt}
  import_output $w 1

  FrameCreate $w.opts -text "Options" -font $Font(bold)
  pack $w.opts -side top -padx 5 -pady 2 -fill x
  set opts [FrameGet $w.opts]

  checkbutton $opts.dl -text "Include Distributed Loads" \
    -variable Import(distloads) -onvalue 1 -offvalue 0
  pack $opts.dl -side left

  frame $opts.dup
  pack $opts.dup -side right
  import_dup_check $opts.dup

  if {[import_buttons $w patran_check]} {
    if {$Import(distloads)} {
      lappend cmd -l
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
    lappend cmd $Import(patfile) $Import(cgnsfile)
    import_run "Patran Import" $cmd $Import(cgnsfile)
  }
}

proc patran_check {w} {
  global Import
  if {[string trim $Import(patfile)] == "" ||
      [string trim $Import(cgnsfile)] == ""} {
    errormsg "must specify a Patran and a CGNS file" $w
    return
  }
  if {![file exists $Import(patfile)]} {
    errormsg "Patran Neutral file doesn't exist" $w
    return
  }
  set Import(done) 1
}

