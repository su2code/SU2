
array set HelpData {
  browser ""
  htmlopts ""
  filelist ""
  winhtml 0
  chmfile ""
  cgns,http "http://www.grc.nasa.gov/www/cgns/CGNS_docs_current/index.html"
}

proc get_browser {} {
  global tcl_platform env

  if {$tcl_platform(platform) != "windows"} {
    foreach name {firefox netscape mozilla} {
      set browser [find_file executable $name]
      if {$browser != ""} {return $browser}
    }
    return ""
  }

  # get browser from registry

  if {![catch {registry get HKEY_CLASSES_ROOT\\.htm {}} htmlfile] ||
      ![catch {registry get HKEY_CLASSES_ROOT\\.html {}} htmlfile]} {
    if {![catch {registry get \
      HKEY_CLASSES_ROOT\\$htmlfile\\shell\\open\\command {}} browser]} {
      return [join [split $browser \\] /]
    }
  }

  # use assoc and ftype commands

  if {[info exists env(COMSPEC)]} {
    set comspec $env(COMSPEC)
  } elseif {[string match "*9?" $tcl_platform(os)]} {
    set comspec command
  } else {
    set comspec cmd
  }
  if {[catch {exec $comspec /c assoc .html} assoc]} {return ""}
  set ftype [string trim [lindex [split $assoc =] 1]]
  if {$ftype == ""} {return ""}
  if {[catch {exec $comspec /c ftype $ftype} cmdstr]} {return ""}
  set cmd [string trim [lindex [split $cmdstr =] 1]]
  return [join [split $cmd \\] /]
}

proc help_find {name} {
  global cmd_dir tcl_platform
  set hlpfile [find_file exists $name $cmd_dir $cmd_dir/help]
  if {$hlpfile == "" && [file isdirectory $cmd_dir/help]} {
    foreach i [glob $cmd_dir/help/*] {
      if {[file isdirectory $i]} {
        set hlpfile [find_file exists $name $i]
        if {$hlpfile != ""} break
      }
    }
  }
  if {$hlpfile != ""} {
    if {$tcl_platform(platform) == "windows"} {
      set hlpfile [file attributes $hlpfile -shortname]
    }
    if {[string tolower [file extension $hlpfile]] != ".chm"} {
      set hlpfile "file://$hlpfile"
    }
  }
  return $hlpfile
}

proc help_defaults {} {
  global HelpData HelpNew cmd_dir tcl_platform
  array set HelpNew {
    browser ""
    htmlopts ""
  }
  set browser [get_browser]
  if {$browser != "" && $tcl_platform(platform) == "windows"} {
    if {[catch {file attributes $browser -shortname} name]} {
      if {[catch {file attributes [lindex $browser 0] -shortname} name]} {
        set name $browser
      } else {
        for {set n 1} {$n < [llength $browser]} {incr n} {
          set opt [lindex $browser $n]
          if {$opt != "%1"} {
            lappend HelpNew(htmlopts) $opt
          }
        }
      }
    }
    set browser $name
  }
  set HelpNew(browser) $browser

  foreach i $HelpData(filelist) {
    if {[info exists HelpData($i,http)]} {
      set HelpNew($i) $HelpData($i,http)
    } else {
      set HelpNew($i) [help_find $HelpData($i,file).html]
      if {$HelpData(winhtml)} {
        set chmfile [help_find $HelpData($i,file).chm]
        if {$chmfile != ""} {
          set HelpNew($i) $chmfile
        }
      }
    }
  }
}

proc help_init {args} {
  global HelpData HelpNew
  set HelpData(filelist) ""
  foreach a $args {
    set i [lindex $a 0]
    set HelpData($i,label) [lindex $a 1]
    if {$HelpData($i,label) == ""} {set HelpData($i,label) $i}
    set HelpData($i,file) [lindex $a 2]
    if {$HelpData($i,file) == ""} {set HelpData($i,file) $i}
    lappend HelpData(filelist) $i
  }
  set data [list browser htmlopts]
  eval lappend data $HelpData(filelist)
  help_defaults
  foreach i $data {
    set HelpData($i) $HelpNew($i)
  }
  if {[info procs tclreg_get] != ""} {
    foreach i $data {
      if {![catch {tclreg_get Help $i} val]} {
        set HelpData($i) $val
      }
    }
  }
  if {[info commands WinHtml] != ""} {
    set HelpData(winhtml) 1
  }
  catch help_menu
}

proc help_setup {} {
  global HelpData HelpNew tcl_platform Font

  set w .hlpsetup
  catch {destroy $w}
  toplevel $w
  wm title $w "Help Setup"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW "destroy $w"

  set lw 7
  set ew 50

  foreach i {browser htmlopts} {
    set HelpNew($i) $HelpData($i)
  }
  foreach i $HelpData(filelist) {
    set HelpNew($i) $HelpData($i)
    set len [string length $HelpData($i,label)]
    if {$len > $lw} {set lw $len}
  }

  FrameCreate $w.hb -text "HTML Browser" -font $Font(bold)
  pack $w.hb -side top -padx 5 -pady 5
  set hb [FrameGet $w.hb]

  foreach j {{browser Path} {htmlopts Options}} {
    set i [lindex $j 0]
    set f [frame $hb.$i]
    pack $f -side top -fill x -expand 1
    label $f.lab -text [lindex $j 1] -width $lw -anchor w
    pack $f.lab -side left
    entry $f.ent -width $ew -textvariable HelpNew($i) -highlightthickness 0
    pack $f.ent -side left -fill x -expand 1 -padx 2
    $f.ent xview end
    if {$i != "htmlopts"} {
      button $f.but -text Browse -padx 0 -pady 0 -command "help_browse $i"
      pack $f.but -side right
    }
  }

  FrameCreate $w.hf -text "Documentation URL" -font $Font(bold)
  pack $w.hf -side top -padx 5
  set hf [FrameGet $w.hf]

  foreach i $HelpData(filelist) {
    set f [frame $hf.$i]
    pack $f -side top -fill x -expand 1
    label $f.lab -text $HelpData($i,label) -width $lw -anchor w
    pack $f.lab -side left
    entry $f.ent -width $ew -textvariable HelpNew($i) -highlightthickness 0
    pack $f.ent -side left -fill x -expand 1 -padx 2
    $f.ent xview end
    button $f.but -text Browse -padx 0 -pady 0 -command "help_browse $i"
    pack $f.but -side right
  }

  set b [frame $w.but]
  pack $b -side top -pady 5
  button $b.accept -text Accept -width 8 -default active \
    -command "help_check $w"
  button $b.default -text Defaults -width 8 -command help_defaults
  button $b.cancel -text Cancel -width 8 -command "destroy $w"
  pack $b.accept $b.default $b.cancel -side left -padx 5

  bind $w <Return> "help_check $w"

  center_window $w .
  set oldFocus [focus]
  set oldGrab [grab current $w]
  if {$oldGrab != ""} {
    set grabStatus [grab status $oldGrab]
  }
  catch {grab $w}
  tkwait visibility $w
  focus $w
  tkwait window $w
  catch {focus $oldFocus}
  if {$oldGrab != ""} {
    if {$grabStatus == "global"} {
      grab -global $oldGrab
    } else {
      grab $oldGrab
    }
  }
  catch help_menu
}

proc help_browse {what} {
  global HelpData HelpNew tcl_platform
  if {$what == "browser"} {
    if {$tcl_platform(platform) == "windows"} {
      set filelist {{{Executable Files} {.exe .com .bat}}}
    } else {
      set filelist {}
    }
    lappend filelist {{All Files} {*}}
    set fname [FileOpen "Select HTML Browser" $HelpNew(browser) . $filelist]
  } else {
    set filelist {{{HTML Files} {.html .htm}}}
    if {$HelpData(winhtml)} {
      lappend filelist {{CHM Files} .chm}
    }
    if [string match "file://*" $HelpNew($what)] {
      set oldname [string range $HelpNew($what) 7 end]
    } else {
      set oldname $HelpNew($what)
    }
    set fname [FileOpen "$HelpData($what,label) Documentation" \
      $oldname . $filelist]
  }
  if {$fname != ""} {
    if {$tcl_platform(platform) == "windows" &&
      ![catch {file attributes $fname -shortname} name]} {
      set fname $name
    }
    if {[string tolower [file extension $fname]] == ".chm"} {
      set HelpNew($what) $fname
    } else {
      set HelpNew($what) file://$fname
    }
  }
}

proc help_check {w} {
  global HelpData HelpNew
  set browser [find_file executable $HelpNew(browser)]
  if {$browser == ""} {
    errormsg "can't find HTML browser or it's not executable" $w
    return
  }
  set data [list browser htmlopts]
  eval lappend data $HelpData(filelist)
  foreach i $data {
    set HelpData($i) $HelpNew($i)
  }
  if {[info procs tclreg_set] != ""} {
    foreach i $data {
      catch {tclreg_set Help $i $HelpData($i)}
    }
  }
  destroy $w
}

proc help_valid {what} {
  global HelpData
  if {![info exists HelpData($what)] || $HelpData($what) == ""} {
    return 0
  }
  if {$HelpData(winhtml) &&
      [string tolower [file extension $HelpData($what)]] == ".chm"} {
    return 1
  }
  if {$HelpData(browser) != "" &&
      [file executable $HelpData(browser)]} {
    return 1
  }
  return 0
}

proc help_show {topic {tag ""} {chmhtml ""}} {
  global HelpData tcl_platform
  set htmlfile $HelpData($topic)
  if {$htmlfile == ""} return
  if {$HelpData(winhtml) &&
      [string tolower [file extension $htmlfile]] == ".chm"} {
    if {$HelpData(chmfile) != $htmlfile} {
      if {[catch {WinHtml file $htmlfile} msg]} {
        errormsg $msg
        return
      }
      set HelpData(chmfile) $htmlfile
    }
    if {$chmhtml == ""} {
      if {[catch {WinHtml index} msg]} {
        errormsg $msg
      }
    } else {
      if {[catch {eval WinHtml topic $chmhtml $tag} msg]} {
        errormsg $msg
      }
    }
    return
  }
  if {$HelpData(browser) == "" || $htmlfile == ""} {
    errormsg "browser and/or HTML file not set up"
  } else {
    set cmd "$HelpData(browser) $HelpData(htmlopts) "
    if {$tcl_platform(platform) == "windows"} {
      append cmd "\"$htmlfile\""
    } else {
      append cmd $htmlfile
    }
    if {$tag != ""} {
      append cmd "\#$tag"
    }
    if {[catch {eval exec $cmd &} msg]} {
      errormsg $msg
    }
  }
}

