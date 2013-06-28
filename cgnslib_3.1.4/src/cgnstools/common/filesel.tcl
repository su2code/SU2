
if {![info exists UseNativeDialogs] || $UseNativeDialogs == ""} {
  if {$tcl_platform(platform) == "windows"} {
    set UseNativeDialogs 1
  } else {
    set UseNativeDialogs 0
  }
}

#----- file selection

set fileselect(toplevel) .fileselect
set fileselect(button)   ""
set fileselect(filter)   *
set fileselect(filename) ""
set fileselect(xpos) -1
set fileselect(ypos) -1
set fileselect(down) [image create bitmap -data "
#define down_width 17
#define down_height 10
static unsigned char down_bits[] = {
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf8, 0x3f, 0x00, 0xf0, 0x1f, 0x00,
   0xe0, 0x0f, 0x00, 0xc0, 0x07, 0x00, 0x80, 0x03, 0x00, 0x00, 0x01, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00};"]

proc fileselect {title defname {x -1} {y -1} {filter *}} {
  global fileselect
  set curdir [pwd]

  set top $fileselect(toplevel)
  catch {destroy $top}
  toplevel $top
  wm title $top $title
  wm transient $top [winfo toplevel [winfo parent $top]]

  #--- buttons

  set w [frame $top.button -relief sunken -bd 2]
  pack $w -side bottom -fill x -padx 2 -pady 2 -anchor s -expand 1
  button $w.accept -relief raised -text Accept \
    -command "set fileselect(button) 1"
  button $w.cancel -relief raised -text Cancel \
    -command "set fileselect(button) 0"
  pack $w.accept $w.cancel -side left -expand 1 -fill x -padx 5 -pady 2

  $fileselect(down) configure -foreground [$w.accept cget -foreground]

  #--- filter entry

  set w [frame $top.filter]
  pack $w -side top -fill x -padx 2 -pady 2 -anchor n
  label $w.label -text "Filter:"
  pack $w.label -side left -anchor w
  entry $w.entry -relief sunken
  pack $w.entry -side left -fill x -expand 1
  menubutton $w.menubut -relief raised -image $fileselect(down) \
    -menu $w.menubut.m -direction left
  menu $w.menubut.m -tearoff 0
  pack $w.menubut -side right -fill y

  bind $w.entry <Escape> {filesel:entry %W [pwd]/$fileselect(filter)}
  bind $w.entry <Return> {
    global fileselect
    set dir [%W get]
    set fileselect(filter) "*"
    set fileselect(filename) ""
    if {![file isdirectory $dir]} {
      set fileselect(filter) [file tail $dir]
      set dir [file dirname $dir]
    }
    catch {cd $dir}
    filesel:Select
  }
  bind $w.entry <Control-u> {%W delete 0 end}

  #--- file selection entry

  set w [frame $top.filename]
  pack $w -side bottom -fill x -padx 2 -pady 2 -anchor s
  label $w.label -text "Filename:"
  pack $w.label -side left -anchor w
  entry $w.entry -relief sunken
  pack $w.entry -side left -fill x -expand 1

  bind $w.entry <Escape> {filesel:entry %W [pwd]/$fileselect(filename)}
  bind $w.entry <Return> {
    global fileselect
    if [file isdirectory [%W get]] {
      set fileselect(filename) ""
      catch {cd [%W get]}
      filesel:Select
    } else {
      set fileselect(button) 1
    }
  }
  bind $w.entry <Control-u> {%W delete 0 end}

  #--- directory list

  set w [frame $top.dirlist]
  pack $w -side left -fill both -padx 2 -expand 1
  label $w.label -text "Directories:"
  pack $w.label -side top -anchor w
  scrollbar $w.scroll -relief sunken -command "$w.list yview"
  pack $w.scroll -side right -fill y
  listbox $w.list -relief sunken -yscroll "$w.scroll set" \
    -selectmode browse
  pack $w.list -side left -fill both -expand 1

  bind $w.list <Double-Button-1> {filesel:Directory %W}
  bind $w.list <Return> {filesel:Directory %W}

  #--- file list

  set w [frame $top.filelist]
  pack $w -side left -fill both -padx 2 -expand 1
  label $w.label -text "Files:"
  pack $w.label -side top -anchor w
  scrollbar $w.scroll -relief sunken -command "$w.list yview"
  pack $w.scroll -side right -fill y
  listbox $w.list -relief sunken -yscroll "$w.scroll set" \
    -selectmode browse
  pack $w.list -side left -fill both -expand 1

  bind $w.list <ButtonRelease-1> {filesel:Filename %W}
  bind $w.list <Return> {filesel:Filename %W}
  bind $w.list <Double-Button-1> {
    global fileselect
    if [filesel:Filename %W] {set fileselect(button) 1}
  }

  #--- initialize

  set fileselect(filter) ""
  foreach filt $filter {
    if {[llength $filt] == 2} {
      set extlist {}
      foreach ext [lindex $filt 1] {
        if {[string index $ext 0] != "*"} {
          set e "*$ext"
        } else {
          set e $ext
        }
        if {$extlist == ""} {
          set extlist $e
        } else {
          append extlist ",$e"
        }
      }
      $top.filter.menubut.m add command \
        -label "[lindex $filt 0] \($extlist\)" \
        -command "set fileselect(filter) $extlist; filesel:Select"
      if {$fileselect(filter) == ""} {
        set fileselect(filter) $extlist
      }
    }
  }
  if {$fileselect(filter) == ""} {
    $top.filter.menubut configure -state disabled
    set fileselect(filter) $filter
  }
  set fileselect(filename) ""
  if {$defname != ""} {
    if [file isdirectory $defname] {
      catch {cd $defname}
    } else {
      catch {cd [file dirname $defname]}
      if [string match "*\[*?\]*" [file tail $defname]] {
        set fileselect(filter) [file tail $defname]
      } else {
        set fileselect(filename) [file tail $defname]
      }
    }
  }
  filesel:Select 0

  #--- locate the window

  if {$x != "" && [winfo exists $x]} {
    center_window $top $x
  } else {
    wm withdraw $top
    update idletasks
    set ww [winfo reqwidth  $top]
    set wh [winfo reqheight $top]

    if {$x < 0} {
      if {$fileselect(xpos) < 0} {
        if [winfo ismapped .] {
          set x [expr [winfo rootx .] + ([winfo width .] - $ww) / 2]
        } else {
          set x [expr ([winfo screenwidth $top] - $ww) / 2]
        }
      } else {
        set x $fileselect(xpos)
      }
    } else {
      set x [expr $x - $ww / 2]
    }
    if {$x < 0} {
      set pos +0
    } elseif {[expr $x + $ww] > [winfo screenwidth $top]} {
      set pos -0
    } else {
      set pos +$x
    }

    if {$y < 0} {
      if {$fileselect(ypos) < 0} {
        if [winfo ismapped .] {
          set y [expr [winfo rooty .] + ([winfo height .] - $wh) / 2]
        } else {
          set y [expr ([winfo screenheight $top] - $wh) / 2]
        }
      } else {
        set y $fileselect(ypos)
      }
    } else {
      set y [expr $y - $wh / 2]
    }
    if {$y < 0} {
      set pos $pos+0
    } elseif {[expr $y + $wh] > [winfo screenheight $top]} {
      set pos $pos-0
    } else {
      set pos $pos+$y
    }

    wm geom $top $pos
    wm deiconify $top
  }

  #--- wait for button click

  set oldFocus [focus]
  set oldGrab [grab current $w]
  if {$oldGrab != ""} {
    set grabStatus [grab status $oldGrab]
  }
  grab $top
  $top.filename.entry selection range 0 end
  focus $top.filename.entry
  tkwait variable fileselect(button)
  if $fileselect(button) {
    set filename [$top.filename.entry get]
  } else {
    set filename ""
  }
  set fileselect(xpos) [winfo rootx $top]
  set fileselect(ypos) [winfo rooty $top]
  destroy $top
  catch {focus $oldFocus}
  if {$oldGrab != ""} {
    if {$grabStatus == "global"} {
      grab -global $oldGrab
    } else {
      grab $oldGrab
    }
  }
  catch {cd $curdir}
  return $filename
}

proc filesel:Directory {w} {
  global fileselect
  if {[$w curselection] != ""} {
    set dir [$w get [$w curselection]]
    if {$dir == ".. <up>"} {
      set dir ".."
    }
    catch {cd $dir}
    if {$fileselect(filename) != "" &&
      ![file isfile $fileselect(filename)]} {
      set fileselect(filename) ""
    }
    filesel:Select
  }
}

proc filesel:Filename {w} {
  global fileselect
  if {[$w curselection] != ""} {
    set fileselect(filename) [$w get [$w curselection]]
    filesel:entry $fileselect(toplevel).filename.entry \
      [pwd]/$fileselect(filename)
    return 1
  }
  return 0
}

proc filesel:Select {{vis 1}} {
  global fileselect tcl_platform
  set curdir [pwd]
  if {[info exists tcl_platform(platform)] &&
    $tcl_platform(platform) == "windows" &&
    [string index $curdir 1] == ":"} {
    set curdir [string range $curdir 2 end]
  }
  if $vis {
    $fileselect(toplevel) configure -cursor watch
    update idletasks
  }

  #--- fill in directories

  set w $fileselect(toplevel).dirlist.list
  catch {set files [lsort [glob -nocomplain *]]}
  $w delete 0 end
  if {$curdir != "/"} {
    $w insert 0 ".. <up>"
  }
  foreach f $files {
    if [file isdirectory "./$f"] {
      $w insert end $f
    }
  }

  #--- fill in files

  set w $fileselect(toplevel).filelist.list
  catch {set files [lsort [eval glob -nocomplain \
    [split $fileselect(filter) ,]]]}
  $w delete 0 end
  set n 0
  foreach f $files {
    if [file isfile "./$f"] {
      $w insert end $f
      if {$f == $fileselect(filename)} {
          $w selection set $n
      }
      incr n 1
    }
  }

  #--- update filter and filename entries

  if {$curdir == "/"} {
    filesel:entry $fileselect(toplevel).filter.entry \
      [pwd]$fileselect(filter)
    filesel:entry $fileselect(toplevel).filename.entry \
      [pwd]$fileselect(filename)
  } else {
    filesel:entry $fileselect(toplevel).filter.entry \
      [pwd]/$fileselect(filter)
    filesel:entry $fileselect(toplevel).filename.entry \
      [pwd]/$fileselect(filename)
  }
  if $vis {
    $fileselect(toplevel) configure -cursor {}
  }
}

proc filesel:entry {w s} {
  $w delete 0 end
  if {$s != ""} {
    $w insert insert $s
    set c [$w index insert]
    if {($c < [$w index @0]) || ($c > [$w index @[winfo width $w]])} {
	$w xview $c
    }
  }
}

#----- directory selection

proc dirselect {title defname {x -1} {y -1}} {
  global fileselect
  set curdir [pwd]

  set top $fileselect(toplevel)
  catch {destroy $top}
  toplevel $top
  wm title $top $title
  wm transient $top [winfo toplevel [winfo parent $top]]

  #--- filter entry

  set w [frame $top.filter]
  pack $w -side top -fill x -padx 2 -pady 2 -anchor n
  label $w.label -text "Filter:"
  pack $w.label -side left -anchor w
  entry $w.entry -relief sunken
  pack $w.entry -side left -fill x -expand 1

  bind $w.entry <Escape> {filesel:entry %W [pwd]/$fileselect(filter)}
  bind $w.entry <Return> {
    global fileselect
    set dir [%W get]
    set fileselect(filter) "*"
    set fileselect(filename) ""
    if {![file isdirectory $dir]} {
      set fileselect(filter) [file tail $dir]
      set dir [file dirname $dir]
    }
    catch {cd $dir}
    filesel:Select
  }
  bind $w.entry <Control-u> {%W delete 0 end}

  #--- buttons

  set w [frame $top.button -relief sunken -bd 2]
  pack $w -side bottom -fill x -padx 2 -pady 2 -anchor s -expand 1
  button $w.accept -relief raised -text Accept \
    -command "set fileselect(button) 1"
  button $w.cancel -relief raised -text Cancel \
    -command "set fileselect(button) 0"
  pack $w.accept $w.cancel -side left -expand 1 -fill x -padx 5 -pady 2

  #--- file selection entry

  set w [frame $top.filename]
  pack $w -side bottom -fill x -padx 2 -pady 2 -anchor s
  label $w.label -text "Directory:"
  pack $w.label -side left -anchor w
  entry $w.entry -relief sunken
  pack $w.entry -side left -fill x -expand 1

  bind $w.entry <Escape> {filesel:entry %W [pwd]/}
  bind $w.entry <Return> {
    global fileselect
    if [file isdirectory [%W get]] {
      set fileselect(button) 1
    }
  }
  bind $w.entry <Control-u> {%W delete 0 end}

  #--- directory list

  set w [frame $top.dirlist]
  pack $w -side left -fill both -padx 2 -expand 1
  label $w.label -text "Directories:"
  pack $w.label -side top -anchor w
  scrollbar $w.scroll -relief sunken -command "$w.list yview"
  pack $w.scroll -side right -fill y
  listbox $w.list -relief sunken -yscroll "$w.scroll set" \
    -selectmode browse -exportselection 0
  pack $w.list -side left -fill both -expand 1

  bind $w.list <Double-Button-1> {filesel:Directory %W}
  bind $w.list <Return> {filesel:Directory %W}

  #--- file list

  set fg [$top.button.accept cget -disabledforeground]
  set bg [$w.list cget -bg]
  if {[winfo rgb $w.list $fg] == [winfo rgb $w.list $bg]} {
    set fg [$w.list cget -fg]
  }
  set w [frame $top.filelist]
  pack $w -side left -fill both -padx 2 -expand 1
  label $w.label -text "Files:"
  pack $w.label -side top -anchor w
  scrollbar $w.scroll -relief sunken -command "$w.list yview"
  pack $w.scroll -side right -fill y
  listbox $w.list -relief sunken -yscroll "$w.scroll set" \
    -fg $fg -bg $bg -selectforeground $fg -selectbackground $bg \
    -selectborderwidth 0 -takefocus 0 -highlightthickness 0
  pack $w.list -side left -fill both -expand 1

  bind $w.list <Button-1> {}
  bind $w.list <Double-Button-1> {}

  #--- initialize

  set fileselect(filter)   "*"
  set fileselect(filename) ""
  if {$defname != ""} {
    if [file isdirectory $defname] {
      catch {cd $defname}
    } else {
      catch {cd [file dirname $defname]}
      if [string match "*\[*?\]*" [file tail $defname]] {
        set fileselect(filter) [file tail $defname]
      }
    }
  }
  filesel:Select 0

  #--- locate the window

  if {$x != "" && [winfo exists $x]} {
    center_window $top $x
  } else {
    wm withdraw $top
    update idletasks
    set ww [winfo reqwidth  $top]
    set wh [winfo reqheight $top]

    if {$x < 0} {
      if {$fileselect(xpos) < 0} {
        if [winfo ismapped .] {
          set x [expr [winfo rootx .] + ([winfo width .] - $ww) / 2]
        } else {
          set x [expr ([winfo screenwidth $top] - $ww) / 2]
        }
      } else {
        set x $fileselect(xpos)
      }
    } else {
      set x [expr $x - $ww / 2]
    }
    if {$x < 0} {
      set pos +0
    } elseif {[expr $x + $ww] > [winfo screenwidth $top]} {
      set pos -0
    } else {
      set pos +$x
    }

    if {$y < 0} {
      if {$fileselect(ypos) < 0} {
        if [winfo ismapped .] {
          set y [expr [winfo rooty .] + ([winfo height .] - $wh) / 2]
        } else {
          set y [expr ([winfo screenheight $top] - $wh) / 2]
        }
      } else {
        set y $fileselect(ypos)
      }
    } else {
      set y [expr $y - $wh / 2]
    }
    if {$y < 0} {
      set pos $pos+0
    } elseif {[expr $y + $wh] > [winfo screenheight $top]} {
      set pos $pos-0
    } else {
      set pos $pos+$y
    }

    wm geom $top $pos
    wm deiconify $top
  }

  #--- wait for button click

  set oldFocus [focus]
  set oldGrab [grab current $top]
  if {$oldGrab != ""} {
    set grabStatus [grab status $oldGrab]
  }
  grab $top
  $top.filter.entry selection range 0 end
  focus -force $top.filter.entry
  tkwait variable fileselect(button)
  if $fileselect(button) {
    set dirname [$top.filename.entry get]
    if ![catch {cd $dirname}] {
      set dirname [pwd]
    }
  } else {
    set dirname ""
  }
  set fileselect(xpos) [winfo rootx $top]
  set fileselect(ypos) [winfo rooty $top]
  destroy $top
  catch {focus $oldFocus}
  if {$oldGrab != ""} {
    if {$grabStatus == "global"} {
      grab -global $oldGrab
    } else {
      grab $oldGrab
    }
  }
  catch {cd $curdir}
  return $dirname
}

proc FileOpen {title defname {wref .} {types ""}} {
  global UseNativeDialogs
  if {$defname == ""} {
    set dir [pwd]
    set file ""
  } else {
    if [file isdirectory $defname] {
      set dir $defname
      set file ""
    } else {
      set dir [file dirname $defname]
      set file [file tail $defname]
    }
    if {$dir == "."} {
      set dir [pwd]
    }
  }
  if $UseNativeDialogs {
    set opts "-title {$title} -initialdir {$dir}"
    if {$defname != "" && [file exists $defname]} {
      append opts " -initialfile {[file tail $defname]}"
    }
    if {$file != ""} {
      append opts " -initialfile {$file}"
    }
    if {$wref != "" && [winfo exists $wref]} {
      append opts " -parent $wref"
    }
    if {$types != ""} {
      append opts " -filetypes {$types}"
    }
    return [eval tk_getOpenFile $opts]
  }
  if {$types == ""} {
    set types *
  }
  while 1 {
    set fname [fileselect $title $dir/$file $wref {} $types]
    if {$fname == "" || [file exists $fname]} {return $fname}
    dialog .fileopen $wref {} $title "$fname
Cannot find this file. Please verify the path and filename." info 0 Ok
  }
}

proc FileSave {title defname {wref .} {types ""} {defext ""}} {
  global UseNativeDialogs
  if {$defname == ""} {
    set dir [pwd]
    set file ""
  } else {
    if [file isdirectory $defname] {
      set dir $defname
      set file ""
    } else {
      set dir [file dirname $defname]
      set file [file tail $defname]
    }
    if {$dir == "."} {
      set dir [pwd]
    }
  }
  if $UseNativeDialogs {
    set opts "-title {$title} -initialdir {$dir}"
    if {$file != ""} {
      append opts " -initialfile {$file}"
    }
    if {$wref != "" && [winfo exists $wref]} {
      append opts " -parent $wref"
    }
    if {$types != ""} {
      append opts " -filetypes {$types}"
    }
    if {$defext != ""} {
      append opts " -defaultextension $defext"
    }
    return [eval tk_getSaveFile $opts]
  }
  if {$types == ""} {
    set types *
  }
  while 1 {
    set fname [fileselect $title $dir/$file $wref {} $types]
    if {$fname == ""} {return ""}
    if {$defext != "" && [file extension $fname] == ""} {
      append fname ".$defext"
    }
    if {![file exists $fname]} {return $fname}
    if ![dialog .filesave $wref {} $title "$fname
This file already exists. Replace the existing file ?" info 0 Yes No] {
      return $fname
    }
  }
}

proc GetDirectory {title dir {wref .}} {
  global UseNativeDialogs
  if {$dir == "" || $dir == "." || ![file isdirectory $dir]} {
    set dir [pwd]
  }
  if {$UseNativeDialogs} {
    set opts "-title {$title} -initialdir {$dir}"
    if {$wref != "" && [winfo exists $wref]} {
      append opts " -parent $wref"
    }
    return [eval tk_getDirectory $opts]
  }
  return [dirselect $title $dir $wref]
}

proc TempName {{prefix ""} {extlist {}}} {
  if {$prefix == ""} {set prefix temp}
  set name $prefix
  foreach ext $extlist {
    if {[file exists $name.$ext]} {
      set name ""
      break
    }
  }
  if {$name != ""} {return $name}
  for {set n 0} {$n < 1000} {incr n} {
    set name [format "$prefix%3.3d" $n]
    if {[file exists $name]} continue
    foreach ext $extlist {
      if {[file exists $name.$ext]} {
        set name ""
        break
      }
    }
    if {$name != ""} {return $name}
  }
  return ""
}

