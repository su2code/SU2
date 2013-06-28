
#########################################################
# Directory Selector TCL version 1.1
#
# Originally written by:
# Daniel Roche, <dan@lectra.com>
#
# Modified for xmktclapp (and for version of Tk prior to 8.0) by:
# D. Richard Hipp, <drh@hwaci.com>

# tk_getDirectory [option value ...]
#
#  options are :
#   [-parent window]      parent window
#   [-initialdir dir]     display in dir
#   [-title string]       make string title of dialog window
#   [-ok string]          make string the label of OK button
#   [-cancel string]      make string the label of CANCEL button
#   [-label string]       make string the label of the directory message
#   [-nofiles]            don't show files
#

proc tk_getDirectory {args} {
  global tcl_platform tk_getDirectory

  set _titre "Directory Selector"
  set _ldir Directory:
  set _open Ok
  set _cancel Cancel
  set _parent {}
  set tk_getDirectory(curdir) [pwd]
  set tk_getDirectory(showfiles) 1

  set ind 0
  set max [llength $args]
  while { $ind < $max } {
    switch -exact -- [lindex $args $ind] {
      "-parent" {
        incr ind
        set _parent [lindex $args $ind]
        incr ind
      }
      "-initialdir" {
        incr ind
        if {![catch {cd [lindex $args $ind]}]} {
          set dir [pwd]
          catch {cd $tk_getDirectory(curdir)}
          set tk_getDirectory(curdir) $dir
        }
        incr ind
      }
      "-title" {
        incr ind
        set _titre [lindex $args $ind]
        incr ind
      }
      "-ok" {
        incr ind
        set _open [lindex $args $ind]
        incr ind
      }
      "-cancel" {
        incr ind
        set _cancel [lindex $args $ind]
        incr ind
      }
      "-label" {
        incr ind
        set _ldir [lindex $args $ind]
        incr ind
      }
      "-nofiles" {
        set tk_getDirectory(showfiles) 0
        incr ind
      }
      default {
        puts "unknown option [lindex $args $ind]"
        return ""
      }
    }
  }

  set tk_getDirectory(fini) 0

  if {![info exists tk_getDirectory:b_up]} {
    image create photo tk_getDirectory:b_up -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAQwQ4EMlJKwJvZaB7fsD1eSQolmV4kaqFaUA8fuG7\
xTht7S7e9r4gTsKhzTo1mO+YNKKaqFFNSP3xKhEAOw==}

    image create photo tk_getDirectory:b_new -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAQwQ1EMlJKwJv5a0BsNVHiZhmbiKohid5jVzsXl5t\
p+uE569lS5+g7vZ7BU9IIRDJ2QGYph2vEgEAOw==}

    image create photo tk_getDirectory:b_dir -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAA4ALAAAAAAQABAAQwQ+0MlJqzsPLc33PsD1jGQ5gtehruwaXh1yWBIW\
exxq47z80sAaK3i4yVq1R4+nUxo9zWUn+tSgHICsdrsNAiMAOw==}

    image create photo tk_getDirectory:b_file -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAQwQ3EMlJ6wTvkb0zqJjWfR+Ycd1TXkDrrqZIZF6M\
zir47rDlUyGaqicJ0l4x4TGpzAGbyyfU+atGAAA7}

    image create bitmap tk_getDirectory:b_down -data "
#define down_width 13
#define down_height 10
static unsigned char down_bits[] = {
   0x00, 0x00, 0x00, 0x00, 0xfe, 0x0f, 0xfc, 0x07, 0xf8, 0x03, 0xf0, 0x01,
   0xe0, 0x00, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00};"
  }

  set w .dirsel
  catch {destroy $w}

  toplevel $w
  wm geometry $w 420x220
  wm title $w $_titre
  if {$_parent == "" || ![winfo exists $_parent]} {
    set _parent [winfo toplevel [winfo parent $w]]
  } else {
    set _parent [winfo toplevel $_parent]
  }
  wm transient $w $_parent

  frame $w.f1 -relief flat -borderwidth 0
  frame $w.f2 -relief sunken -borderwidth 1
  frame $w.f3 -relief flat -borderwidth 0
  pack $w.f1 -fill x -padx 6 -pady 5
  pack $w.f2 -fill both -expand 1 -padx 6
  pack $w.f3 -fill x -padx 6 -pady 3

  label $w.f1.lab -text "Look In:"
  frame $w.f1.dir -relief sunken -bd 2
  entry $w.f1.dir.ent -relief flat -textvariable tk_getDirectory(curdir) \
    -state disabled -cursor {}
  menubutton $w.f1.dir.but -relief raised -image tk_getDirectory:b_down \
    -menu $w.f1.dir.but.m -direction left
  menu $w.f1.dir.but.m -tearoff 0
  pack $w.f1.dir.ent -side left -fill x -expand 1
  pack $w.f1.dir.but -side right -fill y
  button $w.f1.up -image tk_getDirectory:b_up \
    -command "tk_getDirectory:UpDir $w"
  button $w.f1.new -image tk_getDirectory:b_new \
    -command "tk_getDirectory:NewDir $w"
  pack $w.f1.lab -side left
  pack $w.f1.dir -side left -fill x -expand 1
  pack $w.f1.new -side right
  pack $w.f1.up -side right -padx 8

  bind $w.f1.dir.but <1> "tk_getDirectory:Menu $w"
  bind $w.f1.dir.ent <1> "tk_getDirectory:Menu $w;$w.f1.dir.but.m post %X %Y"

  canvas $w.f2.cv -borderwidth 0 -xscrollcommand "$w.f2.sb set" \
    -height 10 -bg white
  scrollbar $w.f2.sb -command "$w.f2.cv xview" -orient horizontal
  pack $w.f2.cv -side top -fill both -expand 1
  pack $w.f2.sb -side top -fill x

  $w.f2.cv bind TXT <Any-Button> "tk_getDirectory:HighlightItem $w"
  $w.f2.cv bind IMG <Any-Button> "tk_getDirectory:HighlightItem $w"
  $w.f2.cv bind TXT <Any-Double-Button> "tk_getDirectory:ClickItem $w"
  $w.f2.cv bind IMG <Any-Double-Button> "tk_getDirectory:ClickItem $w"

  label $w.f3.lab -text $_ldir
  entry $w.f3.ent -relief sunken -textvariable tk_getDirectory(seldir) \
    -state disabled -cursor {}
  pack $w.f3.lab -side left
  pack $w.f3.ent -side left -fill x -expand 1

  set width [string length $_open]
  if {$width < [string length $_cancel]} {
    set width [string length $_cancel]
  }
  if {$tcl_platform(platform) == "windows"} {
    button $w.f3.open -width $width -text $_open -pady 0 \
      -default active -command {set tk_getDirectory(fini) 1}
    button $w.f3.cancel -width $width -text $_cancel -pady 0 \
      -command {set tk_getDirectory(fini) 0}
  } else {
    button $w.f3.open -width $width -text $_open \
      -default active -command {set tk_getDirectory(fini) 1}
    button $w.f3.cancel -width $width -text $_cancel \
      -command {set tk_getDirectory(fini) 0}
  }
  pack $w.f3.cancel -side right
  pack $w.f3.open -side right -padx 8

  bind $w <Return> "$w.f3.open flash; set tk_getDirectory(fini) 1"

  set oldFocus [focus]
  set oldGrab [grab current $w]
  if {$oldGrab != ""} {
    set grabStatus [grab status $oldGrab]
  }
  catch {grab $w}

  if {[winfo ismapped $_parent]} {
    wm withdraw $w
    wm geometry $w "+[winfo rootx $_parent]+[winfo rooty $_parent]"
    update idletasks
    wm deiconify $w
  }
  update
  tk_getDirectory:ShowDir $w $tk_getDirectory(curdir)
  focus -force $w
  tkwait variable tk_getDirectory(fini)

  if { $tk_getDirectory(fini) == 1 } {
    set retval [eval file join [file split $tk_getDirectory(curdir)] \
      $tk_getDirectory(seldir)]
  } else {
    set retval ""
  }

  catch {focus $oldFocus}
  destroy $w
  if {$oldGrab != ""} {
    if {$grabStatus == "global"} {
      grab -global $oldGrab
    } else {
      grab $oldGrab
    }
  }
  return $retval
}

proc tk_getDirectory:ShowDir {w curdir} {
  global tk_getDirectory

  set tk_getDirectory(curdir) $curdir
  $w.f1.dir.ent xview moveto 1.0

  set wi [expr [image width tk_getDirectory:b_dir] + 3]
  set hi [image height tk_getDirectory:b_dir]
  set maxy [expr [winfo height $w.f2.cv]-$hi]

  set lidir [list]
  set lifile [list]
  foreach file [glob -nocomplain $curdir/*] {
    if [ file isdirectory $file ] {
      lappend lidir [file tail $file]
    } else {
      lappend lifile [file tail $file]
    }
  }
  set sldir [lsort $lidir]

  $w.f2.cv delete all
  set ind 0
  set x 2
  set y 2
  foreach file $sldir {
    $w.f2.cv create image $x $y -anchor nw \
      -image tk_getDirectory:b_dir -tags IMG
    $w.f2.cv create text [expr $x+$wi] $y -anchor nw -text $file -tags TXT
    incr y $hi
    if {$y >= $maxy} {
      set bbox [$w.f2.cv bbox all]
      set x [expr [lindex $bbox 2]+10]
      set y 2
    }
  }
  if {$tk_getDirectory(showfiles)} {
    foreach file [lsort $lifile] {
      $w.f2.cv create image $x $y -anchor nw -image tk_getDirectory:b_file
      $w.f2.cv create text [expr $x+$wi] $y -anchor nw -text $file
      incr y $hi
      if {$y >= $maxy} {
        set bbox [$w.f2.cv bbox all]
        set x [expr [lindex $bbox 2]+10]
        set y 2
      }
    }
  }
  $w.f2.cv configure -scrollregion [$w.f2.cv bbox all]
  set tk_getDirectory(seldir) ""
}

proc tk_getDirectory:UpDir {w} {
  global tk_getDirectory
  set curdir $tk_getDirectory(curdir)
  set curlst [file split $curdir]
  set nbr [llength $curlst]
  if {$nbr < 2} return
  set tmp [expr $nbr - 2]
  set newlst [ lrange $curlst 0 $tmp ]
  set newdir [ eval file join $newlst ]
  tk_getDirectory:ShowDir $w $newdir
}

proc tk_getDirectory:HighlightItem {w} {
  global tk_getDirectory
  catch {$w.f2.cv select clear}
  set id [$w.f2.cv find withtag current]
  if {[$w.f2.cv type $id] != "text"} {incr id}
  $w.f2.cv select from $id 0
  $w.f2.cv select to $id end
  set tk_getDirectory(seldir) [$w.f2.cv itemcget $id -text]
}

proc tk_getDirectory:ClickItem {w} {
  global tk_getDirectory
  set id [$w.f2.cv find withtag current]
  if {[$w.f2.cv type $id] != "text"} {incr id}
  set dir [$w.f2.cv itemcget $id -text]
  if {[string length $dir]==0} return
  tk_getDirectory:ShowDir $w [file join $tk_getDirectory(curdir) $dir]
}

proc tk_getDirectory:Menu {w} {
  global tk_getDirectory tcl_platform
  if {![info exist tk_getDirectory(drives)]} {
    if {$tcl_platform(platform) == "unix" ||
      [catch {file volume} tk_getDirectory(drives)]} {
      set tk_getDirectory(drives) {}
    }
    if {$tcl_platform(platform) == "windows"} {
      set tk_getDirectory(drives) [string toupper $tk_getDirectory(drives)]
    }
  }
  set curlst [file split $tk_getDirectory(curdir)]
  set nbr [llength $curlst]
  $w.f1.dir.but.m delete 0 last
  incr nbr -2
  set tmpdir {}
  for {set ind $nbr} {$ind >= 0} {incr ind -1} {
    set tmplst [ lrange $curlst 0 $ind]
    set tmpdir [ eval file join $tmplst]
    $w.f1.dir.but.m add command -label $tmpdir \
      -command "tk_getDirectory:ShowDir $w [list $tmpdir]"
  }
  set rootdir [string toupper $tmpdir]
  foreach drive $tk_getDirectory(drives) {
    if {$drive != $rootdir} {
      $w.f1.dir.but.m add command -label $drive \
        -command "tk_getDirectory:ShowDir $w [list $drive]"
    }
  }
}

proc tk_getDirectory:NewDir {wref} {
  global tcl_platform tk_getDirectory
  set w .newdir
  catch {destroy $w}
  toplevel $w
  wm title $w "New Folder"
  wm transient $w $wref

  set f [frame $w.name]
  pack $f -side top -fill x -padx 5 -pady 5
  label $f.lab -text "Name:"
  pack $f.lab -side left
  entry $f.ent -width 30
  pack $f.ent -side left -fill x -expand 1
  $f.ent insert 0 "New Folder"

  set f [frame $w.but]
  pack $f -side bottom -fill x -padx 5 -pady 5
  button $f.accept -text Ok -width 6 -default active \
    -command {set tk_getDirectory(newdir) 1}
  bind $w <Return> "
    $f.accept flash
    set tk_getDirectory(newdir) 1
  "
  button $f.cancel -text Cancel -width 6 \
    -command {set tk_getDirectory(newdir) 0}
  pack $f.accept $f.cancel -side left -expand 1
  if {$tcl_platform(platform) == "windows"} {
    $f.accept configure -pady 0
    $f.cancel configure -pady 0
  }

  wm withdraw $w
  wm geometry $w \
    "+[expr [winfo rootx $wref.f2.cv]+20]+[expr [winfo rooty $wref.f2.cv]+5]"
  update idletasks
  wm deiconify $w

  set oldFocus [focus]
  set oldGrab [grab current $w]
  if {$oldGrab != ""} {
    set grabStatus [grab status $oldGrab]
  }
  catch {grab $w}
  tkwait visibility $w
  $w.name.ent selection range 0 end
  focus $w.name.ent
  tkwait variable tk_getDirectory(newdir)

  if {$tk_getDirectory(newdir)} {
    set newdir "$tk_getDirectory(curdir)/[string trim [$w.name.ent get]]"
  } else {
    set newdir {}
  }

  catch {focus $oldFocus}
  destroy $w
  if {$oldGrab != ""} {
    if {$grabStatus == "global"} {
      grab -global $oldGrab
    } else {
      grab $oldGrab
    }
  }

  if {$newdir != ""} {
    if {[catch {file mkdir $newdir} msg]} {
      tk_dialog .error Error $msg error 0 Ok
    } else {
      tk_getDirectory:ShowDir $wref $tk_getDirectory(curdir)
    }
  }
}

