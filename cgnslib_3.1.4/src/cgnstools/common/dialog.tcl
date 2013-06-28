
if {![info exists UseNativeDialogs] || $UseNativeDialogs == ""} {
  if {$tcl_platform(platform) == "windows"} {
    set UseNativeDialogs 1
  } else {
    set UseNativeDialogs 0
  }
}

#----- pop-up dialog window - borrowed from tk_dialog

proc dialog {w x y title text bitmap default args} {
  global _DialogDone
  catch {destroy $w}
  toplevel $w -class Dialog
  wm title $w $title
  if {$x != {} && [winfo exists $x] && [winfo ismapped $x]} {
    wm transient $w [winfo toplevel $x]
  } else {
    wm transient $w [winfo toplevel [winfo parent $w]]
  }

  frame $w.top -relief raised -bd 1
  pack $w.top -side top -fill both

  message $w.msg -width 10c -text $text
  pack $w.msg -in $w.top -side right -expand 1 -fill both -padx 5 -pady 5
  if {$bitmap != ""} {
    if {![catch {image type $bitmap} type] && $type == "photo"} {
      label $w.bitmap -image $bitmap
    } else {
      label $w.bitmap -bitmap $bitmap
    }
    pack $w.bitmap -in $w.top -side left -padx 5 -pady 5
  }

  set nbut 0
  set wbut 0
  set fbut $w
  foreach but $args {
    if {$but != ""} {
      incr nbut
      if {$wbut < [string length $but]} {
        set wbut [string length $but]
      }
    }
  }
  if {$nbut} {
    if {$wbut < 4} {set wbut 4}
    frame $w.bot -relief raised -bd 1
    pack $w.bot -side bottom -fill both
    set i 0
    foreach but $args {
      if {$but != ""} {
        button $w.button$i -text $but -width $wbut -command "set _DialogDone $i"
        pack $w.button$i -in $w.bot -side left -expand 1 -padx 3 -pady 5
        if {$i == $default} {
          set fbut $w.button$i
          $fbut configure -default active
          bind $w <Return> "$fbut flash; set _DialogDone $i"
        }
      }
      incr i
    }
  }

  center_window $w $x $y

  if {!$nbut} {
    update idletasks
    return 0
  }

  set oldFocus [focus]
  set oldGrab [grab current $w]
  if {$oldGrab != ""} {
    set grabStatus [grab status $oldGrab]
  }
  catch {grab $w}
  tkwait visibility $w
  focus $fbut
  tkwait variable _DialogDone
  catch {focus $oldFocus}
  destroy $w
  if {$oldGrab != ""} {
    if {$grabStatus == "global"} {
      grab -global $oldGrab
    } else {
      grab $oldGrab
    }
  }
  return $_DialogDone
}

proc MessageBox {title msg {icon ""} {type ""} {default ""}} {
  global UseNativeDialogs
  if {$UseNativeDialogs} {
    set cmd "tk_messageBox -title {$title} -message {$msg}"
    if {$icon != ""} {
      append cmd " -icon $icon"
    }
    if {$type != ""} {
      append cmd " -type $type"
    }
    if {$default != ""} {
      append cmd " -default $default"
    }
    return [eval $cmd]
  }

  set buttons OK
  case $type {
    abortretryignore {set buttons "Abort Retry Ignore"}
    okcancel {set buttons "OK Cancel"}
    retrycancel {set buttons "Retry Cancel"}
    yesno {set buttons "Yes No"}
    yesnocancel {set buttons "Yes No Cancel"}
  }
  set defbut 0
  if {$default != ""} {
    set n 0
    foreach b [split $buttons] {
      if {$default == [string tolower $b]} {
        set defbut $n
        break
      }
      incr n
    }
  }
  set cmd "dialog .messagewin {} {} {$title} {$msg} \
    {$icon} $defbut $buttons"
  set n [eval $cmd]
  if [catch {lindex $buttons $n} result] {
    return ""
  }
  return [string tolower $result]
}

proc DialogBox {w x y title text bitmap args} {
  global _DialogDone

  # get button width

  set width 0
  foreach btn $args {
    if {[string length [lindex $btn 0]] > $width} {
      set width [string length [lindex $btn 0]]
    }
  }
  if !$width {
    return [dialog $w $x $y $title $text $bitmap {} {}]
  }
  incr width 2

  catch {destroy $w}
  toplevel $w -class Dialog
  wm title $w $title
  if {$x != {} && [winfo exists $x] && [winfo ismapped $x]} {
    wm transient $w [winfo toplevel $x]
  } else {
    wm transient $w [winfo toplevel [winfo parent $w]]
  }

  # top message

  frame $w.msg -relief raised -bd 1
  pack $w.msg -side top -fill both -expand 1
  if {$bitmap != ""} {
    if {![catch {image type $bitmap} type] && $type == "photo"} {
      label $w.msg.icon -image $bitmap
    } else {
      label $w.msg.icon -bitmap $bitmap
    }
    pack $w.msg.icon -side left -padx 10 -pady 3
  }
  message $w.msg.msg -text $text -width 8c
  pack $w.msg.msg -pady 3 -anchor w

  # button selections

  set cnt -1
  foreach btn $args {
    incr cnt
    set f [frame $w.f$cnt -relief raised -bd 1]
    pack $f -side top -fill x
    button $f.btn -width $width -text [lindex $btn 0] \
      -padx 3 -pady 2 -command "set _DialogDone $cnt"
    pack $f.btn -side left -padx 3 -pady 3
    pack $f.btn -side left
    message $f.msg -text [lindex $btn 1] -width 10c
    pack $f.msg -anchor w -padx 3 -pady 3
  }

  center_window $w $x $y

  set oldFocus [focus]
  set oldGrab [grab current $w]
  if {$oldGrab != ""} {
    set grabStatus [grab status $oldGrab]
  }
  catch {grab $w}
  focus $w
  tkwait variable _DialogDone
  catch {focus $oldFocus}
  destroy $w
  if {$oldGrab != ""} {
    if {$grabStatus == "global"} {
      grab -global $oldGrab
    } else {
      grab $oldGrab
    }
  }
  return $_DialogDone
}

#----- error message popup window

proc errormsg {msg {x -1} {y -1}} {
  if {$msg != ""} {
    dialog .errorwin $x $y Error $msg error 0 Dismiss
  }
}

proc ErrorMessage {msg} {
  MessageBox Error "$msg" error
}

#----- center a window

proc center_window {w xref {yref -1}} {
  wm withdraw $w
  update idletasks
  set ww [winfo reqwidth $w]
  set wh [winfo reqheight $w]

  if {$xref == ""} {
    set x [expr ([winfo screenwidth  $w] - $ww) / 2]
    set y [expr ([winfo screenheight $w] - $wh) / 2]
  } elseif {[winfo exists $xref] && [winfo ismapped $xref]} {
    set x [expr [winfo rootx $xref] + ([winfo width  $xref] - $ww) / 2]
    set y [expr [winfo rooty $xref] + ([winfo height $xref] - $wh) / 2]
  } else {
    if [catch {expr $xref < 0} x] {
      set x 1
    }
    if {$x} {
      if [winfo ismapped .] {
        set x [expr [winfo rootx .] + ([winfo width .] - $ww) / 2]
      } else {
        set x [expr ([winfo screenwidth $w] - $ww) / 2]
      }
    } else {
      set x [expr $xref - $ww / 2]
    }
    if [catch {expr $yref < 0} y] {
      set y 1
    }
    if {$y} {
      if [winfo ismapped .] {
        set y [expr [winfo rooty .] + ([winfo height .] - $wh) / 2]
      } else {
        set y [expr ([winfo screenheight $w] - $wh) / 2]
      }
    } else {
      set y [expr $yref - $wh / 2]
    }
  }

  if {$x < 0} {
    set pos +0
  } elseif {[expr $x + $ww] > [winfo screenwidth $w]} {
    set pos -0
  } else {
    set pos +$x
  }

  if {$y < 0} {
    set pos $pos+0
  } elseif {[expr $y + $wh] > [winfo screenheight $w]} {
    set pos $pos-0
  } else {
    set pos $pos+$y
  }

  wm geometry $w $pos
  update idletasks
  wm deiconify $w
}

#----- about message window

proc about {title text {bitmap ""}} {
  dialog .about -1 -1 $title $text $bitmap 0 Close
}

#----- popup message

proc popup_message {msg args} {
  global Font

  set font $Font(normal)
  set bg #ffffcc
  set fg black
  set parent .
  set pos ""
  set width 5c
  set wrap 1
  foreach {opt val} $args {
    switch -glob -- $opt {
      -fon* {set font $val}
      -par* {set parent $val}
      -pos* {set pos $val}
      -for* - -fg {set fg $val}
      -bac* - -bg {set bg $val}
      -wid* {set width $val}
      -wra* {set wrap $val}
    }
  }
  if {$pos == ""} {set pos $parent}

  set w .popup
  catch {destroy $w}
  toplevel $w -bg black
  wm overrideredirect $w 1
  wm transient $w [winfo toplevel $parent]
  if {$wrap} {
    message $w.l -text $msg -font $font -relief flat -bg $bg -fg $fg \
      -padx 2 -pady 0 -anchor w -width $width
  } else {
    label $w.l -text $msg -font $font -relief flat -bg $bg -fg $fg \
      -padx 2 -pady 0 -anchor w -justify left -wraplength 0
  }
  pack $w.l -side left -padx 1 -pady 1
  eval center_window $w $pos

  bind $w <ButtonRelease> {catch {destroy .popup};break}
  bind $w <KeyRelease> {catch {destroy .popup};break}
  bind $w <FocusOut> {catch {destroy .popup}}

  set oldFocus [focus]
  set oldGrab [grab current $w]
  if {$oldGrab != ""} {
    set grabStatus [grab status $oldGrab]
  }
  catch {grab $w}
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
}

