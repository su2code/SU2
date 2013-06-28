
array set Balloon {
  set 0
  first 0
  id ""
}

proc set_balloon {target message} {
  global Balloon
  set tags [bindtags $target]
  set n [lsearch $tags Balloon]
  if {$message == ""} {
    if {$n >= 0} {
      bindtags $target [lreplace $tags $n $n]
    }
  } else {
    if {$n < 0} {
      bindtags $target "Balloon $tags"
    }
  }
  set Balloon($target) $message
}

proc end_balloon {target} {
  set_balloon $target ""
}

bind Balloon <Enter> {
  set Balloon(set) 0
  set Balloon(first) 1
  set Balloon(id) [after 500 {Balloon:show %W $Balloon(%W) %X %Y}]
}

bind Balloon <Button> {
  set Balloon(first) 0
  Balloon:kill
}

bind Balloon <Leave> {
  set Balloon(first) 0
  Balloon:kill
}

bind Balloon <Motion> {
  if {$Balloon(set) == 0} {
    after cancel $Balloon(id)
    set Balloon(id) [after 500 {Balloon:show %W $Balloon(%W) %X %Y}]
  }
}

proc Balloon:kill {} {
  global Balloon
  after cancel $Balloon(id)
  if {[winfo exists .balloon] == 1} {
      destroy .balloon
  }
  set Balloon(set) 0
}

proc Balloon:show {target message {cx 0} {cy 0}} {
  global Balloon
  if {$Balloon(first) == 1 } {
    set Balloon(first) 2
    if {$cx == 0 && $cy == 0} {
      set x [expr [winfo rootx $target] + ([winfo width $target]/2)]
      set y [expr [winfo rooty $target] + [winfo height $target] + 4]
    } else {
      set x [expr $cx + 4]
      set y [expr $cy + 4]
    }
    toplevel .balloon -bg black
    wm overrideredirect .balloon 1
    label .balloon.l -text $message -relief flat \
      -bg #ffffcc -fg black -padx 2 -pady 0 -anchor w
    pack .balloon.l -side left -padx 1 -pady 1
    wm geometry .balloon +$x+$y
    set Balloon(set) 1
  }
}

proc entry_balloon {target} {
  global Balloon
  bind $target <Enter> {
    incr Balloon(first)
    set xv [%W xview]
    if {[string compare [focus] %W] &&
      [expr [lindex $xv 1] - [lindex $xv 0] < 0.999]} {
      set Balloon(id) [after 500 {Balloon:entry %W}]
    }
  }
  bind $target <Leave> {
    catch {after cancel $Balloon(id)}
    if {![winfo exists .balloon]} {
      set Balloon(first) 0
    }
  }
  bind $target <Button> {
    catch {after cancel $Balloon(id)}
  }
}

proc Balloon:entry {target} {
  global Balloon
  if {$Balloon(first) == 1} {
    set x [winfo rootx $target]
    set y [winfo rooty $target]
    toplevel .balloon -bg black
    wm overrideredirect .balloon 1
    label .balloon.l -text [$target get] -relief flat -bg #ffffcc -fg black \
      -padx 2 -pady 0 -anchor w -font [$target cget -font]
    pack .balloon.l -side left -padx 1 -pady 1
    wm withdraw .balloon
    update idletasks
    set w [winfo reqwidth .balloon]
    if {[expr $x + $w > [winfo screenwidth .balloon]]} {
      set x [expr [winfo screenwidth .balloon] - $w]
    }
    wm geometry .balloon +$x+$y
    update idletasks
    wm deiconify .balloon
    bind .balloon <Leave> {
      set Balloon(first) 0
      Balloon:kill
    }
    bind .balloon <Button> {Balloon:kill}
  }
}

