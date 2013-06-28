
#----- Create a new title frame widget

set _Combobox(arrow) [image create bitmap -data {
#define down_width  9
#define down_height 5
static char down_bits[] = {
  0xff, 0x01, 0xfe, 0x00, 0x7c, 0x00, 0x38, 0x00, 0x10, 0x00};}]

proc ComboboxCreate {w args} {
  global _Combobox
  frame $w -highlightthickness 0 -takefocus 0
  entry $w.ent -textvariable _Combobox($w,entry)
  pack $w.ent -side left -fill x -expand 1 -padx 0 -pady 0
  button $w.sel -image $_Combobox(arrow) -width 15 -padx 0 -pady 0 \
    -command "Combobox:post $w"
  pack $w.sel -side right -fill y

  bind $w.ent <Key-Up>    "Combobox:modify $w -"
  bind $w.ent <Key-Down>  "Combobox:modify $w +"
  bind $w.ent <Key-Prior> "Combobox:modify $w 0"
  bind $w.ent <Key-Next>  "Combobox:modify $w end"
  bind $w.ent <Alt-Up>    "Combobox:post $w"
  bind $w.ent <Alt-Down>  "Combobox:post $w"

  set _Combobox($w,bg) [$w.ent cget -bg]
  set _Combobox($w,fg) [$w.ent cget -fg]
  set _Combobox($w,font) [$w.ent cget -font]
  set _Combobox($w,width) ""
  set _Combobox($w,height) 10
  set _Combobox($w,relief) [$w.ent cget -relief]
  set _Combobox($w,bd) [$w.ent cget -bd]
  set _Combobox($w,index) ""
  set _Combobox($w,values) {}
  set _Combobox($w,edit) 1
  set _Combobox($w,listproc) ""
  set _Combobox($w,command) ""
  set _Combobox($w,post) ""
  set _Combobox($w,entry) ""
  set _Combobox($w,variable) _Combobox($w,entry)
  set _Combobox($w,state) normal
  set _Combobox($w,listfg) ""
  set _Combobox($w,listbg) ""

  eval ComboboxConfig $w $args
  return $w
}

#----- Change configuration options for the notebook widget

proc ComboboxConfig {w args} {
  global _Combobox

  # get configuration options

  foreach {tag value} $args {
    switch -- $tag {
      -bg {
        set _Combobox($w,bg) $value
      }
      -fg {
        set _Combobox($w,fg) $value
      }
      -font {
        set _Combobox($w,font) $value
      }
      -relief {
        set _Combobox($w,relief) $value
      }
      -bd {
        set _Combobox($w,bd) $value
      }
      -width {
        set _Combobox($w,width) $value
      }
      -height {
        set _Combobox($w,height) $value
      }
      -edit {
        set _Combobox($w,edit) $value
      }
      -values {
        set _Combobox($w,values) $value
      }
      -index {
        set _Combobox($w,index) $value
      }
      -listproc {
        set _Combobox($w,listproc) $value
      }
      -command {
        set _Combobox($w,command) $value
      }
      -post {
        set _Combobox($w,post) $value
      }
      -state {
        set _Combobox($w,state) $value
      }
      -variable {
        set _Combobox($w,variable) $value
      }
      -listfg {
        set _Combobox($w,listfg) $value
      }
      -listbg {
        set _Combobox($w,listbg) $value
      }
    }
  }

  if {$_Combobox($w,width) == ""} {
    if {$_Combobox($w,values) == {}} {
      set width [$w.ent cget -width]
    } else {
      set width 0
      foreach i $_Combobox($w,values) {
        if {$width < [string length $i]} {
          set width [string length $i]
        }
      }
    }
    set _Combobox($w,width) $width
  }

  $w configure -relief $_Combobox($w,relief) -bd $_Combobox($w,bd) \
    -bg $_Combobox($w,bg)
  $w.ent configure -fg $_Combobox($w,fg) -bg $_Combobox($w,bg) \
    -width $_Combobox($w,width) -font $_Combobox($w,font) \
    -textvariable $_Combobox($w,variable) \
    -relief flat -highlightthickness 0 -bd 0

  if {$_Combobox($w,state) == "disabled"} {
    $w.ent configure -state disabled -cursor {} -takefocus 0
    bind $w.ent <ButtonRelease-1> {}
    $w.sel configure -state disabled
  } else {
    if {$_Combobox($w,edit)} {
      $w.ent configure -state normal -cursor xterm -takefocus 1
      bind $w.ent <ButtonRelease-1> "Combobox:unpost $w"
    } else {
      $w.ent configure -state disabled -cursor hand2 -takefocus 1
      bind $w.ent <ButtonRelease-1> "$w.sel invoke"
    }
    if {$_Combobox($w,values) == {}} {
      $w.sel configure -state disabled
      set $_Combobox($w,variable) ""
    } else {
      $w.sel configure -state normal
      if {$_Combobox($w,index) != ""} {
        set n $_Combobox($w,index)
        if {$n >= 0 && $n < [llength $_Combobox($w,values)]} {
          set $_Combobox($w,variable) [lindex $_Combobox($w,values) $n]
        }
      }
    }
  }

  if {$_Combobox($w,listfg) == ""} {
    set _Combobox($w,listfg) $_Combobox($w,fg)
  }
  if {$_Combobox($w,listbg) == ""} {
    set _Combobox($w,listbg) $_Combobox($w,bg)
  }
}

proc ComboboxGet {w opt} {
  global _Combobox
  switch -- $opt {
    -width    {return $_Combobox($w,width)}
    -height   {return $_Combobox($w,height)}
    -bg       {return $_Combobox($w,bg)}
    -fg       {return $_Combobox($w,fg)}
    -font     {return $_Combobox($w,font)}
    -bd       {return $_Combobox($w,bd)}
    -relief   {return $_Combobox($w,relief)}
    -edit     {return $_Combobox($w,edit)}
    -values   {return $_Combobox($w,values)}
    -index    {return $_Combobox($w,index)}
    -listproc {return $_Combobox($w,listproc)}
    -command  {return $_Combobox($w,command)}
    -post     {return $_Combobox($w,post)}
    -variable {return $_Combobox($w,variable)}
    -state    {return $_Combobox($w,state)}
    -listfg   {return $_Combobox($w,listfg)}
    -listbg   {return $_Combobox($w,listbg)}
  }
  return -code error "combobox option $opt not known"
}

proc ComboboxInsert {w value {where end}} {
  global _Combobox
  set _Combobox($w,values) [linsert $_Combobox($w,values) $where $value]
  $w.sel configure -state normal
}

proc ComboboxDelete {w value} {
  global _Combobox
  set n [lsearch $_Combobox($w,values) $value]
  if {$n >= 0} {
    set _Combobox($w,values) [lreplace $_Combobox($w,values) $n $n]
    if {[llength $_Combobox($w,values)] == 0} {
      $w.sel configure -state disabled
    }
  }
}

proc ComboboxEntry {w} {
  return [$w.ent get]
}

proc ComboboxValue {w {n ""}} {
  global _Combobox
  if {$n == ""} {
    set n $_Combobox($w,index)
  }
  if {$n < 0 || $n >= [llength $_Combobox($w,values)]} {
    return ""
  }
  return [lindex $_Combobox($w,values) $n]
}

proc ComboboxIndex {w} {
  global _Combobox
  return $_Combobox($w,index)
}

proc Combobox:post {w} {
  global _Combobox
  if {[winfo exists $w.top]} {
    wm deiconify $top
    raise $top
    return
  }
  if {$_Combobox($w,listproc) != ""} {
    set _Combobox($w,values) [$_Combobox($w,listproc) $w]
  }
  if {$_Combobox($w,values) == {}} return

  set top [toplevel $w.top -relief flat -bd 0 -cursor left_ptr]
  wm overrideredirect $top 1
  wm withdraw $top
  wm transient $top [winfo toplevel $w]
  wm group $top [winfo toplevel $w]

  bind $top <ButtonRelease-1> "Combobox:unpost $w"
  bind $top <ButtonRelease-3> "Combobox:unpost $w"
  bind $top <Key-Up>          "Combobox:move $w -"
  bind $top <Key-Down>        "Combobox:move $w +"
  bind $top <Key-Home>        "Combobox:move $w 0"
  bind $top <Key-End>         "Combobox:move $w end"
  bind $top <Key-Prior>       "Combobox:page $w -"
  bind $top <Key-Next>        "Combobox:page $w +"
  bind $top <Return>          "Combobox:select $w"
  bind $top <space>           "Combobox:select $w"
  bind $top <Escape>          "Combobox:unpost $w"
  bind $top <FocusOut>        "Combobox:unpost $w"
  bind $top <Key-Tab>         "Combobox:select $w"

  set nv [llength $_Combobox($w,values)]
  set h $_Combobox($w,height)
  if {$h >= $nv} {
    set h $nv
  }

  set f [frame $top.f -relief solid -bd 1]
  pack $f -side top -fill both -expand 1

  listbox $f.list -relief flat -bd 0 -highlightthickness 0 \
    -exportselection 0 -height $h -font $_Combobox($w,font) \
    -fg $_Combobox($w,listfg) -bg $_Combobox($w,listbg) \
    -takefocus 0 -cursor hand2
  pack $f.list -side left -fill both -expand 1

  bind $f.list <ButtonRelease-1> "Combobox:select $w"

  if {$h < $nv} {
    scrollbar $f.ys -orient vertical -highlightthickness 0 \
      -takefocus 0 -command "$f.list yview" -width [$w.sel cget -width]
    pack $f.ys -side right -fill y
    $f.list configure -yscroll "$f.ys set"
    bind $f.ys <ButtonRelease-1> break
  }

  foreach i $_Combobox($w,values) {
    $f.list insert end $i
  }
  set s [$w.ent get]
  set len [llength $_Combobox($w,values)]
  set cmd $_Combobox($w,post)
  if {$cmd == "" || [catch {uplevel #0 $cmd $w $s} n] ||
    [catch {expr $n < 0 || $n >= $len} bad] || $bad} {
    set n [lsearch $_Combobox($w,values) $s]
  }
  if {$n >= 0 && $n < $len} {
    set _Combobox($w,index) $n
    $f.list selection set $n
    $f.list activate $n
    $f.list see $n
  }

  update idletasks
  set wx [winfo rootx $w]
  set wy [expr [winfo rooty $w] + [winfo height $w]]
  set ww [winfo width $w]
  set wh [winfo reqheight $top]
  wm geometry $top $ww\x$wh+$wx+$wy
  wm deiconify $top
  raise $top

  set _Combobox($w,focus) [focus]
  set _Combobox($w,grab) [grab current $top]
  if {$_Combobox($w,grab) != ""} {
    set _Combobox($w,status) [grab status $_Combobox($w,grab)]
  }
  catch {grab $top}
  catch {focus $top}
}

proc Combobox:unpost {w} {
  global _Combobox
  if {[winfo exists $w.top]} {
    focus $_Combobox($w,focus)
    destroy $w.top
    if {$_Combobox($w,grab) != ""} {
      if {$_Combobox($w,status) == "global"} {
        grab -global $_Combobox($w,grab)
      } else {
        grab $_Combobox($w,grab)
      }
    }
  }
}

proc Combobox:move {w pos} {
  global _Combobox
  set n [$w.top.f.list curselection]
  if {$n == ""} {set pos 0}
  set last [expr [llength $_Combobox($w,values)] - 1]
  if {$pos == "-"} {
    incr n -1
    if {$n < 0} {set n 0}
  } elseif {$pos == "+"} {
    incr n 1
    if {$n > $last} {set n $last}
  } elseif {$pos == "end"} {
    set n $last
  } else {
    set n $pos
  }
  $w.top.f.list selection clear 0 end
  $w.top.f.list selection set $n $n
  $w.top.f.list see $n
}

proc Combobox:page {w dir} {
  set sel [$w.top.f.list curselection]
  if {$dir == "-"} {
    set top [$w.top.f.list nearest 0]
    if {$top != $sel} {
      set pos $top
    } else {
      $w.top.f.list yview scroll -1 pages
      set pos [$w.top.f.list nearest 0]
    }
  } else {
    set h [expr [winfo height $w.top.f.list] - 4]
    set bot [$w.top.f.list nearest $h]
    if {$bot != $sel} {
      set pos $bot
    } else {
      $w.top.f.list yview scroll +1 pages
      set pos [$w.top.f.list nearest $h]
    }
  }
  $w.top.f.list selection clear 0 end
  $w.top.f.list selection set $pos $pos
}

proc Combobox:modify {w pos} {
  global _Combobox
  if {$_Combobox($w,values) == {}} return
  set last [expr [llength $_Combobox($w,values)] - 1]
  if {$pos == "-"} {
    set n [expr $_Combobox($w,index) - 1]
    if {$n < 0} {set n 0}
  } elseif {$pos == "+"} {
    set n [expr $_Combobox($w,index) + 1]
    if {$n > $last} {set n $last}
  } elseif {$pos == "end"} {
    set n $last
  } else {
    set n $pos
  }
  set _Combobox($w,index) $n
  set cmd $_Combobox($w,command)
  if {$cmd == "" || [catch {uplevel #0 $cmd $w $n} s] || $s == ""} {
    set s [lindex $_Combobox($w,values) $n]
  }
  set _Combobox(temp) $s
  uplevel #0 set $_Combobox($w,variable) {$_Combobox(temp)}
}

proc Combobox:select {w} {
  global _Combobox
  set sel [$w.top.f.list curselection]
  if {$sel != ""} {
    Combobox:modify $w $sel
  }
  Combobox:unpost $w
  return -code break
}

