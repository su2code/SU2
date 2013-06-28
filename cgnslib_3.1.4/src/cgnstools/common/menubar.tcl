set _MenuData(next) 0

proc menubar_create {menulist {top .}} {
  global _MenuData tcl_platform bgColor

  if {$top == "." || $top == ""} {
    set mw .menubar
    set ms .menusep
  } else {
    set mw $top.menubar
    set ms $top.menusep
  }

  set _MenuData(menu) $mw
  set _MenuData(list) ""

  if {$tcl_platform(platform) == "windows"} {
    menu $mw -tearoff 0 -relief flat -type menubar
    foreach j $menulist {
      set i [string tolower $j]
      lappend _MenuData(list) $i
      set m m$_MenuData(next)
      incr _MenuData(next)
      $mw add cascade -label $j -menu $mw.$m -underline 0
      menu $mw.$m -tearoff 0
      set _MenuData($i,menu) $mw.$m
    }
    . configure -menu $mw
  } else {
    frame $mw -relief flat -bg $bgColor(button)
    pack $mw -side top -fill x
    foreach j $menulist {
      set i [string tolower $j]
      lappend _MenuData(list) $i
      set m m$_MenuData(next)
      incr _MenuData(next)
      menubutton $mw.$m -text $j -menu $mw.$m.menu -pady 0 \
        -highlightthickness 0 -underline 0
      pack $mw.$m -side left -padx 5
      menu $mw.$m.menu -tearoff 0
      set _MenuData($i,menu) $mw.$m.menu
    }
  }

  frame $ms -bd 1 -height 2 -relief sunken
  pack $ms -side top -fill x
}

proc menubar_add {menu {before ""}} {
  global tcl_platform _MenuData bgColor

  set i [string tolower $menu]
  if {[info exists _MenuData($i,menu)]} {
    return $_MenuData($i,menu)
  }
  set where ""
  if {$before != ""} {
    set before [string tolower $before]
    set n [lsearch $_MenuData(list) $before]
    if {$n >= 0} {set where $n}
  }
  set mw $_MenuData(menu)
  set m m$_MenuData(next)
  incr _MenuData(next)

  if {$tcl_platform(platform) == "windows"} {
    if {$where == ""} {
      $mw add cascade -label $menu -menu $mw.$m -underline 0
    } else {
      $mw insert $where cascade -label $menu -menu $mw.$m -underline 0
    }
    set _MenuData($i,menu) $mw.$m
  } else {
    menubutton $mw.$m -text $menu -menu $mw.$m.menu -pady 0 \
      -highlightthickness 0 -underline 0
    if {$where == ""} {
      pack $mw.$m -side left -padx 5
    } else {
      pack $mw.$m -side left -padx 5 \
        -before [winfo parent $_MenuData($before,menu)]
    }
    set _MenuData($i,menu) $mw.$m.menu
  }

  menu $_MenuData($i,menu) -tearoff 0
  if {$where == ""} {
    lappend _MenuData(list) $i
  } else {
    set _MenuData(list) [linsert $_MenuData(list) $where $i]
  }
  return $_MenuData($i,menu)
}

proc menubar_delete {menu} {
  global _MenuData tcl_platform

  set i [string tolower $menu]
  set n [lsearch $_MenuData(list) $i]
  if {$n < 0} return
  if {$tcl_platform(platform) == "windows"} {
    $_MenuData(menu) delete $n $n
  } else {
    set m [winfo parent $_Menudata($i,menu)]
    pack forget $m
    destroy $m
  }
  set _MenuData(list) [lreplace $_MenuData(list) $n $n]
}

proc menubar_get {menu} {
  global _MenuData
  set i [string tolower $menu]
  if {[lsearch $_MenuData(list) $i] >= 0} {
    return $_MenuData($i,menu)
  }
  return ""
}

proc menubar_state {menu state {entry ""}} {
  global _MenuData tcl_platform

  set i [string tolower $menu]
  set n [lsearch $_MenuData(list) $i]
  if {$n < 0} return
  if {$entry != ""} {
    $_MenuData($i,menu) entryconfigure $entry -state $state
  } elseif {$tcl_platform(platform) == "windows"} {
    $_MenuData(menu) entryconfigure $n -state $state
  } else {
    [winfo parent $_MenuData($i,menu)] configure -state $state
  }
}

