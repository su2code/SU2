
proc color_generate {index} {
  set index [expr abs($index) % 132]
  set i [expr $index % 12]
  set h [lindex {0 1 2 3 4 5 .5 1.25 2.65 3.4 4.5 5.5} $i]
  set i [expr int($h)]
  set h [expr $h - double($i)]
  set j [expr $index / 12]
  if {[expr $j % 2] == 0} {
    set v 1.0
    set s [expr 1.0 - sqrt(double($j) / 22.0)]
  } else {
    set v [expr 1.0 - sqrt(double($j) / 44.0)]
    set s 1.0
  }
  switch $i {
    0 {
        set r $v
        set g [expr $v * (1.0 - ($s * (1.0 - $h)))]
        set b [expr $v * (1.0 - $s)]
      }
    1 {
        set r [expr $v * (1.0 - ($s * $h))]
        set g $v
        set b [expr $v * (1.0 - $s)]
      }
    2 {
        set r [expr $v * (1.0 - $s)]
        set g $v
        set b [expr $v * (1.0 - ($s * (1.0 - $h)))]
      }
    3 {
        set r [expr $v * (1.0 - $s)]
        set g [expr $v * (1.0 - ($s * $h))]
        set b $v
      }
    4 {
        set r [expr $v * (1.0 - ($s * (1.0 - $h)))]
        set g [expr $v * (1.0 - $s)]
        set b $v
      }
    5 {
        set r $v
        set g [expr $v * (1.0 - $s)]
        set b [expr $v * (1.0 - ($s * $h))]
      }
    6 {
        set r $v
        set g [expr $v * (1.0 - $s)]
        set b [expr $v * (1.0 - $s)]
      }
  }

  if [expr $r < 0.0] {set r 0.0}
  if [expr $r > 1.0] {set r 1.0}
  if [expr $g < 0.0] {set g 0.0}
  if [expr $g > 1.0] {set g 1.0}
  if [expr $b < 0.0] {set b 0.0}
  if [expr $b > 1.0] {set b 1.0}

  return [list $r $g $b]
}

proc color_value {rgb} {
  return [format "#%2.2x%2.2x%2.2x" \
    [expr int(255.0 * [lindex $rgb 0])] \
    [expr int(255.0 * [lindex $rgb 1])] \
    [expr int(255.0 * [lindex $rgb 2])]]
}

proc color_lighten {rgb {amt 0.25}} {
  set newclr ""
  foreach c $rgb {
    set nc [expr $c + $amt]
    if [expr $nc < 0.0] {set nc 0.0}
    if [expr $nc > 1.0] {set nc 1.0}
    lappend newclr $nc
  }
  return $newclr
}

proc color_darken {rgb {amt 0.25}} {
  set newclr ""
  foreach c $rgb {
    set nc [expr $c - $amt]
    if [expr $nc < 0.0] {set nc 0.0}
    if [expr $nc > 1.0] {set nc 1.0}
    lappend newclr $nc
  }
  return $newclr
}

proc color_gray {rgb} {
  set gray [expr 0.30 * [lindex $rgb 0] + \
                 0.59 * [lindex $rgb 1] + \
                 0.11 * [lindex $rgb 2]]
  if [expr $gray < 0.0] {return 0.0}
  if [expr $gray > 1.0] {return 1.0}
  return $gray
}

array set _ColorData {
  win ""
  done 0
  red 127
  green 127
  blue 127
  red,y 0
  green,y 0
  blue,y 0
  colorbars 16
  arrowsize 5
  hue {0 6 1 2 3 9 10 5}
  sat {10 6 2 0 1 5 9}
}

proc color_select {w title {oldclr ""} {loc ""}} {
  global _ColorData
    
  if {[llength $oldclr] == 3} {
    set r [expr int(255.0 * [lindex $oldclr 0])]
    if {$r < 0} {
      set _ColorData(red) 0
    } elseif {$r > 255} {
      set _ColorData(red) 255
    } else {
      set _ColorData(red) $r
    }
    set g [expr int(255.0 * [lindex $oldclr 1])]
    if {$g < 0} {
      set _ColorData(green) 0
    } elseif {$g > 255} {
      set _ColorData(green) 255
    } else {
      set _ColorData(green) $g
    }
    set b [expr int(255.0 * [lindex $oldclr 2])]
    if {$b < 0} {
      set _ColorData(blue) 0
    } elseif {$b > 255} {
      set _ColorData(blue) 255
    } else {
      set _ColorData(blue) $b
    }
  }

  set _ColorData(win) $w
  catch {destroy $w}
  toplevel $w
  wm protocol $w WM_DELETE_WINDOW {set _ColorData(done) 0}
  wm title $w $title
  wm transient $w [winfo toplevel [winfo parent $w]]
  wm resizable $w 0 0

  frame $w.l
  pack $w.l -side left -fill y -expand yes -padx 5 -pady 5

  frame $w.l.t
  pack $w.l.t -side top -pady 5 -anchor n

  for {set j 0} {$j < 7} {incr j} {
    set f [frame $w.l.t.y$j]
    pack $f -side top -fill x
    set sat [lindex $_ColorData(sat) $j]
    for {set i 0} {$i < 8} {incr i} {
      set hue [lindex $_ColorData(hue) $i]
      set n [expr $hue + 12 * $sat]
      set rgb [color_generate $n]
      set clr [color_value $rgb]
      button $f.b$i -padx 3m -pady 1m -bg $clr -activebackground $clr \
        -command "color:setrgb $rgb"
      pack $f.b$i -side left
    }
  }
  set f [frame $w.l.t.y8]
  pack $f -side top -fill x
  for {set i 0} {$i < 8} {incr i} {
    set g [expr double($i) / 7.0]
    set rgb [list $g $g $g]
    set clr [color_value $rgb]
    button $f.b$i  -padx 3m -pady 1m -bg $clr -activebackground $clr \
      -command "color:setrgb $rgb"
    pack $f.b$i -side left
  }
  
  frame $w.l.b
  pack $w.l.b -side bottom -fill both -expand yes
  
  set f [frame $w.l.b.b]
  pack $f -side left -padx 5 -pady 5 -fill y -expand yes
  button $f.a -text Accept -underline 0 -padx 3m -pady 1m \
    -command {set _ColorData(done) 1}
  button $f.c -text Cancel -underline 0 -padx 3m -pady 1m \
    -command {set _ColorData(done) 0}
  pack $f.a $f.c -side top -expand yes
  
  bind $w <Alt-a> "$f.a flash; set _ColorData(done) 1"
  bind $w <Alt-c> "$f.c flash; set _ColorData(done) 0"

  frame $w.l.b.s -relief sunken -bd 2 -width 100 -height 50 \
    -bg [format "#%2.2x%2.2x%2.2x" \
       $_ColorData(red) $_ColorData(green) $_ColorData(blue)]
  pack $w.l.b.s -side right -expand yes -fill both -padx 2 -pady 5

  frame $w.r
  pack $w.r -side right -fill y -expand yes -padx 5 -pady 10

  set width [expr {[winfo reqwidth $w.l.t.y0.b0] - \
      2*([$w.l.t.y0.b0 cget -highlightthickness] + \
         [$w.l.t.y0.b0 cget -bd])}]
  set height [expr 8 * [winfo reqheight $w.l.t.y0.b0]]
   
  foreach {c l} {red R green G blue B} {
    set f [frame $w.r.$c]
    if {$c == "green"} {
      pack $f -side left -fill y -expand yes -padx 5
    } else {
      pack $f -side left -fill y -expand yes
    }

    set box [frame $f.box]
    pack $box -side bottom -fill x -expand yes

    label $box.lab -text $l
    entry $box.ent -width 4 -textvariable _ColorData($c)
    pack $box.lab $box.ent -side top -fill x -expand yes

    bind $box.ent <Return> [list color:entryvalue $c]

    canvas $f.color -height $height -width $width -relief sunken -bd 2
    pack $f.color -side left -expand yes -fill both
    canvas $f.sel -height $height -width [expr 2 * $_ColorData(arrowsize)] \
      -highlightthickness 0
    pack $f.sel -side right -expand yes -fill y
    
    bind $f.color <Configure> "color:drawscale $c"

    bind $f.color <Enter> "color:selColor $f $c red"
    bind $f.color <Leave> "color:selColor $f $c black"
    bind $f.color <ButtonPress-1> "color:Move $f $c %y"
    bind $f.color <B1-Motion> "color:Move $f $c %y"
    bind $f.color <ButtonRelease-1> "color:endMove $f $c %y"

    bind $f.sel <Enter> "color:selColor $f $c red"
    bind $f.sel <Leave> "color:selColor $f $c black"
    bind $f.sel <ButtonPress-1> "color:Move $f $c %y"
    bind $f.sel <B1-Motion> "color:Move $f $c %y"
    bind $f.sel <ButtonRelease-1> "color:endMove $f $c %y"
  }

  if {$loc != ""} {center_window $w $loc}

  set oldFocus [focus]
  set oldGrab [grab current $w]
  if {$oldGrab != ""} {
    set grabStatus [grab status $oldGrab]
  }
  catch {grab $w}
  tkwait visibility $w
  raise $w
  focus $w
  tkwait variable _ColorData(done)
  catch {focus $oldFocus}
  destroy $w
  if {$oldGrab != ""} {
    if {$grabStatus == "global"} {
      grab -global $oldGrab
    } else {
      grab $oldGrab
    }
  }

  if {$_ColorData(done)} {
    set r [format "%.3f" [expr double($_ColorData(red)) / 255.0]]
    set g [format "%.3f" [expr double($_ColorData(green)) / 255.0]]
    set b [format "%.3f" [expr double($_ColorData(blue)) / 255.0]]
    return [list $r $g $b]
  }
  return ""
}

proc color:setrgb {r g b} {
  global _ColorData

  set _ColorData(red) [expr int(255.0 * $r)]
  set _ColorData(green) [expr int(255.0 * $g)]
  set _ColorData(blue) [expr int(255.0 * $b)]

  color:updateall
}

proc color:entryvalue {clr} {
  global _ColorData
  
  if {[catch {
    set val [expr int($_ColorData($clr))]}
      ]} {set val 0}
  if {$val < 0} {set val 0}
  if {$val > 255} {set val 255}
  set _ColorData($clr) $val

  color:updateall
}

proc color:drawcolor {} {
  global _ColorData

  set color [format "#%2.2x%2.2x%2.2x" \
    $_ColorData(red) $_ColorData(green) $_ColorData(blue)]
  $_ColorData(win).l.b.s configure -bg $color
}

proc color:drawscale {clr} {
  global _ColorData
  
  set dc [expr 255.0 / double($_ColorData(colorbars))]
  set bar $_ColorData(win).r.$clr.color
  $bar delete all
  set height [winfo height $bar]
  set dx [winfo width $bar]
  set dy [expr double($height) / $_ColorData(colorbars)]
  for {set i 0} {$i < $_ColorData(colorbars)} {incr i} {
    set y [expr $i * $dy]
    set c [expr 255 - int($i * $dc)]
    if {$clr == "red"} {
      set color [format "#%2.2x%2.2x%2.2x" \
                 $c $_ColorData(green) $_ColorData(blue)]
    } elseif {$clr == "green"} {
      set color [format "#%2.2x%2.2x%2.2x" \
                 $_ColorData(red) $c $_ColorData(blue)]
    } else {
      set color [format "#%2.2x%2.2x%2.2x" \
                 $_ColorData(red) $_ColorData(green) $c]
    }
    $bar create rect 0 $y $dx [expr $y + $dy] \
      -fill $color -outline $color
  }

  set sel $_ColorData(win).r.$clr.sel
  set y [expr $_ColorData(arrowsize) + \
         ($height - 2 * $_ColorData(arrowsize)) * \
         (1.0 - $_ColorData($clr) / 255.0)]
  set _ColorData($clr,y) $y
  $sel delete all
  set _ColorData($clr,sel) [$sel create polygon 0 $y \
    [expr 2 * $_ColorData(arrowsize)] [expr $y + $_ColorData(arrowsize)] \
    [expr 2 * $_ColorData(arrowsize)] [expr $y - $_ColorData(arrowsize)]]
}

proc color:updateall {} {
  color:drawcolor
  color:drawscale red  
  color:drawscale green  
  color:drawscale blue  
}

proc color:selColor {f clr selclr} {
  global _ColorData

  $f.sel itemconfigure $_ColorData($clr,sel) -fill $selclr
}

proc color:Move {f clr y} {
  global _ColorData

  set bar $f.color
  set height [winfo height $bar]
  if {$y < 0} {set y 0}
  if {$y > $height} {set y $height}
  set c [expr int(255.0 * (1.0 - double($y) / double($height)))]
  if {$c < 0} {set c 0}
  if {$c > 255} {set c 255}
  set _ColorData($clr) $c
  color:drawcolor

  incr height -$_ColorData(arrowsize)
  if {$y < $_ColorData(arrowsize)} {set y $_ColorData(arrowsize)}
  if {$y > $height} {set y $height}
  set diff [expr $y - $_ColorData($clr,y)]
  $f.sel move $_ColorData($clr,sel) 0 $diff
  set _ColorData($clr,y) [expr $_ColorData($clr,y) + $diff]
}

proc color:endMove {f clr y} {
  color:Move $f $clr $y
  color:drawscale red  
  color:drawscale green  
  color:drawscale blue  
}

proc color:allcolors {{w .allcolor}} {
  catch {destroy $w}
  toplevel $w
  wm title $w "All Colors"

  set n 0
  for {set j 0} {$j < 11} {incr j} {
    set f [frame $w.y$j]
    pack $f -side top -fill x
    for {set i 0} {$i < 12} {incr i} {
      set clr [color_value [color_generate $n]]
      button $f.b$i -bg $clr -activebackground $clr
      pack $f.b$i -side left
      incr n
    }
  }
}


