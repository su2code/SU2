
set _Tree(font) {Helvetica 10}

set _Tree(maskdata) {
#define mask_width 9
#define mask_height 9
static unsigned char mask_bits[] = {
   0xff, 0x01, 0xff, 0x01, 0xff, 0x01,
   0xff, 0x01, 0xff, 0x01, 0xff, 0x01,
   0xff, 0x01, 0xff, 0x01, 0xff, 0x01};}

set _Tree(openbm) {
#define open_width 9
#define open_height 9
static unsigned char open_bits[] = {
   0xff, 0x01, 0x01, 0x01, 0x01, 0x01,
   0x01, 0x01, 0x7d, 0x01, 0x01, 0x01,
   0x01, 0x01, 0x01, 0x01, 0xff, 0x01};}

set _Tree(closedbm) {
#define closed_width 9
#define closed_height 9
static unsigned char closed_bits[] = {
   0xff, 0x01, 0x01, 0x01, 0x11, 0x01,
   0x11, 0x01, 0x7d, 0x01, 0x11, 0x01,
   0x11, 0x01, 0x01, 0x01, 0xff, 0x01};}

#----- create a new tree widget

proc TreeCreate {w args} {
  global _Tree

  canvas $w
  bind $w <Destroy> "TreeDestroy $w"

  foreach i {sort icon tags selection selidx build list} {
    set _Tree($w:$i) {}
  }
  set _Tree($w:root) 0
  set _Tree($w:open) 0
  set _Tree($w:padx) 0
  set _Tree($w:pady) 0
  set _Tree($w:indent) 12
  set _Tree($w:space) 17
  set _Tree($w:lines) 1
  set _Tree($w:fill) [$w cget -insertbackground]
  set _Tree($w:font) $_Tree(font)
  set _Tree($w:iconbg) ""

  image create bitmap _Tree($w:openbm) \
    -maskdata $_Tree(maskdata) -data $_Tree(openbm)
  image create bitmap _Tree($w:closedbm) \
    -maskdata $_Tree(maskdata) -data $_Tree(closedbm)

  eval TreeConfig $w $args
}

#----- configure tree widget

proc TreeConfig {w args} {
  global _Tree

  set opts {}
  foreach {op val} $args {
    switch -glob -- $op {
      -padx   {set _Tree($w:padx) $val}
      -pady   {set _Tree($w:pady) $val}
      -fore* -
      -fg -
      -fill   {set _Tree($w:fill) $val}
      -font   {set _Tree($w:font) $val}
      -ind*   {set _Tree($w:indent) $val}
      -spac*  {set _Tree($w:space) $val}
      -sort   {set _Tree($w:sort) $val}
      -icon  -
      -imag*  {set _Tree($w:icon) $val}
      -tag*   {set _Tree($w:tags) $val}
      -open   {set _Tree($w:open) $val}
      -lin*   {set _Tree($w:lines) $val}
      -iconbg {set _Tree($w:iconbg) $val}
      default {lappend opts $op $val}
    }
  }

  eval $w config $opts

  set bg [$w cget -background]
  _Tree($w:openbm) configure -foreground $_Tree($w:fill) -background $bg
  _Tree($w:closedbm) configure -foreground $_Tree($w:fill) -background $bg

  Tree:buildwhenidle $w
}

#----- destroy the tree widget

proc TreeDestroy {w} {
  global _Tree
  catch {destroy $w}
  image delete _Tree($w:openbm) _Tree($w:closedbm)
  foreach t [array names _Tree $w:*] {
    unset _Tree($t)
  }
}

#----- insert a new item in the tree

proc TreeInsert {w v args} {
  global _Tree
  if {$v == "/"} {
    set _Tree($w:root) 1
  } else {
    set dir [file dirname $v]
    set n [file tail $v]
    foreach {op arg} $args {
      switch -exact -- $op {
        -dir {
          set dir $arg
          set n $v
          if {$dir == "/"} {
            set v "/$n"
          } else {
            set v "$dir/$n"
          }
        }
      }
    }
    if {![info exists _Tree($w:$dir:open)]} {
      return -code error "parent item \"$dir\" is missing"
    }
    lappend _Tree($w:$dir:children) $n
    if {$_Tree($w:sort) != ""} {
      set _Tree($w:$dir:children) [eval lsort $_Tree($w:sort) \
        {$_Tree($w:$dir:children)}]
    }
  }

  Tree:dfltconfig $w $v
  foreach {op val} $args {
    switch -glob -- $op {
      -fore* -
      -fg    -
      -fill   {set _Tree($w:$v:fill) $val}
      -font   {set _Tree($w:$v:font) $val}
      -icon  -
      -imag*  {set _Tree($w:$v:icon) $val}
      -tag*   {set _Tree($w:$v:tags) $val}
      -text   {set _Tree($w:$v:text) $val}
      -open   {set _Tree($w:$v:open) $val}
      -iconbg {set _Tree($w:$v:iconbg) $val}
    }
  }
  Tree:buildwhenidle $w
}

#----- delete item from the tree

proc TreeDelete {w v} {
  global _Tree
  if {![info exists _Tree($w:$v:open)]} return
  if {$v == "/"} {
    foreach t [array names _Tree $w:*:*] {
      unset _Tree($t)
    }
    set _Tree($w:root) 0
    set _Tree($w:selection) {}
    set _Tree($w:selidx) {}
    set _Tree($w:list) {}
    Tree:buildwhenidle $w
    return
  }
  foreach c $_Tree($w:$v:children) {
    catch {TreeDelete $w $v/$c}
  }
  foreach t [array names _Tree $w:$v:*] {
    unset _Tree($t)
  }
  set dir [file dirname $v]
  set n [file tail $v]
  set i [lsearch -exact $_Tree($w:$dir:children) $n]
  if {$i >= 0} {
    set _Tree($w:$dir:children) [lreplace $_Tree($w:$dir:children) $i $i]
  }
  if {$_Tree($w:selection) != "" &&
    ![TreeExists $w $_Tree($w:selection)]} {
    set _Tree($w:selection) ""
  }
  Tree:buildwhenidle $w
}

#----- move/rename an item

proc TreeMove {w vold vnew} {
  global _Tree
  if {![TreeExists $w $vold]} return
  set pold [file dirname $vold]
  set pnew [file dirname $vnew]
  if {$pnew == "."} {
    if {$vnew == [file tail $vold]} return
  } else {
    if {![TreeExists $w $pnew]} return
  }
  set n [lsearch -exact $_Tree($w:$pold:children) [file tail $vold]]
  if {$pnew == "."} {
    set _Tree($w:$pold:children) \
      [lreplace $_Tree($w:$pold:children) $n $n $vnew]
    if {$pold == "/"} {
      set vnew "/$vnew"
    } else {
      set vnew "$pold/$vnew"
    }
    set pnew $pold
  } else {
    set _Tree($w:$pold:children) \
      [lreplace $_Tree($w:$pold:children) $n $n]
    lappend _Tree($w:$pnew:children) [file tail $vnew]
  }
  if {$_Tree($w:sort) != ""} {
    set _Tree($w:$pnew:children) [eval lsort $_Tree($w:sort) \
      {$_Tree($w:$pnew:children)}]
  }
  set n [string length "$w:$vold"]
  foreach told [array names _Tree "$w:$vold\[/:\]*"] {
    set _Tree($w:$vnew[string range $told $n end]) $_Tree($told)
    unset _Tree($told)
  }
  set _Tree($w:$vnew:text) [file tail $vnew]
  Tree:buildwhenidle $w
}

#----- check if item exists in tree

proc TreeExists {w v} {
  global _Tree
  return [info exists _Tree($w:$v:open)]
}

#----- check if item is visible in tree

proc TreeVisible {w v} {
  global _Tree
  if {$_Tree($w:build) != ""} {
    Tree:build $w
  }
  return [info exists _Tree($w:$v:tag)]
}

#----- set/get item

proc TreeSet {w v args} {
  global _Tree
  if {![TreeExists $w $v]} return
  set redraw 0
  foreach {op val} $args {
    switch -glob -- $op {
      -fore* -
      -fg    -
      -fill   {set tag fill}
      -font   {set tag font}
      -icon  -
      -imag*  {set tag icon}
      -tag*   {set tag tags}
      -open   {set tag open}
      -iconbg {set tag iconbg}
      default {set tag ""}
    }
    if {$tag != ""} {
      if {$val == ""} {set val $_Tree($w:$tag)}
      if {$val != $_Tree($w:$v:$tag)} {
        incr redraw
        set _Tree($w:$v:$tag) $val
      }
    }
  }
  if {$redraw} {
    Tree:buildwhenidle $w
  }
}

proc TreeGet {w v op} {
  global _Tree
  if {[TreeExists $w $v]} {
    switch -glob -- $op {
      -fore* -
      -fg    -
      -fill   {return $_Tree($w:$v:fill)}
      -font   {return $_Tree($w:$v:font)}
      -icon  -
      -imag*  {return $_Tree($w:$v:icon)}
      -tag*   {return $_Tree($w:$v:tags)}
      -open   {return $_Tree($w:$v:open)}
      -child* {return $_Tree($w:$v:children)}
      -iconbg {return $_Tree($w:$v:iconbg)}
    }
  }
  return ""
}

#----- open/close a branch of a tree

proc TreeOpen {w v} {
  global _Tree
  if {[info exists _Tree($w:$v:open)] &&
      !$_Tree($w:$v:open) &&
      [llength $_Tree($w:$v:children)] > 0} {
    set _Tree($w:$v:open) 1
    Tree:build $w
  }
}

proc TreeClose {w v} {
  global _Tree
  if {[info exists _Tree($w:$v:open)] && $_Tree($w:$v:open)} {
    set _Tree($w:$v:open) 0
    Tree:build $w
  }
}

proc TreeToggle {w v} {
  global _Tree
  if {[info exists _Tree($w:$v:open)]} {
    if {$_Tree($w:$v:open)} {
      TreeClose $w $v
    } else {
      TreeOpen $w $v
    }
  }
}

#----- open/close one level below a branch

proc TreeOpenLevel {w v} {
  global _Tree
  if {$v == "/"} {
    set cnt 0
    foreach c $_Tree($w:$v:children) {
      incr cnt [Tree:openlevel $w "/$c"]
    }
  } else {
    set cnt [Tree:openlevel $w $v]
  }
  if $cnt {
    Tree:build $w
  }
}

proc TreeCloseLevel {w v} {
  global _Tree
  if {$v == "/"} {
    set cnt 0
    foreach c $_Tree($w:$v:children) {
      incr cnt [Tree:closelevel $w "/$c"]
    }
  } else {
    set cnt [Tree:closelevel $w $v]
  }
  if $cnt {
    Tree:build $w
  }
}

#----- expand/collapse a branch of the tree

proc TreeExpand {w v} {
  global _Tree
  if {$_Tree($w:root)} {
    Tree:setopen $w $v 1
    Tree:build $w
  }
}

proc TreeCollapse {w v} {
  global _Tree
  if {$_Tree($w:root)} {
    Tree:setopen $w $v 0
    Tree:build $w
  }
}

#----- set/get selection

proc TreeSelectionSet {w v} {
  global _Tree
  if {[TreeExists $w $v]} {
    set _Tree($w:selection) $v
  } else {
    set _Tree($w:selection) ""
  }
  if {$_Tree($w:build) == ""} {
    Tree:drawselection $w
  }
}

proc TreeSelectionGet {w} {
  global _Tree
  return $_Tree($w:selection)
}

#----- get item at specified coordinates

proc TreeAt {w x y} {
  global _Tree
  set x [$w canvasx $x]
  set y [$w canvasy $y]
  foreach m [$w find overlapping $x $y $x $y] {
    if {[info exists _Tree($w:tag:$m)]} {
      return $_Tree($w:tag:$m)
    }
  }
  return ""
}

proc TreeTypeAt {w x y} {
  global _Tree
  set x [$w canvasx $x]
  set y [$w canvasy $y]
  foreach m [$w find overlapping $x $y $x $y] {
    if {[info exists _Tree($w:tag:$m)]} {
      return [list [$w type $m] $_Tree($w:tag:$m)]
    }
  }
  return ""
}

#----- search for an item

proc TreeFind {w v pat {case 0}} {
  global _Tree
  if {![TreeExists $w $v]} {return ""}
  if {$case} {
    set opts ""
  } else {
    set opts -nocase
  }
  if {$v == "/"} {
    return [Tree:match $w / 0 $pat $opts]
  }
  if {$_Tree($w:$v:children) != {}} {
    set node [Tree:match $w $v 0 $pat $opts]
    if {$node != ""} {return $node}
  }
  while {$v != "/"} {
    set parent [file dirname $v]
    set node [file tail $v]
    set n [lsearch -exact $_Tree($w:$parent:children) $node]
    if {$n >= 0} {
      incr n
      set node [Tree:match $w $parent $n $pat $opts]
      if {$node != ""} {return $node}
    }
    set v $parent
  }
  return ""
}

#----- force item to be visible

proc TreeSee {w v} {
  global _Tree
  if {![info exists _Tree($w:$v:open)]} return
  set redraw 0
  set node {}
  foreach p [split [file dirname $v] /] {
    if {$p != {}} {
      append node "/$p"
      if {!$_Tree($w:$node:open)} {
        set _Tree($w:$node:open) 1
        incr redraw
      }
    }
  }
  if {$redraw || $_Tree($w:build) != ""} {
    Tree:build $w
  }

  set bb [$w cget -scrollregion]
  set wd [expr [lindex $bb 2] - [lindex $bb 0]]
  set ht [expr [lindex $bb 3] - [lindex $bb 1]]

  set bb [$w bbox $_Tree($w:$v:tag)]
  set yt [$w canvasy 0]
  set yb [$w canvasy [winfo height $w]]
  if {[expr [lindex $bb 1] <= $yt || [lindex $bb 3] >= $yb]} {
    set sr [$w yview]
    set ds [expr [lindex $sr 1] - [lindex $sr 0]]
    set yp [expr 0.5 * (double([lindex $bb 1] + [lindex $bb 3]) / \
      double($ht) - $ds)]
    if {$yp < 0.0} {set yp 0.0}
    $w yview moveto $yp
  }
  set xl [$w canvasx 0]
  set xr [$w canvasx [winfo width $w]]
  if {[expr [lindex $bb 0] <= $xl || [lindex $bb 2] >= $xr]} {
    set sr [$w xview]
    set ds [expr [lindex $sr 1] - [lindex $sr 0]]
    set xp [expr 0.5 * (double([lindex $bb 0] + [lindex $bb 2]) / \
      double($wd) - $ds)]
    if {$xp < 0.0} {set xp 0.0}
    $w xview moveto $xp
  }
}

#----- edit item

proc TreeEdit {w v} {
  global _Tree
  if {![info exists _Tree($w:$v:open)]} return
  TreeSee $w $v
  set tag $_Tree($w:$v:tag)
  set coords [$w coords $tag]
  set x [lindex $coords 0]
  set y [lindex $coords 1]

  set fg $_Tree($w:$v:fill)
  set bg [$w cget -bg]
  set wd [expr {[winfo width $w] - [$w cget -bd] - \
    [$w cget -highlightthickness]}]
  set wmax [expr {[$w canvasx $wd] - $x}]

  set _Tree($w:edit) [file tail $v]
  set selection [TreeSelectionGet $w]
  TreeSelectionSet $w ""
  $w itemconfigure $tag -fill $bg

  set frame [frame $w.edit -relief flat -bd 0 \
    -highlightthickness 0 -bg $bg]
  set ent [entry $frame.ent -relief solid -bd 0 \
    -highlightthickness 1 -fg $fg -bg $bg -width 0 \
    -font [$w itemcget $tag -font] -textvariable _Tree($w:edit)]
  pack $ent -side left -anchor w -ipadx 5
  bind $ent <Escape> "set _Tree($w:done) 0"
  bind $ent <Return> "set _Tree($w:done) 1"
  set id [$w create window $x $y -window $frame -anchor w]

  trace variable _Tree($w:edit) w "Tree:edit_size $w $ent $id $wmax"
  set oldFocus [focus]
  set oldGrab [grab current $w]
  if {$oldGrab != ""} {
    set grabStatus [grab status $oldGrab]
  }
  catch {grab $frame}
  tkwait visibility $ent
  focus $ent
  $ent selection range 0 end
  $ent icursor end
  $ent xview end
  tkwait variable _Tree($w:done)
  trace vdelete _Tree($w:edit) w "Tree:edit_size $w $ent $id $wmax"

  catch {focus $oldFocus}
  if {$oldGrab != ""} {
    if {$grabStatus == "global"} {
      grab -global $oldGrab
    } else {
      grab $oldGrab
    }
  }
  destroy $frame
  $w delete $id
  $w itemconfigure $tag -fill $fg
  TreeSelectionSet $w $selection
  if {$_Tree($w:done)} {
    return $_Tree($w:edit)
  }
  return ""
}

#----- get next/previous visible item

proc TreeNext {w v} {
  global _Tree
  if {![TreeVisible $w $v]} {return ""}
  set n [lsearch -exact $_Tree($w:list) $_Tree($w:$v:tag)]
  if {$n < 0} {return ""}
  incr n
  if {$n == [llength $_Tree($w:list)]} {set n 0}
  return $_Tree($w:tag:[lindex $_Tree($w:list) $n])
}

proc TreePrev {w v} {
  global _Tree
  if {![TreeVisible $w $v]} {return ""}
  set n [lsearch -exact $_Tree($w:list) $_Tree($w:$v:tag)]
  if {$n < 0} {return ""}
  incr n -1
  if {$n < 0} {
    set n [expr [llength $_Tree($w:list)] - 1]
  }
  return $_Tree($w:tag:[lindex $_Tree($w:list) $n])
}

#----------------------------------------------------------------
# these are for internal use only
#----------------------------------------------------------------

proc Tree:dfltconfig {w v} {
  global _Tree
  set _Tree($w:$v:children) {}
  if {![info exists _Tree($w:$v:open)]} {
    set _Tree($w:$v:open) $_Tree($w:open)
  }
  foreach i {fill font icon tags iconbg} {
    set _Tree($w:$v:$i) $_Tree($w:$i)
  }
  if {$v == "/"} {
    set _Tree($w:$v:text) $v
  } else {
    set _Tree($w:$v:text) [file tail $v]
  }
}

proc Tree:buildwhenidle {w} {
  global _Tree
  if {$_Tree($w:build) == ""} {
    set _Tree($w:build) [after idle "Tree:build $w"]
  }
}

proc Tree:build {w} {
  global _Tree
  if {$_Tree($w:build) != ""} {
    catch {after cancel $_Tree($w:build)}
    set _Tree($w:build) ""
  }
  $w delete all
  foreach t [array names _Tree $w:*:tag] {
    unset _Tree($t)
  }
  set _Tree($w:list) {}
  if {!$_Tree($w:root)} return
  set icon $_Tree($w:/:icon)
  set taglist x
  foreach tag $_Tree($w:/:tags) {
    lappend taglist $tag
  }
  set x 5
  set y 5
  if {$icon != ""} {
    set k [$w create image $x $y -image $icon -anchor w -tags $taglist]
    incr x 20
    set _Tree($w:tag:$k) /
  }
  set j [$w create text $x $y -text $_Tree($w:/:text) \
    -font $_Tree($w:/:font) -fill $_Tree($w:/:fill) \
    -anchor w -tags $taglist]
  lappend _Tree($w:list) $j
  set _Tree($w:tag:$j) /
  set _Tree($w:/:tag) $j
  set _Tree($w:y) [expr $y + 17]
  Tree:buildlayer $w / $_Tree($w:indent)
  set bb [$w bbox all]
  $w config -scrollregion [list [expr [lindex $bb 0] - 5] \
    [expr [lindex $bb 1] - 5] [expr [lindex $bb 2] + 5] \
    [expr [lindex $bb 3] + 5]]
  Tree:drawselection $w
}

proc Tree:buildlayer {w v in} {
  global _Tree
  if {$v=="/"} {
    set vx {}
  } else {
    set vx $v
  }
  set start [expr $_Tree($w:y)-10]
  set y $start
  foreach c $_Tree($w:$v:children) {
    set y $_Tree($w:y)
    incr _Tree($w:y) 17
    if {$_Tree($w:lines)} {
      $w create line $in $y [expr $in+10] $y -fill $_Tree($w:fill)
    }
    set icon $_Tree($w:$vx/$c:icon)
    set taglist x
    foreach tag $_Tree($w:$vx/$c:tags) {
      lappend taglist $tag
    }
    set x [expr $in+12]
    if {$icon != ""} {
      set k [$w create image $x $y -image $icon -anchor w -tags $taglist]
      incr x [image width $icon]
      if {$_Tree($w:$vx/$c:iconbg) != ""} {
        set bb [$w bbox $k]
        set r [$w create rectangle $bb -fill $_Tree($w:$vx/$c:iconbg) -outline ""]
        $w lower $r
      }
      set _Tree($w:tag:$k) $vx/$c
    }
    incr x $_Tree($w:padx)
    set j [$w create text $x $y -text $_Tree($w:$vx/$c:text) \
      -font $_Tree($w:$vx/$c:font) -fill $_Tree($w:$vx/$c:fill) \
      -anchor w -tags $taglist]
    lappend _Tree($w:list) $j
    set _Tree($w:tag:$j) $vx/$c
    set _Tree($w:$vx/$c:tag) $j
    if {[string length $_Tree($w:$vx/$c:children)]} {
      if {$_Tree($w:$vx/$c:open)} {
         set j [$w create image $in $y -image _Tree($w:openbm)]
         Tree:buildlayer $w $vx/$c [expr $in+$_Tree($w:indent)]
      } else {
         set j [$w create image $in $y -image _Tree($w:closedbm)]
      }
      $w bind $j <1> "TreeToggle $w {$vx/$c}"
    } else {
      set _Tree($w:$vx/$c:open) 0
    }
  }
  if {$_Tree($w:lines)} {
    set j [$w create line $in $start $in [expr $y+1] -fill $_Tree($w:fill)]
    $w lower $j
  }
}

proc Tree:setopen {w v open} {
  global _Tree
  if {[info exists _Tree($w:$v:children)] &&
      [string length $_Tree($w:$v:children)]>0} {
    set _Tree($w:$v:open) $open
    if {$v=="/"} {
      set vx {}
    } else {
      set vx $v
    }
    foreach c $_Tree($w:$v:children) {
      Tree:setopen $w $vx/$c $open
    }
  }
}

proc Tree:openlevel {w v} {
  global _Tree
  set cnt 0
  if {[info exists _Tree($w:$v:children)] &&
      [llength $_Tree($w:$v:children)] > 0} {
    if {!$_Tree($w:$v:open)} {
      set _Tree($w:$v:open) 1
      incr cnt
    } else {
      foreach c $_Tree($w:$v:children) {
        incr cnt [Tree:openlevel $w "$v/$c"]
      }
    }
  }
  return $cnt
}

proc Tree:closelevel {w v} {
  global _Tree
  set cnt 0
  if {[info exists _Tree($w:$v:children)] &&
      [llength $_Tree($w:$v:children)] > 0 &&
      $_Tree($w:$v:open)} {
    foreach c $_Tree($w:$v:children) {
      if {$_Tree($w:$v/$c:open)} {
        incr cnt [Tree:closelevel $w "$v/$c"]
      }
    }
    if !$cnt {
      set _Tree($w:$v:open) 0
      incr cnt
    }
  }
  return $cnt
}

proc Tree:drawselection {w} {
  global _Tree
  set selfg [$w cget -selectforeground]
  set selbg [$w cget -selectbackground]
  if {[string length $_Tree($w:selidx)]} {
    set id [$w find withtag sel]
    if {$id != ""} {
      set tag [lindex [$w gettags $id] 1]
      if {[info exists _Tree($w:tag:$tag)]} {
        set v $_Tree($w:tag:$tag)
        $w itemconfigure $tag -fill $_Tree($w:$v:fill)
      }
    }
    $w delete $_Tree($w:selidx)
  }
  set v $_Tree($w:selection)
  if {![string length $v] ||![info exists _Tree($w:$v:tag)]} return
                                                                                                                                                                                                                                                                                                                                                                                                                                                                         set bbox [$w bbox $_Tree($w:$v:tag)]
  if {[llength $bbox]==4} {
    set tag $_Tree($w:$v:tag)
    $w itemconfigure $tag -fill $selfg
    set i [eval $w create rectangle $bbox -fill $selbg \
      -outline $selbg -tags [list "sel $tag"]]
    set _Tree($w:selidx) $i
    $w lower $i
  } else {
    set _Tree($w:selidx) {}
  }
}

proc Tree:match {w v n pat args} {
  global _Tree
  if {$v == "/"} {
    set vx ""
  } else {
    set vx $v
  }
  set len [llength $_Tree($w:$v:children)]
  while {$n < $len} {
    set c [lindex $_Tree($w:$v:children) $n]
    if {[eval string match $args {$pat} {$c}] ||
        [eval string match $args {$pat} {$vx/$c}]} {return "$vx/$c"}
    set node [Tree:match $w "$vx/$c" 0 $pat $args]
    if {$node != ""} {return $node}
    incr n
  }
  return ""
}

proc Tree:edit_size {path entry id wmax args} {
  set entw [winfo reqwidth $entry]
  if {$entw + 5 >= $wmax} {
    $path itemconfigure $id -width $wmax
  } else {
    $path itemconfigure $id -width 0
  }
}

