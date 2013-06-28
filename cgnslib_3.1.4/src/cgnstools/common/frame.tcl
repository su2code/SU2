
#----- Create a new title frame widget

proc FrameCreate {w args} {
  global _Frame

  button $w
  set _Frame($w,bg) [$w cget -bg]
  set _Frame($w,fg) [$w cget -fg]
  set _Frame($w,font) [$w cget -font]
  set _Frame($w,text) ""
  set _Frame($w,side) left
  set _Frame($w,base) center
  set _Frame($w,padx) 5
  set _Frame($w,pady) 5
  set _Frame($w,ipad) 2
  set _Frame($w,bd) 2
  set _Frame($w,indent) 10
  set _Frame($w,relief) groove
  destroy $w

  frame $w -bd 0 -highlightthickness 0 -relief flat
  frame $w.p -bd 0 -highlightthickness 0 -relief flat
  frame $w.b -bd 0 -highlightthickness 0 -relief flat
  pack $w.b -side bottom -fill both -expand 1
  frame $w.b.p -bd 0 -highlightthickness 0 -relief flat
  label $w.l -bd 0 -highlightthickness 0 -relief flat -padx 0 -pady 0
  frame $w.f

  eval FrameConfig $w $args
  return $w
}

#----- Change configuration options for the notebook widget

proc FrameConfig {w args} {
  global _Frame

  # get configuration options

  foreach {tag value} $args {
    switch -- $tag {
      -text {
        set _Frame($w,text) $value
      }
      -padx {
        set _Frame($w,padx) $value
      }
      -pady {
        set _Frame($w,pady) $value
      }
      -ipad {
        set _Frame($w,ipad) $value
      }
      -bg {
        set _Frame($w,bg) $value
      }
      -fg {
        set _Frame($w,fg) $value
      }
      -font {
        set _Frame($w,font) $value
      }
      -side {
        set _Frame($w,side) $value
      }
      -base {
        set _Frame($w,base) $value
      }
      -bd {
        set _Frame($w,bd) $value
      }
      -relief {
        set _Frame($w,relief) $value
      }
      -indent {
        set _Frame($w,indent) $value
      }
    }
  }

  catch {place forget $w.l}
  catch {pack forget $w.f}

  $w.b configure -bd $_Frame($w,bd) -relief $_Frame($w,relief)
  pack $w.f -in $w.b -side bottom -fill both -padx $_Frame($w,padx) \
    -pady $_Frame($w,pady) -expand 1

  if {$_Frame($w,text) == ""} {
    catch {pack forget $w.p}
    catch {pack forget $w.b.p}
    return
  }

  catch {pack $w.p -side top -fill x}
  catch {pack $w.b.p -side top -fill x}

  $w.l configure -text $_Frame($w,text) -font $_Frame($w,font) \
    -fg $_Frame($w,fg) -padx $_Frame($w,ipad)

  set height [winfo reqheight $w.l]
  switch $_Frame($w,base) {
    top {set htop $height; set hbot 1; set y 0}
    center {set htop [expr {$height / 2}]
            set hbot [expr {$height - $htop}]
            set y 0}
    bottom {set htop 1; set hbot $height
            set y [expr {$_Frame($w,bd) + 1}]}
  }
  $w.p configure -height $htop
  $w.b.p configure -height $hbot

  switch $_Frame($w,side) {
    left {set relx 0.0
          set x $_Frame($w,indent)
          set anchor nw}
    center {set relx 0.5;
            set x 0
            set anchor n}
    right {set relx 1.0
           set x -$_Frame($w,indent)
           set anchor ne}
  }
  place $w.l -relx $relx -x $x -y $y -anchor $anchor
}

proc FrameGet {w} {
  return $w.f
}

proc FrameLabel {w} {
  return $w.l
}

