#!/bin/sh
# the next line restarts using wish \
exec plotwish -f "$0" "$@"

proc error_exit {msg} {
  wm withdraw .
  tk_dialog .error Error $msg error 0 Exit
  exit 1
}

if {[catch {package require Tk 8.0} msg]} {
  error_exit $msg
}

#---------- get platform

set platform $tcl_platform(platform)
if {$platform == "windows" && [info exists env(TERM)] &&
    ($env(TERM) == "cygwin" || $env(TERM) == "xterm")} {
  set platform cygwin
}

#----- get startup directory and name

set cmd_name [file tail $argv0]
set cmd_dir  [file dirname $argv0]
if {![file exists $argv0] || [file isdirectory $argv0]} {
  if {$platform == "windows"} {
    set sep ";"
  } else {
    set sep ":"
  }
  foreach i [split $env(PATH) $sep] {
    if {$sep == ";"} {
      set i [join [split $i \\] /]
    }
    if {[file exists $i/$cmd_name] && ![file isdirectory $i/$cmd_name]} {
      set cmd_dir $i
      break;
    }
  }
}
set curdir [pwd]
if ![catch {cd $cmd_dir}] {
  set cmd_dir [pwd]
  cd $curdir
}
if {$tcl_platform(platform) == "windows"} {
  set cmd_dir [file attributes $cmd_dir -shortname]
}

#----- set path to tcl scripts

set auto_path "$cmd_dir $cmd_dir/../common $auto_path"
if {[info exists env(TCL_PROC_DIR)]} {
  set auto_path "$env(TCL_PROC_DIR) $auto_path"
}

#----- initialize

if {$platform == "windows"} {
  set vers [join [split $tcl_version .] {}]
  if {[info commands CGNSopen] == {}} {
    if {[catch {load cgnstcl$vers} msg]} {
      error_exit $msg
    }
  }
  catch {load tclreg$vers registry}
} else {
  if {[info commands CGNSopen] == {}} {
    error_exit "need to run script with plotwish"
  }
}

if {[catch {package require Tkogl} msg]} {
  error_exit $msg
}

if [catch {config_defaults 1} msg] {error_exit $msg}

array set ProgData {
  winwidth 700
  winheight 500
  seppos 0.25
  sepwd 7
  cgnsfile ""
  displaylist ""
  axislist ""
  planelist ""
  cutlist ""
  background {0.0 0.0 0.0}
  autocolor 1
  showcolors 1
  showlines 1
  twosided 0
  culling disable
  revnorm 0
  axis 0
  meshvis 1
  meshmode 1
  regvis 0
  regmode 2
  edgeangle 60
  fov 30
  np 0.025
  fp 2.0
  bases {}
  nzones 0
  zones {}
  curbase ""
  curnode ""
  curdim 0
  curclr ""
  curmode 0
  cutplane ""
  cutcolor {0.7 0.7 0.4 0.5}
  cutmode 1
  usecutclr 0
  ignorevis 0
  drawID ""
  dotrace 0
  reg,file ".cgnstools"
  reg,base "HKEY_CURRENT_USER/Software/CGNS"
  reg,key  "CGNSplot"
  reg,vals {cgnsfile background autocolor twosided culling \
            axis edgeangle meshvis meshmode regvis regmode fov \
            winwidth winheight seppos cutcolor usecutclr ignorevis \
            showcolors showlines}
}

#----- read registry

if {[tclreg_init -base $ProgData(reg,base) -fname $ProgData(reg,file)]} {
  foreach i $ProgData(reg,vals) {
    if {![catch {tclreg_get $ProgData(reg,key) $i} val] && $val != ""} {
      set ProgData($i) $val
    }
  }
}

#---------- main window

wm title . "CGNSplot"
wm protocol . WM_DELETE_WINDOW do_quit

proc do_quit {} {
  global ProgData
  catch CGNSclose
  set ProgData(winwidth) [winfo width .main]
  set ProgData(winheight) [winfo height .main]
  foreach i $ProgData(reg,vals) {
    catch {tclreg_set $ProgData(reg,key) $i $ProgData($i)}
  }
  catch tclreg_close
  catch {WinHtml close}
  destroy .
  exit 0
}

#----- menu

menubar_create {File Options View Display Help}

#----- file menu

set m [menubar_get File]
$m add command -label "Load CGNS..." -command load_cgns
$m add separator
$m add command -label "Quit" -command do_quit

#----- options menu

set m [menubar_get Options]
$m add command -label "Set Defaults..." -command set_defaults
$m add command -label "Set Perspective..." -command set_perspective

#----- view menu

set m [menubar_get View]
$m add command -label "Open Level" -command tree_open
$m add command -label "Close Level" -command tree_close
$m add separator
$m add command -label "Reset View" -command {reset_view; $OGLwin redraw}
$m add command -label "Center View" -command Center

#----- display menu

set m [menubar_get Display]
$m add checkbutton -label "Two Sided" -command set_twosided \
  -variable ProgData(twosided) -onvalue 1 -offvalue 0
$m add checkbutton -label "Backface Culling" -command set_culling \
  -variable ProgData(culling) -onvalue enable -offvalue disable
#$m add checkbutton -label "Reverse Normals" -command set_normals \
#  -variable ProgData(revnorm) -onvalue 1 -offvalue 0
$m add checkbutton -label "Display Axis" -command set_axis \
  -variable ProgData(axis) -onvalue 1 -offvalue 0

if {[info commands OGLcutplane] != ""} {
  $m add command -label "Cutting Plane..." -command cutplane_control
}

proc set_twosided {} {
  global ProgData OGLwin
  if {$ProgData(twosided)} {
    .toolbar.but.b3.twosided configure -relief sunken
    set_balloon .toolbar.but.b3.twosided "Disable Two-Sided Lighting"
  } else {
    .toolbar.but.b3.twosided configure -relief raised
    set_balloon .toolbar.but.b3.twosided "Enable Two-Sided Lighting"
  }
  $OGLwin eval -lightmodel lightmodeltwoside $ProgData(twosided)
  $OGLwin redraw
}

proc set_culling {} {
  global ProgData OGLwin
  if {$ProgData(culling) == "enable"} {
    .toolbar.but.b3.culling configure -relief sunken
    set_balloon .toolbar.but.b3.culling "Disable Backface Culling"
  } else {
    .toolbar.but.b3.culling configure -relief raised
    set_balloon .toolbar.but.b3.culling "Enable Backface Culling"
  }
  $OGLwin eval -$ProgData(culling) cullface
  $OGLwin redraw
}

proc set_normals {} {
  global ProgData OGLwin
  if {$ProgData(revnorm)} {
    .toolbar.but.b3.revnorm configure -relief sunken
    set_balloon .toolbar.but.b3.revnorm "Default Face Normals"
  } else {
    .toolbar.but.b3.revnorm configure -relief raised
    set_balloon .toolbar.but.b3.revnorm "Reverse Face Normals"
  }
}

proc set_axis {} {
  global ProgData OGLwin
  if {$ProgData(axis)} {
    .toolbar.but.b3.axis configure -relief sunken
    set_balloon .toolbar.but.b3.axis "Disable Axis Display"
  } else {
    .toolbar.but.b3.axis configure -relief raised
    set_balloon .toolbar.but.b3.axis "Enable Axis Display"
  }
  OGLaxis $ProgData(axis)
  $OGLwin redraw
}

#----- help menu

set m [menubar_get Help]
$m add command -label "CGNSplot Help..." \
  -command {help_show cgnsplot "" cgnstools/cgnsplot/index.html}
$m add command -label "CGNS Help..." -command {help_show cgns}
$m add separator
$m add command -label "Configure Help..." -command help_setup
$m add separator
$m add command -label "About..." -underline 0 -command do_about

image create photo img_about -data {\
R0lGODlhIAAgALMAAAAA/wD//4QAAMbGxv8AAP//////////////////////////////////\
/////////ywAAAAAIAAgAAAEv3AMIGelEttbN8+fp3mkGGpj2plgawluiZ6XAJMooO98v9sC\
H881cwkINtyqyBocgaoOs3hE2qI0oqQKvMWas6o1eRKaxd2gT0sUj8my7MSdvo3AFfebfJ/m\
CYCAdTA5Zj+BgYM9AQEAjI6MjYx6e1CRko+RkJSVNpiYkJmTiIJ1mp+hnwE2iHWphmeJXUKi\
p5eSXEC2l6mqorO3r7e+tQK1tb3ExMmOsM7Ow8rSydHU09HD1tfYx9nYvc/hQhEAADs=}

proc do_about {} {
  global ProgData
  dialog .about -1 -1 "About CGNSplot" \
"CGNSplot Version 3.2

Bruce Wedan
leavingdust@gmail.com" img_about 0 Close
}

#---------- toolbar

frame .toolbar
pack .toolbar -side top -pady 2 -fill x

set f [frame .toolbar.but]
pack $f -side left

#---

image create photo img_open -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAQwQ4EMlJKwJvZcC7BxdnSV04nCgKjtR6vZgmZ49L\
bmlus7xV9j4QQPYRtWbI3RCXU10WgKaTVPQAexEAOw==}

button $f.open -image img_open -takefocus 0 \
  -command {load_cgns ""}
pack $f.open -side left -padx 5
set_balloon $f.open "Open a CGNS File..."

#---

image create photo img_colors -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAQwRFEMlJ63TYsc2ab1QiJku5KKhCPewDvMAhH+Fo\
nimVadwHTi0XbEazGI/IyohkSqkmy5tTl+F0PlSM1bdqwWKzLutLTB4jADs=}

image create photo img_viewing -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAQwQqEMlJqwT2uo0z6JkFht44fdSnVuipmiXpyUhL\
wiwMbA5onzzcj3bB0ToRADs=}

set b [frame $f.b1]
pack $b -side left -padx 5

button $b.colors -image img_colors -takefocus 0 \
  -command set_defaults
set_balloon $b.colors "Set Defaults..."

button $b.viewing -image img_viewing -takefocus 0 \
  -command set_perspective
set_balloon $b.viewing "Set Perspective..."

pack $b.colors $b.viewing -side left -padx 1

#---

image create photo img_expand -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAAAQvEMlJq712DJzX5pP2VQBIWuWVUitBTCsMlMTz\
IrHl3jnYUy7J7Cep8UA7k3IpiQAAOw==}

image create photo img_collapse -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAAAQkEMlJq7046wrwGBQgXsMCTh35hRvbSqm1urJJ\
V98Jvzvv/5UIADs=}

image create photo img_reset -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAAAQiEMlJKwXYWqA17xXGjQgZguA3qWhbni58vaQY\
vzHb2XdvRQA7}

image create photo img_center -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAQwQmEMmJAL22Ypmp/dolduJYcqiXqiB5umcFwG9s\
YzOd07fZy7yfJAIAOw==}

set b [frame $f.b2]
pack $b -side left -padx 5

button $b.expand -image img_expand -takefocus 0 -command tree_open
set_balloon $b.expand "Open One Level"

button $b.collapse -image img_collapse -takefocus 0 -command tree_close
set_balloon $b.collapse "Close One Level"

button $b.reset -image img_reset -takefocus 0 \
  -command {reset_view; $OGLwin redraw}
set_balloon $b.reset "Reset View"

button $b.center -image img_center -takefocus 0 \
  -command Center
set_balloon $b.center "Center View"

pack $b.expand $b.collapse $b.reset $b.center -side left

#---

image create photo img_twosided -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAQwRGEMlJKwLnUMZoBtSyPNwjWpJIMuaCbh2KTQvn\
Sgd4aak95S/OyyIc0nyyG6KIWOgwB9GoJP1ApavW4iMBKJnOV42hNDIpEQA7}

image create photo img_culling -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAAAQqEMlJq5UgX6oxuF0VhuCHfKR3YlM3spzZwpxI\
z/Vq56nc+zcRKlXaGC0RADs=}

image create photo img_axis -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAQwQwEEkJUJ04o5V6p9o2hJlAnteJLaOKKIqJujNt\
Z+ydY1UqwYKAMGAp/mIqgNKm9J0iADs=}

image create photo img_revnorm -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAAAQ5EMlJq6XsIsAn+1rnhdsEMMCVStzpiCwLODRD\
v2s634595zIeDVYK4iw600q1RGSYTc3l85SOrNgIADs=}

image create photo img_cutting -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAQwQwEMlJq5UEkFsBQB6HbGNpXp+UimDIErBJWjFL\
rSuH22DP+7wNbFYhCTEXI4b4m0QAADs=}

set b [frame $f.b3]
pack $b -side left -padx 5

button $b.twosided -image img_twosided -takefocus 0 \
  -command toggle_twosided
if {$ProgData(twosided)} {
  $b.twosided configure -relief sunken
  set_balloon $b.twosided "Disable Two-Sided Lighting"
} else {
  set_balloon $b.twosided "Enable Two-Sided Lighting"
}

button $b.culling -image img_culling -takefocus 0 \
  -command toggle_culling
if {$ProgData(culling) == "enable"} {
  $b.culling configure -relief sunken
  set_balloon $b.culling "Disable Backface Culling"
} else {
  set_balloon $b.culling "Enable Backface Culling"
}

button $b.revnorm -image img_revnorm -takefocus 0 \
  -command toggle_normals
if {$ProgData(revnorm)} {
  $b.revnorm configure -relief sunken
  set_balloon $b.revnorm "Default Face Normals"
} else {
  set_balloon $b.revnorm "Reverse Face Normals"
}

button $b.axis -image img_axis -takefocus 0 \
  -command toggle_axis
if {$ProgData(axis)} {
  $b.axis configure -relief sunken
  set_balloon $b.axis "Disable Axis Display"
} else {
  set_balloon $b.axis "Enable Axis Display"
}

#pack $b.twosided $b.culling $b.revnorm $b.axis -side left -padx 1
pack $b.twosided $b.culling $b.axis -side left -padx 1

if {[info commands OGLcutplane] != ""} {
  button $b.cutting -image img_cutting -takefocus 0 \
    -command cutplane_control
  set_balloon $b.cutting "Cutting Plane..."
  pack $b.cutting -side left -padx 1
}

proc toggle_twosided {} {
  global ProgData
  if {$ProgData(twosided)} {
    set ProgData(twosided) 0
  } else {
    set ProgData(twosided) 1
  }
  set_twosided
}

proc toggle_culling {} {
  global ProgData
  if {$ProgData(culling) == "enable"} {
    set ProgData(culling) disable
  } else {
    set ProgData(culling) enable
  }
  set_culling
}

proc toggle_normals {} {
  global ProgData
  if {$ProgData(revnorm)} {
    set ProgData(revnorm) 0
  } else {
    set ProgData(revnorm) 1
  }
  set_normals
}

proc toggle_axis {} {
  global ProgData
  if {$ProgData(axis)} {
    set ProgData(axis) 0
  } else {
    set ProgData(axis) 1
  }
  set_axis
}

#---

image create photo img_help -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAQwQiEMlJq50kX5kJ1hvShd+4mSJ4qmTrXl28ehw7\
t+j75joVAQA7}

button $f.help -image img_help -takefocus 0 \
  -command {help_show cgnsplot "" cgnstools/cgnsplot/index.html}
set_balloon $f.help "CGNSplot Help..."

pack $f.help -side left -padx 5

frame .toolsep -bd 1 -height 2 -relief sunken
pack .toolsep -side top -fill x

#----- status bar

set ProgData(status) ""
label .status -textvariable ProgData(status) -relief sunken -anchor w
pack .status -side bottom -fill x -padx 5 -pady 3

proc display_message {msg} {
  global ProgData
  set ProgData(status) $msg
  update idletasks
}

#---------- main window

frame .main -width $ProgData(winwidth) -height $ProgData(winheight)
pack .main -side top -fill both -expand 1 -padx 5 -pady 2

#--- window seperator

frame .main.sep -width 4 -bd 2 -relief raised -cursor sb_h_double_arrow

bind .main.sep <ButtonPress-1> sep_begin_move

proc sep_begin_move {} {
  set width [winfo width .main]
  set x [winfo rootx .main.sep]
  set y [winfo rooty .main.sep]
  set h [winfo height .main.sep]
  set xmin [expr [winfo rootx .main] + 100]
  set xmax [expr [winfo rootx .main] + $width - 100]

  set top [toplevel .main.move -borderwidth 1 \
    -relief raised -cursor sb_h_double_arrow]
  wm overrideredirect $top 1
  wm geom $top "4x$h+$x+$y"

  update idletasks
  grab set $top
  bind $top <ButtonRelease-1> "sep_end_move $top $xmin $xmax %X"
  bind $top <Motion> "sep_move $top $xmin $xmax %X $y"
}

proc sep_move {top xmin xmax x y} {
  if {$x < $xmin} {
    set x $xmin
  } elseif {$x > $xmax} {
    set x $xmax
  }
  wm geom $top "+$x+$y"
}

proc sep_end_move {top xmin xmax x} {
  global ProgData
  destroy $top
  if {$x < $xmin} {
    set x $xmin
  } elseif {$x > $xmax} {
    set x $xmax
  }
  set s [expr double($x - [winfo rootx .main]) / \
    double([winfo width .main])]
  sep_locate $s
}

proc sep_locate {s} {
  global ProgData
  if [expr $s < 0.5] {
    place .main.sel -relx 0 -x 0 \
      -relwidth $s -width -$ProgData(sepwd) \
      -rely 0 -relheight 1
    place .main.plot -relx $s -x $ProgData(sepwd) \
      -relwidth [expr 1.0 - $s] -width -$ProgData(sepwd) \
      -rely 0 -relheight 1
  } else {
    place .main.plot -relx 0 -x 0 \
      -relwidth $s -width -$ProgData(sepwd) \
      -rely 0 -relheight 1
    place .main.sel -relx $s -x $ProgData(sepwd) \
      -relwidth [expr 1.0 - $s] -width -$ProgData(sepwd) \
      -rely 0 -relheight 1
  }
  place .main.sep -relx $s -x -2 -rely 0 -relheight 1
  set ProgData(seppos) $s
}

frame .main.plot
frame .main.sel

sep_locate $ProgData(seppos)

#----- OpenGL window

set OGLwin .main.plot.gl
if {[catch {OGLwin $OGLwin} msg]} {
  error_exit $msg
}
pack $OGLwin -side left -fill both -expand 1

#---------- selections

# base selection

set f [frame .main.sel.base]
pack $f -side top -fill x -pady 5
label $f.lab -text Base
pack $f.lab -side left
set ProgData(baselist) [ComboboxCreate $f.cb -width 10 \
  -edit 0 -state disabled -command set_base]
pack $f.cb -side left -fill x -expand 1

# object tree

image create photo box_empty -data {
R0lGODlhDAAMALMAAAAAAMbGxv//////////////////////////////////////////////\
/////////yH5BAEAAAEALAAAAAAMAAwAAAQaEMhJZ7g4y8zD7tgHesB4iSDaqRyrlWYlSxEA\
Ow==}

image create photo box_checked -data {
R0lGODlhDAAMALMAAAAAAMbGxv//////////////////////////////////////////////\
/////////yH5BAEAAAEALAAAAAAMAAwAAAQjEMhJabj4CrmznkGHddJnBuUnokC2siOwpqJY\
vqTcelvlSxEAOw==}

image create photo box_question -data {
R0lGODlhDAAMALMAAIQAAMbGxv8AAP//////////////////////////////////////////\
/////////yH5BAEAAAEALAAAAAAMAAwAAAQjMMggap1TgL0FrsAEXJIWSqNXbiKJoePbqnIn\
i/Gtv7beSxEAOw==}

set f [frame .main.sel.objs]
pack $f -side top -fill both -expand 1

set ProgData(tree) $f.tree

scrollbar $f.ys -orient vertical -command "$ProgData(tree) yview" \
  -takefocus 0 -highlightthickness 0
pack $f.ys -side right -fill y

scrollbar $f.xs -orient horizontal -command "$ProgData(tree) xview" \
  -takefocus 0 -highlightthickness 0
pack $f.xs -side bottom -fill x

TreeCreate $ProgData(tree) -width 150 -height 300 -relief sunken \
  -bd 2 -highlightthickness 1 -yscrollcommand "$f.ys set" -takefocus 1 \
  -xscrollcommand "$f.xs set" -font $Font(normal) -padx 4 \
  -lines $ProgData(showlines)
pack $ProgData(tree) -side left -fill both -expand 1

bind $ProgData(tree) <1> {tree_show %W %x %y}
bind $ProgData(tree) <ButtonRelease-2> {tree_info %W %x %y}
bind $ProgData(tree) <3> {tree_menu %W %x %y}
bind $ProgData(tree) <Double-1> {tree_at %W %x %y Toggle}
bind $ProgData(tree) <Shift-1> {tree_at %W %x %y Expand}
bind $ProgData(tree) <Control-1> {tree_at %W %x %y Collapse}

set ProgData(menu) [menu .nodemenu -tearoff 0]
$ProgData(menu) add radiobutton -label "Outline" \
  -variable ProgData(curmode) -value 1 -command update_mode
$ProgData(menu) add radiobutton -label "Mesh" \
  -variable ProgData(curmode) -value 2 -command update_mode
$ProgData(menu) add radiobutton -label "Shaded" \
  -variable ProgData(curmode) -value 3 -command update_mode
$ProgData(menu) add separator
$ProgData(menu) add command -label "Color..." \
  -command "update_color $ProgData(tree)"
$ProgData(menu) add command -label "AutoColor" \
  -command "auto_color $ProgData(tree)"
$ProgData(menu) add command -label "Info..."

proc tree_show {w x y} {
  set value [TreeTypeAt $w $x $y]
  if {$value == ""} return
  set type [lindex $value 0]
  set node [lindex $value 1]
  if {$type == "image"} {
    toggle_visibility $node
  } else {
    if {$node != [TreeSelectionGet $w]} {
      set_node $node
    }
  }
}

proc tree_info {w x y} {
  set value [TreeTypeAt $w $x $y]
  if {$value == ""} return
  if {[lindex $value 0] == "text"} {
    show_info $w $x $y [lindex $value 1]
  }
}

proc tree_menu {w x y} {
  global ProgData
  set value [TreeTypeAt $w $x $y]
  if {$value == ""} return
  if {[lindex $value 0] == "image"} return
  set node [lindex $value 1]
  if {$node != [TreeSelectionGet $w]} {
    set_node $node
  }
  if {$ProgData(curdim) == 0} {
    set state disabled
  } else {
    set state normal
  }
  foreach n {0 1 2 4 5} {
    $ProgData(menu) entryconfigure $n -state $state
  }
  set cnt [llength [TreeGet $w $node -tag]]
  if {$cnt < 1 || $cnt > 2} {
    $ProgData(menu) entryconfigure 6 -state disabled
  } else {
    $ProgData(menu) entryconfigure 6 -state normal \
      -command "show_info $w $x $y {$node}"
  }
  $ProgData(menu) post [expr [winfo rootx $ProgData(tree)] + $x] \
    [expr [winfo rooty $ProgData(tree)] + $y]
}

proc tree_at {w x y mode} {
  set value [TreeTypeAt $w $x $y]
  if {$value == ""} return
  set type [lindex $value 0]
  set node [lindex $value 1]
  if {$type == "image"} {
    if {$mode == "Expand"} {
      set_visibility $node 1 1
      if {$node != [TreeSelectionGet $w]} {
        set_node $node
      }
    } elseif {$mode == "Collapse"} {
      set_visibility $node 0 1
    } else {
      return
    }
    update_node $node
  } else {
    if {$node != [TreeSelectionGet $w]} {
      set_node $node
    }
    Tree$mode $w $node
  }
}

proc tree_open {} {
  global ProgData
  catch {TreeOpenLevel $ProgData(tree) /}
}

proc tree_close {} {
  global ProgData
  catch {TreeCloseLevel $ProgData(tree) /}
}

# options

frame .main.sel.opts
pack .main.sel.opts -side top -fill x

set f [frame .main.sel.opts.top]
pack $f -side top -fill x

radiobutton $f.o -text "Outline" -variable ProgData(curmode) -value 1 \
   -command update_mode -state disabled
pack $f.o -side left
radiobutton $f.m -text "Mesh" -variable ProgData(curmode) -value 2 \
   -command update_mode -state disabled
pack $f.m -side right

set f [frame .main.sel.opts.bot]
pack $f -side top -fill x

radiobutton $f.s -text "Shaded" -variable ProgData(curmode) -value 3 \
   -command update_mode -state disabled
pack $f.s -side left

set ProgData(clrbut) [button $f.c -relief solid -text Color... \
  -borderwidth 1 -command "update_color $f.c" -state disabled]
pack $f.c -side right -padx 5

proc options_state {dim} {
  if {$dim} {
    set state normal
  } else {
    set state disabled
  }
  foreach i {top.o top.m bot.s} {
    .main.sel.opts.$i configure -state $state
  }
}

#----- help

proc help_menu {} {
  if {[help_valid cgnsplot]} {
    menubar_state Help normal 0
    .toolbar.but.help configure -state normal
  } else {
    menubar_state Help disabled 0
    .toolbar.but.help configure -state disabled
  }
  if {[help_valid cgns]} {
    menubar_state Help normal 1
  } else {
    menubar_state Help disabled 1
  }
}

help_init {cgnsplot CGNSplot} {cgns CGNS}

#----- set default colors and zone/region viewing

proc set_defaults {{loc .}} {
  global ProgData Defaults OGLwin fgColor Font
  foreach i {showcolors showlines meshvis meshmode regvis regmode \
             background autocolor} {
    set Defaults($i) $ProgData($i)
  }

  set w .defaults
  toplevel $w
  wm title $w "Defaults"
  wm transient $w .
  wm protocol $w WM_DELETE_WINDOW {set ProgData(done) 0}

  FrameCreate $w.tree -text "Display Tree" -font $Font(bold) -padx 0 -pady 0
  pack $w.tree -side top -padx 5 -pady 2 -fill x -expand 1
  set tree [FrameGet $w.tree]

  checkbutton $tree.colors -text "Show Colors" -onvalue 1 -offvalue 0 \
    -variable Defaults(showcolors) -highlightthickness 0
  checkbutton $tree.lines -text "Show Lines" -onvalue 1 -offvalue 0 \
    -variable Defaults(showlines) -highlightthickness 0
  pack $tree.colors $tree.lines -side left -expand 1

  FrameCreate $w.zone -text "Volume Mesh" -font $Font(bold) -padx 0 -pady 0
  pack $w.zone -side top -padx 5 -pady 2 -fill x -expand 1
  set zone [FrameGet $w.zone]

  checkbutton $zone.vis -text "Visible as" -onvalue 1 -offvalue 0 \
    -variable Defaults(meshvis) -highlightthickness 0
  pack $zone.vis -side top -anchor w
  set f [frame $zone.mode]
  pack $f -side top -fill x -expand 1
  radiobutton $f.o -text "Outline" -variable Defaults(meshmode) -value 1 \
     -highlightthickness 0
  radiobutton $f.m -text "Mesh" -variable Defaults(meshmode) -value 2 \
     -highlightthickness 0
  radiobutton $f.s -text "Shaded" -variable Defaults(meshmode) -value 3 \
     -highlightthickness 0
  pack $f.o $f.m $f.s -side left -expand 1

  FrameCreate $w.reg -text Regions -font $Font(bold) -padx 0 -pady 0
  pack $w.reg -side top -padx 5 -pady 2 -fill x -expand 1
  set reg [FrameGet $w.reg]

  checkbutton $reg.vis -text "Visible as" -onvalue 1 -offvalue 0 \
    -variable Defaults(regvis) -highlightthickness 0
  pack $reg.vis -side top -anchor w
  set f [frame $reg.mode]
  pack $f -side top -fill x -expand 1
  radiobutton $f.o -text "Outline" -variable Defaults(regmode) -value 1 \
     -highlightthickness 0
  radiobutton $f.m -text "Mesh" -variable Defaults(regmode) -value 2 \
     -highlightthickness 0
  radiobutton $f.s -text "Shaded" -variable Defaults(regmode) -value 3 \
     -highlightthickness 0
  pack $f.o $f.m $f.s -side left -expand 1

  FrameCreate $w.colors -text Colors -font $Font(bold)
  pack $w.colors -side top -padx 5 -pady 2 -fill x -expand 1
  set colors [FrameGet $w.colors]

  set f [frame $colors.bg]
  pack $f -side top -fill x -expand 1
  label $f.lab -text "Background Color"
  pack $f.lab -side left
  button $f.but -text Select -relief solid -bd 1 \
    -command "background_color $f.but"
  color_button $f.but $Defaults(background)
  pack $f.but -side right

  set auto [frame $colors.auto]
  pack $auto -side top -fill x -expand 1
  label $auto.lab -text "Assign Color Based On"
  pack $auto.lab -side top -anchor w
  set f [frame $auto.mode]
  pack $f -side top -fill x -expand 1
  radiobutton $f.z -text "Zone" -variable Defaults(autocolor) -value 0 \
     -highlightthickness 0
  radiobutton $f.r -text "Region" -variable Defaults(autocolor) -value 1 \
     -highlightthickness 0
  pack $f.z $f.r -side left -expand 1

  set f [frame $w.but]
  pack $f -side top -padx 5 -pady 2 -fill x -expand 1
  button $f.default -text Defaults -width 8 -underline 0 \
    -command "default_values $w"
  button $f.accept -text Accept -width 8 -underline 0 \
    -command {set ProgData(done) 1}
  button $f.cancel -text Cancel -width 8 -underline 0 \
    -command {set ProgData(done) 0}
  pack $f.default $f.accept $f.cancel -side left -pady 5

  bind $w <Alt-d> "$f.default flash; default_values $w"
  bind $w <Alt-a> "$f.accept flash; set ProgData(done) 1"
  bind $w <Alt-c> "$f.cancel flash; set ProgData(done) 0"

  center_window $w $loc
  set oldFocus [focus]
  set oldGrab [grab current $w]
  if {$oldGrab != ""} {
    set grabStatus [grab status $oldGrab]
  }
  catch {grab $w}
  tkwait visibility $w
  focus $w
  tkwait variable ProgData(done)
  catch {focus $oldFocus}
  destroy $w
  if {$oldGrab != ""} {
    if {$grabStatus == "global"} {
      grab -global $oldGrab
    } else {
      grab $oldGrab
    }
  }

  if {$ProgData(done)} {
    foreach i {meshvis meshmode regvis regmode background autocolor} {
      set ProgData($i) $Defaults($i)
    }
    eval $OGLwin eval -clearcolor $ProgData(background)
    init_display
    Center
    if {$Defaults(showlines) != $ProgData(showlines)} {
      set ProgData(showlines) $Defaults(showlines)
      TreeConfig $ProgData(tree) -lines $ProgData(showlines)
    }
    if {$Defaults(showcolors) != $ProgData(showcolors)} {
      set ProgData(showcolors) $Defaults(showcolors)
      build_tree
    }
  }
}

proc default_values {w} {
  global Defaults
  array set Defaults {
    showcolors 1
    showlines 1
    meshvis 1
    meshmode 1
    regvis 0
    regmode 2
    autocolor 0
    background {0.0 0.0 0.0}
  }

  set cf [FrameGet $w.colors]
  color_button $cf.bg.but $Defaults(background)
}

proc background_color {but} {
  global Defaults
  set clr [select_color $Defaults(background) "Background Color" $but]
  if {$clr != ""} {
    set Defaults(background) $clr
    color_button $but $clr
  }
}

proc color_button {but clr} {
  global fgColor bgColor
  if {$clr == ""} {
    $but configure -fg $fgColor(normal) -activeforeground $fgColor(active) \
      -bg $bgColor(normal) -activebackground $bgColor(active) \
      -state disabled
  } else {
    if {[color_gray $clr] > 0.5} {
      set fg black
    } else {
      set fg white
    }
    $but configure -fg $fg -activeforeground $fg \
      -bg [color_value $clr] \
      -activebackground [color_value [color_lighten $clr]] \
      -state normal
  }
}

proc select_color {oldclr title {loc .}} {
  return [color_select .color $title $oldclr $loc]
}

#----- set perspective -----

set Perspective(fov) $ProgData(fov)
set Perspective(fd) [expr 0.5 / tan(0.00872665 * $Perspective(fov))]
set Perspective(np) [expr $Perspective(fd) * $ProgData(np)]
set Perspective(fp) [expr $Perspective(fd) * $ProgData(fp)]
set Perspective(ar) 1.0

set Viewpoint(x) 0
set Viewpoint(y) 0
set Viewpoint(z) $Perspective(fd)

proc apply_perspective {what val} {
  global ProgData Perspective OGLwin Viewpoint
  set Perspective(fov) $ProgData(fov)
  set fd [expr 0.5 / tan(0.00872665 * $Perspective(fov))]
  set Viewpoint(z) [expr $fd * $Viewpoint(z) / $Perspective(fd)]
  set Perspective(fd) $fd
  if {[expr $ProgData(fp) <= $ProgData(np)]} {
    if {$what == "np"} {
      set ProgData(fp) [expr 0.01 * (100.0 * $ProgData(np) + 1.0)]
    } else {
      set ProgData(np) [expr $ProgData(fp) - 0.001]
    }
  }
  set Perspective(np) [expr $fd * $ProgData(np)]
  set Perspective(fp) [expr $fd * $ProgData(fp)]
  compute_viewport $OGLwin
  set_view $OGLwin
}

proc set_perspective {{loc .}} {
  global ProgData fgColor Font

  set w .perspective
  if [winfo exists $w] {
      wm deiconify $w
      return
  }

  toplevel $w
  wm title $w "Perspective"

  FrameCreate $w.opts -font $Font(bold)
  pack $w.opts -side top -padx 5 -pady 5
  set opts [FrameGet $w.opts]

  set f [frame $opts.fov]
  pack $f -side top
  label $f.lab -width 18 -anchor w -text "Field of View"
  scale $f.scl -length 200 -from 5 -to 85 -variable ProgData(fov) \
    -orient horizontal -command {apply_perspective fov}
  pack $f.lab $f.scl -side left

  set f [frame $opts.np]
  pack $f -side top
  label $f.lab -width 18 -anchor w -text "Near Clipping Plane"
  scale $f.scl -length 200 -from 0.001 -to 1 -resolution .001 \
    -variable ProgData(np) -orient horizontal -command {apply_perspective np}
  pack $f.lab $f.scl -side left

  set f [frame $opts.fp]
  pack $f -side top
  label $f.lab -width 18 -anchor w -text "Far Clipping Plane"
  scale $f.scl -length 200 -from .01 -to 5 -resolution .01 \
    -variable ProgData(fp) -orient horizontal -command {apply_perspective fp}
  pack $f.lab $f.scl -side left

  set f [frame $w.but]
  pack $f -side top -fill x -expand 1
  button $f.default -text Defaults -width 8 -command perspective_defaults
  button $f.close -text Close -width 8 -command "destroy $w"
  pack $f.default $f.close -side left -pady 2 -expand 1

  center_window $w $loc
}

proc perspective_defaults {} {
  global ProgData

  array set ProgData {
    fov 30
    np 0.025
    fp 2
  }
  apply_perspective all 0
}

#----- OGL window controls

set ProgData(axislist) [OGLaxis $ProgData(axis)]
$OGLwin main -clear colorbuffer depthbuffer -call $ProgData(axislist)

if {[info commands OGLcutplane] != ""} {
  OGLcutconfig $ProgData(cutcolor) $ProgData(usecutclr) $ProgData(ignorevis)
  set ProgData(cutlist) [OGLcutplane]
  set ProgData(planelist) [OGLdrawplane]
}

eval $OGLwin eval \
  -clearcolor $ProgData(background) \
  -$ProgData(culling) cullface \
  -matrixmode projection \
  -loadidentity \
  -perspective $Perspective(fov) $Perspective(ar) \
    $Perspective(np) $Perspective(fp) \
  -matrixmode modelview \
  -loadidentity \
  -light light0 position 0 1 5 0 \
  -light light0 diffuse .8 .8 .8 \
  -light light0 specular 1 1 1 \
  -light light0 ambient .2 .2 .2 \
  -lightmodel lightmodeltwoside $ProgData(twosided) \
  -enable lighting -enable light0 \
  -enable depthtest -enable normalize -shademodel flat

set ViewMatrix [$OGLwin get modelviewmatrix]

bind $OGLwin <Configure> {compute_viewport %W}
bind $OGLwin <Any-ButtonPress> {set xlast %x ; set ylast %y}
bind $OGLwin <B1-Motion> {Rotate %W %x %y}
bind $OGLwin <Shift-B1-Motion> {Pan %W %x %y}
bind $OGLwin <Control-B1-Motion> {Zoom %W %x %y}
bind $OGLwin <B2-Motion> {Zoom %W %x %y}
bind $OGLwin <B3-Motion> {Pan %W %x %y}

catch {bind . <MouseWheel> {WheelZoom %D}}
bind . <KeyPress-c> {Center}
bind . <KeyPress-r> {reset_view ; $OGLwin redraw}
bind . <KeyPress-t> toggle_twosided
bind . <KeyPress-b> toggle_culling
bind . <KeyPress-a> toggle_axis
if {[info commands OGLcutplane] != ""} {
  bind . <KeyPress-p> cutplane_control
  bind . <Return> {if [winfo exists .cutplane] draw_cut}
}

proc compute_viewport {w} {
  global ProgData Perspective
  set ProgData(width) [winfo width $w]
  set ProgData(height) [winfo height $w]
  set ProgData(scale) [expr 1.0 / double($ProgData(height))]
  set Perspective(ar) [expr double($ProgData(width)) / \
      double($ProgData(height))]
  $w eval \
    -matrixmode projection \
    -loadidentity \
    -perspective $Perspective(fov) $Perspective(ar) \
      $Perspective(np) $Perspective(fp) \
    -matrixmode modelview
}

proc set_view {w} {
  global Viewpoint ViewMatrix
  eval $w eval \
    -loadidentity \
    -lookat $Viewpoint(x) $Viewpoint(y) $Viewpoint(z) \
      $Viewpoint(x) $Viewpoint(y) 0 0 1 0 \
    -multmatrix $ViewMatrix
  $w redraw
}

proc Rotate {w x y} {
  global xlast ylast ViewMatrix Viewpoint Perspective
  set scl [expr $Viewpoint(z) / $Perspective(fd)]
  if {[expr $scl > 1.0]} {
    set dx [expr $x - $xlast]
    set dy [expr $y - $ylast]
  } else {
    if {[expr $scl < 0.1]} {set scl 0.1}
    set dx [expr $scl * ($x - $xlast)]
    set dy [expr $scl * ($y - $ylast)]
  }
  set xlast $x
  set ylast $y
  eval $w eval \
    -loadidentity \
    -rotate $dy 1 0 0 \
    -rotate $dx 0 1 0 \
    -multmatrix $ViewMatrix
  set ViewMatrix [$w get modelviewmatrix]
  set_view $w
}

proc Pan {w x y} {
  global xlast ylast Viewpoint Perspective ProgData
  set scl [expr $Viewpoint(z) / $Perspective(fd)]
  if {[expr $scl > 1.0]} {
    set scl $ProgData(scale)
  } else {
    if {[expr $scl < 0.1]} {set scl 0.1}
    set scl [expr $scl * $ProgData(scale)]
  }
  set Viewpoint(x) [expr $scl * ($xlast - $x) + $Viewpoint(x)]
  set Viewpoint(y) [expr $scl * ($y - $ylast) + $Viewpoint(y)]
  set xlast $x
  set ylast $y
  set_view $w
}

proc Zoom {w x y} {
  global xlast ylast Viewpoint
  set scl [expr 0.01 * $Viewpoint(z)]
  if {[expr $scl < 0.001]} {set scl 0.001}
  set Viewpoint(z) [expr $scl * ($ylast - $y) + $Viewpoint(z)]
  if {[expr $Viewpoint(z) < 0.001]} {
    set Viewpoint(z) 0.001
  }
  set xlast $x
  set ylast $y
  set_view $w
}

proc WheelZoom {d} {
  global Viewpoint OGLwin
  set scl [expr 0.001 * $Viewpoint(z)]
  if {[expr $scl < 0.0001]} {set scl 0.0001}
  set Viewpoint(z) [expr $scl * $d + $Viewpoint(z)]
  if {[expr $Viewpoint(z) < 0.001]} {
    set Viewpoint(z) 0.001
  }
  set_view $OGLwin
}

proc reset_clipping {} {
  global ProgData Perspective OGLwin
  array set ProgData {
    np 0.025
    fp 2.0
  }
  set Perspective(np) [expr $Perspective(fd) * $ProgData(np)]
  set Perspective(fp) [expr $Perspective(fd) * $ProgData(fp)]
  compute_viewport $OGLwin
}

proc Center {} {
  global ProgData Viewpoint ViewMatrix Perspective OGLwin
  set bounds [CGNSbounds 0 $ViewMatrix]
  set b [lindex $bounds 0]
  set xt [expr -0.5 * ([lindex $b 1] + [lindex $b 0])]
  set xd [expr double([lindex $b 1]) - double([lindex $b 0])]
  set b [lindex $bounds 1]
  set yt [expr -0.5 * ([lindex $b 1] + [lindex $b 0])]
  set yd [expr double([lindex $b 1]) - double([lindex $b 0])]
  set b [lindex $bounds 2]
  set zt [expr -0.5 * ([lindex $b 1] + [lindex $b 0])]
  set zd [expr double([lindex $b 1]) - double([lindex $b 0])]
  set dd [expr $xd * $xd + $yd * $yd + $zd * $zd]
  if {[expr $dd < 1.0e-15]} {
    set scl 1.0
  } else {
    set scl [expr 1.0 / sqrt($dd)]
  }

#  reset_clipping
  eval $OGLwin eval \
    -matrixmode modelview \
    -loadidentity \
    -scale $scl $scl $scl \
    -translate $xt $yt $zt \
    -multmatrix $ViewMatrix
  set ViewMatrix [$OGLwin get modelviewmatrix]

  set Viewpoint(x) 0
  set Viewpoint(y) 0
  set Viewpoint(z) $Perspective(fd)
  set_view $OGLwin
}

proc reset_view {} {
  global ProgData Viewpoint ViewMatrix Perspective OGLwin
  set bounds [CGNSbounds 0]
  set b [lindex $bounds 0]
  set xt [expr -0.5 * ([lindex $b 1] + [lindex $b 0])]
  set xd [expr double([lindex $b 1]) - double([lindex $b 0])]
  set b [lindex $bounds 1]
  set yt [expr -0.5 * ([lindex $b 1] + [lindex $b 0])]
  set yd [expr double([lindex $b 1]) - double([lindex $b 0])]
  set b [lindex $bounds 2]
  set zt [expr -0.5 * ([lindex $b 1] + [lindex $b 0])]
  set zd [expr double([lindex $b 1]) - double([lindex $b 0])]
  set dd [expr $xd * $xd + $yd * $yd + $zd * $zd]
  if {[expr $dd < 1.0e-15]} {
    set scl 1.0
  } else {
    set scl [expr 1.0 / sqrt($dd)]
  }

  set Viewpoint(x) 0
  set Viewpoint(y) 0
  set Viewpoint(z) $Perspective(fd)

  reset_clipping
  $OGLwin eval \
    -matrixmode modelview \
    -loadidentity \
    -scale $scl $scl $scl \
    -translate $xt $yt $zt
  set ViewMatrix [$OGLwin get modelviewmatrix]
  eval $OGLwin eval \
    -loadidentity \
    -lookat $Viewpoint(x) $Viewpoint(y) $Viewpoint(z) \
      $Viewpoint(x) $Viewpoint(y) 0 0 1 0 \
    -multmatrix $ViewMatrix
}

proc init_data {} {
  global ProgData
  set ProgData(/,dim) 0
  set ProgData(/,mode) 0
  set ProgData(/,vis) 1
  set ProgData(/,clr) ""
  set ProgData(/,msg) [CGNSsummary]
  set nz 0
  foreach z $ProgData(zones) {
    set clr [OGLcolor $nz]
    set zone /$z
    set ProgData($zone,dim) 0
    set ProgData($zone,vis) 1
    set ProgData($zone,mode) 0
    set ProgData($zone,clr) $clr
    set ProgData($zone,msg) [CGNSsummary $nz]
    set subdirs {}
    set nr 0
    foreach r $ProgData($nz,regions) {
      set reg $zone/$r
      if [catch {CGNSregiondim $nz $nr} dim] {
        set msg $dim
        set dim 0
      } else {
        set msg ""
      }
      if {$dim < 1} {
        set ProgData($reg,dim) 0
        set ProgData($reg,vis) 0
        set ProgData($reg,mode) 0
        set ProgData($reg,clr) ""
        if {$msg == ""} {
          set ProgData($reg,msg) "couldn't process for exterior faces"
        } else {
          set ProgData($reg,msg) $msg
        }
      } else {
        set ProgData($reg,dim) $dim
        if {$nr == 0 || $dim == 3} {
          set ProgData($reg,vis) $ProgData(meshvis)
          set ProgData($reg,mode) $ProgData(meshmode)
        } else {
          set ProgData($reg,vis) $ProgData(regvis)
          set ProgData($reg,mode) $ProgData(regmode)
        }
if {$ProgData(autocolor)} {
set clr [OGLcolor [expr $nz + $nr]]
}
        set ProgData($reg,clr) $clr
        set ProgData($reg,msg) [CGNSsummary $nz $nr]
      }
      set dir [file dirname $r]
      if {$dir != "" && $dir != "."} {
        if {[lsearch $subdirs $dir] < 0} {
          lappend subdirs $dir
          set ProgData($zone/$dir,dim) $ProgData($reg,dim)
          set ProgData($zone/$dir,vis) 1
          set ProgData($zone/$dir,mode) $ProgData($reg,mode)
          set ProgData($zone/$dir,clr) $clr
          set ProgData($zone/$dir,msg) ""
        } else {
          if {$ProgData($zone/$dir,dim) < $ProgData($reg,dim)} {
            set ProgData($zone/$dir,dim) $ProgData($reg,dim)
          }
          if {$ProgData($zone/$dir,mode) != $ProgData($reg,mode)} {
            set ProgData($zone/$dir,mode) 0
          }
        }
      }
      if {$nr == 0} {
          set ProgData($zone,dim) $ProgData($reg,dim)
          set ProgData($zone,mode) $ProgData($reg,mode)
      } else {
        if {$ProgData($zone,dim) < $ProgData($reg,dim)} {
          set ProgData($zone,dim) $ProgData($reg,dim)
        }
        if {$ProgData($zone,mode) != $ProgData($reg,mode)} {
          set ProgData($zone,mode) 0
        }
      }
      incr nr
    }
    if {$nz == 0} {
      set ProgData(/,dim) $ProgData($zone,dim)
      set ProgData(/,mode) $ProgData($zone,mode)
    } else {
      if {$ProgData(/,dim) < $ProgData($zone,dim)} {
        set ProgData(/,dim) $ProgData($zone,dim)
      }
      if {$ProgData(/,mode) != $ProgData($zone,mode)} {
        set ProgData(/,mode) 0
      }
    }
    incr nz
  }
}

proc init_display {} {
  global ProgData OGLwin
  if {$ProgData(displaylist) != ""} {
    $OGLwin deletelist $ProgData(displaylist)
  }
  if {$ProgData(nzones) < 1} {
    set ProgData(displaylist) ""
    $OGLwin main -clear colorbuffer depthbuffer -call $ProgData(axislist)
  } else {
    set calllist {}
    set nz 0
    foreach z $ProgData(zones) {
      set nr 0
      foreach r $ProgData($nz,regions) {
        if {$ProgData(/$z/$r,dim)} {
          set clr $ProgData(/$z/$r,clr)
          if {$ProgData(/$z/$r,vis)} {
            set dl [OGLregion $nz $nr $ProgData(/$z/$r,mode) $clr]
          } else {
            set dl [OGLregion $nz $nr 0 $clr]
          }
          lappend calllist -call $dl
        }
        incr nr
      }
      incr nz
    }
    if {$ProgData(cutlist) != ""} {
      lappend calllist -call $ProgData(cutlist)
    }
    if {$ProgData(planelist) != ""} {
      lappend calllist -call $ProgData(planelist)
    }
    set ProgData(displaylist) [eval $OGLwin newlist $calllist]
    $OGLwin main -clear colorbuffer depthbuffer \
      -call $ProgData(axislist) -call $ProgData(displaylist)
  }
  OGLaxis $ProgData(axis)
}

proc cleanup_display {} {
  global ProgData OGLwin
  display_message ""
  if {$ProgData(cgnsfile) == ""} return
  catch CGNSclose
  close_cutplane
  if {$ProgData(displaylist) != ""} {
    $OGLwin deletelist $ProgData(displaylist)
    set ProgData(displaylist) ""
  }
  $OGLwin main -clear colorbuffer depthbuffer -call $ProgData(axislist)
  TreeDelete $ProgData(tree) /
  array set ProgData {
    bases {}
    nzones 0
    zones {}
    curnode ""
    curdim 0
    curclr ""
    curmode 0
  }
  OGLaxis $ProgData(axis)
  $OGLwin redraw
}

#---------- set base

proc build_tree {} {
  global ProgData
  set_node ""
  TreeDelete $ProgData(tree) /
  TreeInsert $ProgData(tree) / -text $ProgData(curbase)
  set nz 0
  foreach z $ProgData(zones) {
    set zone /$z
    if {$ProgData($zone,dim) == 0} {
      set icon box_question
      set clr ""
    } else {
      if {$ProgData($zone,vis)} {
        set icon box_checked
      } else {
        set icon box_empty
      }
      if {$ProgData(showcolors)} {
        set clr [color_value $ProgData($zone,clr)]
      } else {
        set clr ""
      }
    }
    TreeInsert $ProgData(tree) $zone -tag $nz -icon $icon -iconbg $clr
    set nr 0
    foreach r $ProgData($nz,regions) {
      set dir [file dirname $r]
      if {$dir != "" && $dir != "." &&
        ![TreeExists $ProgData(tree) $zone/$dir]} {
        if {$ProgData($zone/$dir,dim) == 0} {
          set icon box_question
          set clr ""
        } else {
          if {$ProgData($zone/$dir,vis)} {
            set icon box_checked
          } else {
            set icon box_empty
          }
          if {$ProgData(showcolors)} {
            set clr [color_value $ProgData($zone/$dir,clr)]
          } else {
            set clr ""
          }
        }
        TreeInsert $ProgData(tree) $zone/$dir -icon $icon -iconbg $clr
      }
      set reg $zone/$r
      if {$ProgData($reg,dim) == 0} {
        set icon box_question
        set clr ""
      } else {
        if {$ProgData($reg,vis)} {
          set icon box_checked
        } else {
          set icon box_empty
        }
        if {$ProgData(showcolors)} {
          set clr [color_value $ProgData($reg,clr)]
        } else {
          set clr ""
        }
      }
      TreeInsert $ProgData(tree) $reg -tag "$nz $nr" -icon $icon -iconbg $clr
      incr nr
    }
    incr nz
  }
}

proc set_base {w nb} {
  global ProgData
  display_message "reading zone information..."
  if {[catch {CGNSbase $nb} zones]} {
    display_message $zones
    errormsg $zones
    return ""
  }
  set ProgData(zones) $zones
  set ProgData(nzones) [llength $zones]
  . configure -cursor watch
  update

  for {set nz 0} {$nz < $ProgData(nzones)} {incr nz} {
    display_message "reading zone [expr $nz + 1]..."
    if {[catch {CGNSzone $nz} regs]} {
      errormsg $regs
      set regs ""
    }
    set ProgData($nz,regions) $regs
  }

  set ProgData(curbase) [lindex $ProgData(bases) $nb]
  init_data
  display_message "building display lists..."
  init_display
  reset_view
  build_tree
  if [winfo exists .cutplane] {
    OGLcutplane 0
    center_cutplane
  }

  . configure -cursor {}
  display_message ""
  return ""
}

proc show_info {w x y node} {
  global ProgData Font
  if {$node == ""} return
  if {$node == "/"} {
    catch CGNSgetbase msg
  } else {
    set tag [TreeGet $w $node -tag]
    if {[llength $tag] == 1} {
      catch {CGNSgetzone $tag} msg
    } elseif {[llength $tag] == 2} {
      catch {CGNSgetregion [lindex $tag 0] [lindex $tag 1]} msg
    } else {
      set msg ""
    }
  }
  if {$msg != ""} {
    incr x [winfo rootx $w]
    incr y [winfo rooty $w]
    popup_message $msg -parent $w -position "$x $y" \
      -font $Font(fixed) -wrap 0
  }
}

proc set_node {node} {
  global ProgData
  TreeSelectionSet $ProgData(tree) $node
  set ProgData(curnode) $node
  if {$node == "" || ![info exists ProgData($node,dim)]} {
    set ProgData(curdim) 0
    set ProgData(curmode) 0
    set ProgData(curclr) ""
    set msg ""
  } else {
    set ProgData(curdim) $ProgData($node,dim)
    set ProgData(curmode) $ProgData($node,mode)
    set ProgData(curclr) $ProgData($node,clr)
    set msg $ProgData($node,msg)
  }
  color_button $ProgData(clrbut) $ProgData(curclr)
  options_state $ProgData(curdim)
  display_message $msg
}

proc set_visibility {node vis all} {
  global ProgData
  if {$ProgData($node,dim) == 0} return
  if {$vis} {
    set ProgData($node,vis) 1
    TreeSet $ProgData(tree) $node -icon box_checked
  } else {
    set ProgData($node,vis) 0
    TreeSet $ProgData(tree) $node -icon box_empty
  }
  if {$all} {
    foreach c [TreeGet $ProgData(tree) $node -children] {
      set_visibility $node/$c $vis 1
    }
  }
}

proc toggle_visibility {node {all 0}} {
  global ProgData
  if {$node == "" || ![info exists ProgData($node,dim)] ||
      $ProgData($node,dim) == 0} return
  if {$ProgData($node,vis)} {
    set vis 0
  } else {
    set vis 1
  }
  set_visibility $node $vis $all
  update_node $node
  if {$vis && $node != [TreeSelectionGet $ProgData(tree)]} {
    set_node $node
  }
}

proc node_visible {node} {
  global ProgData
  if {$ProgData($node,vis) == 0} {return 0}
  set node [file dirname $node]
  while {$node != "" && $node != "." && $node != "/"} {
    if {$ProgData($node,vis) == 0} {return 0}
    set node [file dirname $node]
  }
  return 1
}

proc update_children {node mode clr} {
  global ProgData
  if {$mode != ""} {
    set ProgData($node,mode) $mode
  }
  if {$clr != ""} {
    if {[llength $clr] == 1} {
      set ProgData($node,clr) [OGLcolor $clr]
    } else {
      set ProgData($node,clr) $clr
    }
    if {$ProgData(showcolors)} {
      TreeSet $ProgData(tree) $node -iconbg [color_value $ProgData($node,clr)]
    }
  }
  set children [TreeGet $ProgData(tree) $node -children]
  if {[llength $children] > 0} {
    if {$node == "/"} {set node ""}
    foreach c $children {
      update_children $node/$c $mode $clr
      if {[llength $clr] == 1} {incr clr}
    }
  } else {
    set tag [TreeGet $ProgData(tree) $node -tag]
    if {[llength $tag] == 2} {
      if [node_visible $node] {
        set mode $ProgData($node,mode)
      } else {
        set mode 0
      }
      OGLregion [lindex $tag 0] [lindex $tag 1] $mode $ProgData($node,clr)
    }
  }
}

proc update_node {node {mode ""} {clr ""}} {
  global ProgData OGLwin
  if {$node == "" || $ProgData($node,dim) == 0} return
  update_children $node $mode $clr
  OGLaxis $ProgData(axis)
  $OGLwin redraw
}

proc parent_mode {node} {
  global ProgData
  if {$node == "" || $node == "."} return
  set mode ""
  set par $node
  if {$par == "/"} {set par ""}
  foreach c [TreeGet $ProgData(tree) $node -children] {
    if {$mode == ""} {
      set mode $ProgData($par/$c,mode)
    } else {
      if {$mode != $ProgData($par/$c,mode)} {
        set mode 0
        break
      }
    }
  }
  if {$mode == ""} return
  set ProgData($node,mode) $mode
  if {$node != "/"} {
    parent_mode [file dirname $node]
  }
}

proc update_mode {} {
  global ProgData
  set node $ProgData(curnode)
  if {$node == "" || $ProgData($node,dim) == 0} return
  update_node $node $ProgData(curmode)
  if {$node != "/"} {
    parent_mode [file dirname $node]
  }
}

proc update_color {w} {
  global ProgData
  set node $ProgData(curnode)
  if {$node == "" || $ProgData($node,dim) == 0} return
  set newclr [select_color $ProgData(curclr) "Set Color" $w]
  if {$newclr != ""} {
    set ProgData(curclr) $newclr
    color_button $ProgData(clrbut) $newclr
    update_node $node "" $newclr
  }
}

proc auto_color {w} {
  global ProgData
  set node $ProgData(curnode)
  if {$node == "" || $ProgData($node,dim) == 0} return
  set clr [expr [llength [split $node /]] - 2]
  set newclr [OGLcolor $clr]
  set ProgData(curclr) $newclr
  color_button $ProgData(clrbut) $newclr
  update_node $node "" $clr
}

#---------- cutting plane

proc normalize_cutplane {} {
  global ProgData
  set r 0.0
  foreach i {x y z} {
    set r [expr $r + $ProgData(cut$i) * $ProgData(cut$i)]
  }
  if [expr $r == 0.0] {
    set r 1.0
  } else {
    set r [expr sqrt($r)]
  }
  foreach i {x y z} {
    set ProgData(cut$i) [expr $ProgData(cut$i) / $r]
  }
}

proc reset_cutplane {} {
  global ProgData
  set bounds [CGNSbounds 0]
  array set ProgData {
    rotx 0
    roty 0
    rotz 0
    lastx 0
    lasty 0
    lastz 0
  }
  set n 0
  set dist 0.0
  set diag 0.0
  foreach i {x y z} {
    set b [lindex $bounds $n]
    set r [expr [lindex $b 1] - [lindex $b 0]]
    set diag [expr $diag + $r * $r]
    set r [expr 0.5 * ([lindex $b 0] + [lindex $b 1])]
    set dist [expr $dist + $r * $ProgData(cut$i)]
    incr n
  }
  set ProgData(translate) [expr ($ProgData(cutd) - $dist) / sqrt($diag)]
}

proc draw_cut {{recompute 1}} {
  global ProgData OGLwin
  set ProgData(cutplane) \
    [list $ProgData(cutx) $ProgData(cuty) $ProgData(cutz) $ProgData(cutd)]
  if {$recompute} {
    OGLcutplane $ProgData(cutmode) $ProgData(cutplane)
  } else {
    OGLcutplane $ProgData(cutmode)
  }
  if {$ProgData(cutmode)} OGLdrawplane
  $OGLwin redraw
  reset_cutplane
}

proc draw_cutplane {} {
  global ProgData OGLwin
  if {$ProgData(drawID) == ""} {
    set ProgData(drawID) [after idle draw_cutplane]
    return
  }
  set ProgData(drawID) ""
  set plane {}
  foreach i {x y z d} {
    if [catch {expr double($ProgData(cut$i))} val] return
    lappend plane $val
  }
  if [catch {OGLdrawplane $plane} msg] {
    OGLdrawplane
  } else {
    set msg ""
  }
  $OGLwin redraw
  display_message $msg
}

proc update_cutplane {data var op} {
  global ProgData
  if {$ProgData(dotrace)} draw_cutplane
}

proc center_cutplane {} {
  global ProgData OGLwin
  set ProgData(dotrace) 0
  normalize_cutplane
  set bounds [CGNSbounds 0]
  set d 0.0
  set n 0
  foreach i {x y z} {
    set b [lindex $bounds $n]
    set r [expr 0.5 * ([lindex $b 0] + [lindex $b 1])]
    set d [expr $d + $r * $ProgData(cut$i)]
    incr n
  }
  set ProgData(cutd) [format "%g" $d]
  set ProgData(translate) 0
  draw_cutplane
  set ProgData(dotrace) 1
}

proc set_cutplane {dir} {
  global ProgData OGLwin
  set ProgData(dotrace) 0
  set bounds [CGNSbounds 0]
  if {$dir == "x"} {
    array set ProgData {
      cutx 1
      cuty 0
      cutz 0
    }
    set b [lindex $bounds 0]
  } elseif {$dir == "y"} {
    array set ProgData {
      cutx 0
      cuty 1
      cutz 0
    }
    set b [lindex $bounds 1]
  } else {
    array set ProgData {
      cutx 0
      cuty 0
      cutz 1
    }
    set b [lindex $bounds 2]
  }
  set ProgData(cutd) [format "%g" \
    [expr 0.5 * ([lindex $b 0] + [lindex $b 1])]]
  draw_cutplane
  reset_cutplane
  set ProgData(dotrace) 1
}

proc rotate_cutplane {axis amt} {
  global ProgData
  set ProgData(dotrace) 0
  normalize_cutplane
  set bounds [CGNSbounds 0]
  set n 0
  set d $ProgData(cutd)
  foreach i {x y z} {
    set b [lindex $bounds $n]
    set p$i [expr 0.5 * ([lindex $b 0] + [lindex $b 1])]
    set d [expr $d - [set p$i] * $ProgData(cut$i)]
    incr n
  }
  foreach i {x y z} {
    set r$i [expr [set p$i] + $d * $ProgData(cut$i)]
  }
  set ca [expr cos(0.0174533 * ($amt - $ProgData(last$axis)))]
  set sa [expr sin(0.0174533 * ($amt - $ProgData(last$axis)))]
  set ProgData(last$axis) $amt
  set x $ProgData(cutx)
  set y $ProgData(cuty)
  set z $ProgData(cutz)
  if {$axis == "x"} {
    set ProgData(cuty) [format "%g" [expr $ca * $y - $sa * $z]]
    set ProgData(cutz) [format "%g" [expr $ca * $z + $sa * $y]]
  } elseif {$axis == "y"} {
    set ProgData(cutx) [format "%g" [expr $ca * $x + $sa * $z]]
    set ProgData(cutz) [format "%g" [expr $ca * $z - $sa * $x]]
  } else {
    set ProgData(cutx) [format "%g" [expr $ca * $x - $sa * $y]]
    set ProgData(cuty) [format "%g" [expr $ca * $y + $sa * $x]]
  }
  set d 0.0
  foreach i {x y z} {
    set d [expr $d + $ProgData(cut$i) * [set r$i]]
  }
  set ProgData(cutd) [format "%g" $d]
  draw_cutplane
  set ProgData(dotrace) 1
}

proc translate_cutplane {amt} {
  global ProgData
  set ProgData(dotrace) 0
  normalize_cutplane
  set bounds [CGNSbounds 0]
  set n 0
  set dist 0.0
  set diag 0.0
  foreach i {x y z} {
    set b [lindex $bounds $n]
    set r [expr [lindex $b 1] - [lindex $b 0]]
    set diag [expr $diag + $r * $r]
    set r [expr 0.5 * ([lindex $b 0] + [lindex $b 1])]
    set dist [expr $dist + $r * $ProgData(cut$i)]
    incr n
  }
  set ProgData(cutd) [format "%g" [expr $dist + $amt * sqrt($diag)]]
  draw_cutplane
  set ProgData(dotrace) 1
}

proc cutplane_color {but {trans ""}} {
  global ProgData
  if {$trans == ""} {
    set clr [select_color [lrange $ProgData(cutcolor) 0 2] "Select Color" $but]
    if {$clr == ""} return
    set trans [lindex $ProgData(cutcolor) 3]
    set ProgData(cutcolor) $clr
    lappend ProgData(cutcolor) $trans
    color_button $but $clr
  } else {
    set ProgData(cutcolor) [lreplace $ProgData(cutcolor) 3 3 $trans]
  }
  OGLcutconfig $ProgData(cutcolor)
  draw_cutplane
}

proc cutplane_config {recompute} {
  global ProgData
  OGLcutconfig $ProgData(cutcolor) $ProgData(usecutclr) $ProgData(ignorevis)
  draw_cut $recompute
}

proc cutplane_defaults {but} {
  global ProgData
  set ProgData(cutcolor) {0.7 0.7 0.4 0.5}
  set ProgData(transparency) 0.5
  color_button $but $ProgData(cutcolor)
  set ProgData(cutmode) 1
  set ProgData(usecutclr) 0
  set ProgData(ignorevis) 0
  OGLcutconfig $ProgData(cutcolor) 0 0
  set_cutplane x
}

proc cutplane_control {} {
  global ProgData Font
  set w .cutplane
  if [winfo exists $w] {
    wm deiconify $w
    raise $w
    focus $w
    return
  }
  if {$ProgData(cutplane) == ""} {
    array set ProgData {
      cutx 1
      cuty 0
      cutz 0
    }
  } else {
    set ProgData(cutx) [lindex $ProgData(cutplane) 0]
    set ProgData(cuty) [lindex $ProgData(cutplane) 1]
    set ProgData(cutz) [lindex $ProgData(cutplane) 2]
  }

  set ProgData(transparency) [lindex $ProgData(cutcolor) 3]

  toplevel $w
  wm title $w "Cutting Plane"
  wm protocol $w WM_DELETE_WINDOW close_cutplane

  FrameCreate $w.clr -text "Color" -font $Font(bold) -padx 0 -pady 0
  pack $w.clr -side top -padx 5 -pady 2 -fill x -expand 1
  set clr [FrameGet $w.clr]

  label $clr.lab -text Transparency
  pack $clr.lab -side left
  scale $clr.scl -length 50 -from 0 -to 1 -resolution 0.01 \
    -orient horizontal -label "" -showvalue 0 -width 12 \
    -highlightthickness 0 -variable ProgData(transparency) \
    -command "cutplane_color $clr.but"
  pack $clr.scl -side left -fill x -expand 1
  button $clr.but -text Select -relief solid -bd 1 \
    -command "cutplane_color $clr.but"
  pack $clr.but -side right
  color_button $clr.but $ProgData(cutcolor)

  FrameCreate $w.pos -text "Cut Plane" -font $Font(bold) -padx 0 -pady 0
  pack $w.pos -side top -padx 5 -pady 2 -fill x -expand 1
  set pos [FrameGet $w.pos]

  foreach r {X Y Z} {
    set i [string tolower $r]
    set f [frame $pos.d$i]
    pack $f -side top -fill x
    label $f.lab -text "$r Direction" -width 10 -anchor w
    pack $f.lab -side left
    entry $f.ent -width 10 -textvariable ProgData(cut$i)
    pack $f.ent -side left -fill x -expand 1 -padx 2
    button $f.but -text "$r Cut" -width 6 -command "set_cutplane $i"
    pack $f.but -side right
    set f [frame $pos.r$i]
    pack $f -side top -fill x
    label $f.lab -width 10 -text "$r Rotate" -anchor e
    pack $f.lab -side left
    scale $f.scl -length 50 -from -90 -to 90 -resolution 1 \
      -orient horizontal -label "" -showvalue 0 -width 12 \
      -highlightthickness 0 -variable ProgData(rot$i) \
      -command "rotate_cutplane $i"
    pack $f.scl -side left -fill x -expand 1
  }

  set f [frame $pos.dd]
  pack $f -side top -fill x
  label $f.lab -text "Offset" -width 10 -anchor w
  pack $f.lab -side left
  entry $f.ent -width 10 -textvariable ProgData(cutd)
  pack $f.ent -side left -fill x -expand 1 -padx 2
  button $f.but -text Center -width 6 -command center_cutplane
  pack $f.but -side right
  set f [frame $pos.rd]
  pack $f -side top -fill x
  label $f.lab -width 10 -text "Translate" -anchor e
  pack $f.lab -side left
  scale $f.scl -length 50 -from -.5 -to 0.5 -resolution 0.001 \
    -orient horizontal -label "" -showvalue 0 -width 12 \
    -highlightthickness 0 -variable ProgData(translate) \
    -command translate_cutplane
  pack $f.scl -side left -fill x -expand 1

  FrameCreate $w.mode -text "Display" -font $Font(bold) -padx 0 -pady 0
  pack $w.mode -side top -padx 5 -pady 2 -fill x -expand 1
  set mode [FrameGet $w.mode]
  set f [frame $mode.d]
  pack $f -side top -fill x
  set n 0
  foreach i {Off Planar Elements Shaded} {
    radiobutton $f.f$n -text $i -variable ProgData(cutmode) \
      -value $n -command {draw_cut 0}
    pack $f.f$n -side left -expand 1
    incr n
  }
  set f [frame $mode.b]
  pack $f -side top -fill x
  checkbutton $f.cc -text "Use Cutplane Color" \
    -variable ProgData(usecutclr) -onvalue 1 -offvalue 0 \
    -command {cutplane_config 0}
  checkbutton $f.iv -text "Ignore Region Visibility" \
    -variable ProgData(ignorevis) -onvalue 1 -offvalue 0 \
    -command {cutplane_config 1}
  pack $f.cc $f.iv -side left -expand 1

  frame $w.but
  pack $w.but -side top -fill x -padx 5 -expand 1
  button $w.but.r -text Compute -width 8 -default active \
    -command draw_cut
  button $w.but.d -text Defaults -width 8 \
    -command "cutplane_defaults $clr.but"
  button $w.but.c -text Close -width 8 -command close_cutplane
  pack $w.but.r $w.but.d $w.but.c -side left -expand 1

  bind $w <Return> draw_cut

  center_cutplane
  reset_cutplane
  foreach i {x y z d} {
    trace variable ProgData(cut$i) w update_cutplane
  }
  set ProgData(dotrace) 1
}

proc close_cutplane {} {
  global ProgData OGLwin
  if {[winfo exists .cutplane]} {
    set ProgData(dotrace) 1
    foreach i {x y z d} {
      trace vdelete ProgData(cut$i) w update_cutplane
    }
    OGLcutplane 0
    OGLdrawplane
    $OGLwin redraw
    destroy .cutplane
  }
}

#---------- set variable for plotting

proc set_variable {w n} {
  global ProgData
  return ""
}

#---------- read cgns file

set FileList {
  {{CGNS Files} {.cgns .cga .cgh .cgx}}
  {{All Files} {*}}
}

proc load_cgns {{filename ""}} {
  global ProgData FileList ModelData
  global Translate Scale

  if {$filename == ""} {
    set filename [FileOpen "Open CGNS File" \
      $ProgData(cgnsfile) . $FileList]
  }
  if {$filename == ""} return
  cleanup_display
  display_message "reading base information..."
  . configure -cursor watch
  update
  set msg ""
  if {[catch {CGNSopen $filename} bases]} {
    set msg $bases
  } else {
    set msg ""
  }
  . configure -cursor {}
  if {$msg != ""} {
    display_message $msg
    errormsg $msg
    return
  }
  set ProgData(bases) $bases
  ComboboxConfig $ProgData(baselist) -values $bases -index 0 -state normal
  if {[llength $bases] < 2} {
    ComboboxConfig $ProgData(baselist) -state disabled
  }
  set_base $ProgData(baselist) 0

  wm title . "CGNSplot : [file tail $filename]"
  set ProgData(cgnsfile) $filename
}

catch {
  config_icon . [list cgnsplot cgns] \
    [list $cmd_dir $cmd_dir/images $cmd_dir/../common]
}
wm deiconify .

if {$argc} {
  set file [lindex $argv [expr $argc - 1]]
  if {[string index $file 0] != "-" && [file exists $file]} {
    load_cgns $file
  }
}

