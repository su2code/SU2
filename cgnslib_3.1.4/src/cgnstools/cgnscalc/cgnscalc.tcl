#!/bin/sh
# the next line restarts using wish \
exec calcwish -f "$0" "$@"

proc error_exit {msg} {
  wm withdraw .
  tk_dialog .error Error $msg error 0 Exit
  exit 1
}

if {[catch {package require Tk 8.0} msg]} {
  error_exit $msg
}

#----- get startup directory and name

set cmd_name [file tail $argv0]
set cmd_dir  [file dirname $argv0]
if {![file exists $argv0] || [file isdirectory $argv0]} {
  if {$tcl_platform(platform) == "windows"} {
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

#---------- initialize

if {$tcl_platform(platform) == "windows"} {
  set vers [join [split $tcl_version .] {}]
  if {[info commands CalcInit] == {}} {
    if {[catch {load calctcl$vers} msg]} {
      error_exit $msg
    }
  }
  catch {load tclreg$vers registry}
  if {[info exists env(USERPROFILE)]} {
    set home [join [split $env(USERPROFILE) \\] /]
  } elseif {[info exists env(HOME)]} {
    set home [join [split $env(HOME) \\] /]
  } else {
    set home ~
  }
  set InitCmds "$home/_cgnsCalc.clc"
} else {
  if {[info commands CalcInit] == {}} {
    error_exit "need to run script with calcwish"
  }
  set InitCmds "~/.cgnsCalc.clc"
}

if [catch {config_defaults 1} msg] {error_exit $msg}

set butw 5

#----- initialize

set Prefix ""
set Suffix ""

set Count 0
set Index 0
set History(-1) ""

array set ProgData {
  CMDfile {}
  CGNSfile {}
  base ""
  zone ""
  solution ""
  regs ""
  vars ""
  syms ""
  cmds ""
  symlist {}
  readonly 0
  reg,file ".cgnstools"
  reg,base "HKEY_CURRENT_USER/Software/CGNS"
  reg,key  "CGNScalc"
  reg,vals {CMDfile CGNSfile}
}

if {[tclreg_init -base $ProgData(reg,base) -fname $ProgData(reg,file)]} {
  foreach i $ProgData(reg,vals) {
    if {![catch {tclreg_get $ProgData(reg,key) $i} val] && $val != ""} {
      set ProgData($i) $val
    }
  }
  catch units_read
}

#----- window title

wm title . "CGNScalc"
wm protocol . WM_DELETE_WINDOW do_quit

bind . <Up> {select_history -1}
bind . <Down> {select_history +1}
bind . <Delete> delete_symbol

#----- pull-down menu

menubar_create {File Options Commands Plot Help}

# file menu

set m [menubar_get File]
$m add command -label "Load CGNS File..." -command {load_cgns ""}
$m add command -label "Save Variables..." -command save_cgns \
  -state disabled
$m add command -label "Edit File..." -command {do_edit ""}
$m add separator
$m add command -label "Quit" -command do_quit

# options menu

set m [menubar_get Options]
$m add command -label "Delete Symbols" -command delete_symbol
$m add command -label "Reset Calculator" -command do_reset

# commands menu

set m [menubar_get Commands]
$m add command -label "Command History..." -command do_history
$m add separator
$m add command -label "Load Command File..." -command {do_load ""}
$m add cascade -label "Save Command File..." -command do_save

# plot menu

# help menu

set m [menubar_get Help]
$m add command -label "CGNScalc Help..." -command {help_show cgnscalc}
$m add command -label "CGNS Help..." -command {help_show cgns}
$m add separator
$m add command -label "Configure Help..." -command help_setup
$m add separator
$m add command -label "Unit Conversions..." -command units_convert
$m add command -label "About..." -underline 0 -command do_about

image create photo img_about -data {\
R0lGODlhIAAgALMAAAAAAP///8DAwP//////////////////////////////////////////\
/////////yH5BAEAAAIALAAAAAAgACAAAASoUAhAq70YyM1n+GAoilTHAWOaoprpqQEKUvPX\
dvKay/xssrEYLUiyhW6S3igntGVORqGS6KSxoMAiq1JzQmHg5TdMjo3LYCQTrVKXhzWuOzyF\
XzfALXNK5GGpTn1HUW55cIJdfnhRbUqOX2tsPouRkjZnliSYmZNJMFx9a3eebTNZO6MeclY7\
g3FjlY6iZpSmrUGxtJ6gpnFyRnOcmovCw0kZyMkVEhEAADs=}

proc do_about {} {
  global ProgData
  dialog .about -1 -1 "About CGNScalc" \
"CGNScalc Version 1.1

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

image create photo img_save -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAA4ALAAAAAAQABAAQwQ90MlJqwRjgM13BpeGjOSIgQ6mdYCphW1Jtugp\
z2/6sVye8rwLMKiL3Tiwm6smUp5Cmaj0A+Utq6yrZTuJAAA7}

image create photo img_edit -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAA4ALAAAAAAQABAAQwRT0DkwK6VSHsR7H8sAjMdjnigCEMxTouOpZiXw\
2Ph9S8SLpiNG72cCPUSYjWc5y9hKUNeORgQ6WEpmZ/Vg+KrI1RdlHBwo2c9oTa2ampK1XJ5x\
RAAAOw==}

set b [frame $f.b1]
pack $b -side left -padx 5

button $b.open -image img_open -command {load_cgns ""}
set_balloon $b.open "Open a CGNS File..."

button $b.save -image img_save -command save_cgns -state disabled
set_balloon $b.save "Save Selected Results to CGNS File..."

button $b.edit -image img_edit -command {do_edit ""}
set_balloon $b.edit "Edit a File..."

pack $b.open $b.save $b.edit -side left -padx 1

#---

image create photo img_delete -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAQwQoEMlJqwVYZqv3phjAXSE5IuWYgmfbiSq4ojAX\
fup9zjjda5WazWWJAAA7}

image create photo img_reset -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAAAQwEMlJqwUY24mT91rVfV/IkShAAWirSmxLvkgs\
g+eNw/rO3zSYbRYUpoqrDHLDbDIjADs=}

set b [frame $f.b2]
pack $b -side left -padx 5

button $b.delete -image img_delete -command delete_symbol
set_balloon $b.delete "Delete Selected Symbols"

button $b.reset -image img_reset -command do_reset
set_balloon $b.reset "Reset Calculator"

pack $b.delete $b.reset -side left -padx 1

#---

image create photo img_history -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAQwQeEMlJq7046827/xogjqRIleSFltXKTi6golk8\
glIEADs=}

image create photo img_read -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgMDAwICAgP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAA4ALAAAAAAQABAAAARBEMgJnL0XjMNHpUi2dRVWOtqgehgqrmwmWdoU\
Z6J9tig19xNeC/Ao/oSoBzGIJBaNyOREGXU+d0OlhNp8Qps+QAQAOw==}

image create photo img_write -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgMDAwICAgP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAA4ALAAAAAAQABAAAAQ/EEjpqr0OvE3zRJYmdmG4nQCWVqN4TezJqeb0\
rFgs4/nksz5Sb3AoDngqohH5GjidzNATmpNOjlXWUYLNBiURADs=}

set b [frame $f.b3]
pack $b -side left -padx 5

button $b.history -image img_history -command do_history
set_balloon $b.history "Command History..."

button $b.read -image img_read -command {do_load ""}
set_balloon $b.read "Read a Command File..."

button $b.write -image img_write -command do_save
set_balloon $b.write "Write Commands to File..."

pack $b.history $b.read $b.write -side left -padx 1

#---

image create photo img_plot -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAAAQyEMlJq70W4KsvQwnVVUwZTiN5otsqpZT7Wow8\
Ix+e5JVW1juMJkHc3IwogHLJhCGflggAOw==}

button $f.plot -image img_plot -state disabled
set_balloon $f.plot "Plot Data..."

pack $f.plot -side left -padx 6

#---

image create photo img_convert -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAA4ALAAAAAAQABAAQwQt0MlJq70VSO349B0ghhwYntuIXZ7YUmaWxeu3\
3eTb4Sld8yngztYDuYbC3yQCADs=}

button $f.convert -image img_convert -command "units_convert $f.convert"
set_balloon $f.convert "Unit Conversions..."
pack $f.convert -side left -padx 6

#---

image create photo img_help -data {\
R0lGODlhEAAQALMAAAAAAIAAAACAAICAAAAAgIAAgACAgICAgMDAwP8AAAD/AP//AAAA//8A\
/wD//////yH5BAEAAAgALAAAAAAQABAAQwQiEMlJq50kX5kJ1hvShd+4mSJ4qmTrXl28ehw7\
t+j75joVAQA7}

button $f.help -image img_help -command {help_show cgnscalc}
set_balloon $f.help "CGNScalc Help..."

pack $f.help -side left -padx 6

frame .toolsep -bd 1 -height 2 -relief sunken
pack .toolsep -side top -fill x

#----- base/zone/solution selection

frame .loc
pack .loc -side top -fill x -expand 1 -pady 3

foreach i {Base Zone Solution} {
  set j [string tolower $i]
  set f [frame .loc.$j]
  pack $f -side left -fill x -expand 1 -padx 5
  label $f.lab -text $i
  pack $f.lab -side left
  set ProgData($j) [ComboboxCreate $f.cb -width 10 -edit 0 -command "set_$j"]
  pack $f.cb -side left -fill x -expand 1
  bind $f.cb.ent <ButtonRelease-3> "location_info $i $f.cb.ent"
}

proc location_info {type w} {
  global Font
  if {$type == "Solution"} {
    set msg [CalcSolnInfo]
  } else {
    set msg [Calc${type}Info]
  }
  if {$msg != ""} {
    popup_message $msg -parent $w -font $Font(fixed) -wrap 0
  }
}

frame .locsep -bd 1 -height 2 -relief sunken
pack .locsep -side top -fill x

#----- top panel

set top [frame .top]
pack $top -side top -fill x -anchor n

#----- create region list

set list [frame $top.regs]
pack $list -side left -fill both -expand 1 -padx 2 -pady 2

label $list.lab -text "Regions"
pack $list.lab -side top -fill x

scrollbar $list.vs -relief sunken -command "$list.list yview"
pack $list.vs -side right -fill y
scrollbar $list.hs -relief sunken -command "$list.list xview" -orient horiz
pack $list.hs -side bottom -fill x

set ProgData(regs) [listbox $list.list -relief sunken -width 15 -height 5 \
  -exportselection 0 -yscroll "$list.vs set" -xscroll "$list.hs set"]
pack $list.list -side left -fill both -expand 1

bind $list.list <ButtonRelease-1> load_variables
bind $list.list <ButtonRelease-3> {region_properties %W %X %y}

proc load_regions {} {
  global ProgData
  set sel [$ProgData(regs) curselection]
  if {$sel != ""} {
    set sel [$ProgData(regs) get $sel]
  }
  set index ""
  $ProgData(regs) delete 0 end
  set n 0
  foreach v [CalcRegList] {
    $ProgData(regs) insert end $v
    if {$v == $sel} {
      set index $n
    }
    incr n
  }
  if {$index != ""} {
    $ProgData(regs) selection set $index $index
    $ProgData(regs) see $index
  }
}

proc region_properties {w x y} {
  global Font
  set reg [lindex [split [$w get [$w nearest $y]] "\["] 0]
  if {$reg == ""} return
  set msg [CalcRegInfo $reg]
  if {$msg != ""} {
    incr y [winfo rooty $w]
    popup_message $msg -parent $w -position "$x $y" -font $Font(fixed) -wrap 0
  }
}

#----- create variable list

set list [frame $top.vars]
pack $list -side left -fill both -expand 1 -padx 2 -pady 2

label $list.lab -text "Variables"
pack $list.lab -side top -fill x

scrollbar $list.vs -relief sunken -command "$list.list yview"
pack $list.vs -side right -fill y
scrollbar $list.hs -relief sunken -command "$list.list xview" -orient horiz
pack $list.hs -side bottom -fill x

set ProgData(vars) [listbox $list.list -relief sunken -width 15 -height 5 \
  -exportselection 0 -yscroll "$list.vs set" -xscroll "$list.hs set"]
pack $list.list -side left -fill both -expand 1

bind $list.list <Double-Button-1> add_variable
bind $list.list <ButtonRelease-3> {variable_properties %W %X %y}

proc load_variables {} {
  global ProgData
  set sel [$ProgData(vars) curselection]
  if {$sel != ""} {
    set sel [$ProgData(vars) get $sel]
  }
  set index ""
  $ProgData(vars) delete 0 end
  set n 0
  foreach v [CalcVarList] {
    $ProgData(vars) insert end $v
    if {$v == $sel} {
      set index $n
    }
    incr n
  }
  if {$index != ""} {
    $ProgData(vars) selection set $index $index
    $ProgData(vars) see $index
  }
}

proc variable_properties {w x y} {
  global Font
  set var [lindex [split [$w get [$w nearest $y]] "\["] 0]
  if {$var == ""} return
  set msg [CalcVarInfo $var]
  if {$msg != ""} {
    incr y [winfo rooty $w]
    popup_message $msg -parent $w -position "$x $y" -font $Font(fixed) -wrap 0
  }
}

#----- symbol list

set list [frame $top.syms]
pack $list -side left -fill both -expand 1 -padx 2 -pady 2

label $list.lab -text "Symbols"
pack $list.lab -side top -fill x

scrollbar $list.vs -relief sunken -command "$list.list yview"
pack $list.vs -side right -fill y
scrollbar $list.hs -relief sunken -command "$list.list xview" -orient horiz
pack $list.hs -side bottom -fill x

set ProgData(syms) [listbox $list.list -relief sunken -width 15 -height 5 \
  -exportselection 0 -yscroll "$list.vs set" -xscroll "$list.hs set"]
pack $list.list -side left -fill both -expand 1

bind $list.list <Double-Button-1> {
  set n [%W curselection]
  if {$n != ""} {add_symbol [%W get $n]}
}
bind $list.list <ButtonRelease-3> {symbol_properties %W %X %y}

proc sort_symbols {sym1 sym2} {
  if {[string first "(" $sym1] >= 0} {
    if {[string first "(" $sym2] >= 0} {
      return [string compare $sym1 $sym2]
    }
    return 1
  }
  if {[string first "\[" $sym1] >= 0} {
    if {[string first "(" $sym2] >= 0} {
      return -1
    }
    if {[string first "\[" $sym2] >= 0} {
      return [string compare $sym1 $sym2]
    }
    return 1
  }
  if {[string first "(" $sym2] >= 0 ||
      [string first "\[" $sym2] >= 0} {
    return -1
  }
  return [string compare $sym1 $sym2]
}

proc load_symbols {} {
  global ProgData
  $ProgData(syms) delete 0 end
  set ProgData(symlist) [list]
  foreach s [lsort -command sort_symbols [CalcSymList]] {
    $ProgData(syms) insert end $s
    lappend ProgData(symlist) [lindex [split $s "\["] 0]
  }
}

proc delete_symbol {} {
  global ProgData
  set syms {}
  foreach i [$ProgData(syms) curselection] {
    lappend syms [lindex [split [lindex \
      [split [$ProgData(syms) get $i] "\("] 0] "\["] 0]
  }
  if {$syms != {}} {
    if {![dialog .delsym . {} "Delete Symbols" \
      "Delete the symbols:\n$syms" question 0 Yes No]} {
      catch {CalcDelete $syms}
      load_symbols
    }
  }
}

proc symbol_properties {w x y} {
  global Font
  set sym [lindex [split [lindex \
      [split [$w get [$w nearest $y]] "\("] 0] "\["] 0]
  if {$sym == ""} return
  set msg [CalcSymInfo $sym]
  if {$msg != ""} {
    incr y [winfo rooty $w]
    popup_message $msg -parent $w -position "$x $y" -font $Font(fixed) -wrap 0
  }
}

#----- create function buttons

set func [frame $top.func]
pack $func -side right -padx 2 -pady 2

set f [frame $func.s]
pack $f -side top -fill x -padx 2 -pady 2

label $f.lab -text "Scalar Functions"
pack $f.lab -side top -fill x

foreach s {avg sum min max len mag} {
  button $f.$s -relief raised -text $s -width $butw \
    -command "add_function $s"
  pack $f.$s -side left -expand 1
}

set f [frame $func.v]
pack $f -side top -fill x -padx 2 -pady 2

label $f.lab -text "Vector Functions"
pack $f.lab -side top -fill x

set i 0
set j 0
foreach v { \
  {deg deg add_function} \
  {arc a set_prefix} \
  {exp exp add_function} \
  {rad rad add_function} \
  {hyp h set_suffix} \
  {log log add_function} \
  {abs abs add_function} \
  {sin sin add_trigfunc} \
  {pow10 pow10 add_function} \
  {round round add_function} \
  {cos cos add_trigfunc} \
  {log10 log10 add_function} \
  {ceil ceil add_function} \
  {tan tan add_trigfunc} \
  {sqr sqr add_function} \
  {floor floor add_function} \
  {fact fact add_function} \
  {sqrt sqrt add_function}} {
  if !$i {
    incr j 1
    set ff [frame $f.$j]
    pack $ff -side left -fill y
  }
  incr i 1
  button $ff.$i -relief raised -text [lindex $v 0] -width $butw \
    -command "[lindex $v 2] [lindex $v 1]"
  pack $ff.$i -side top -fill x  -expand 1
  if {$i == 3} {
    set i 0
  }
}

#----- intrinsic constants

set f [frame $func.c]
pack $f -side top -fill x -padx 2 -pady 2

label $f.lab -text "Intrinsic Constants"
pack $f.lab -side top -fill x

set n 0
foreach b { \
  {PI pi} \
  {E e} \
  {TOL tol} \
  {DIM dim} \
  {index index} \
  {rand rand()}} {
  incr n 1
  button $f.$n -relief raised -text [lindex $b 0] -width $butw \
    -command "add_name [lindex $b 1]"
  pack $f.$n -side left -expand 1
}

#----- create operator buttons

set oper [frame .oper]
pack $oper -side top -fill x -padx 2 -pady 2 -anchor n

set f [frame $oper.but]
pack $f -side left -fill both -expand 1
button $f.1 -relief raised -text "Bksp" -width 4 -command command_backspace
button $f.2 -relief raised -text "Clear" -width 4 -command command_clear
pack $f.1 $f.2 -side top -fill x  -expand 1

set operators {
  {0 5} {1 6} {2 7} {3 8} {4 9} {" " .}
  {= ,} {+ -} {* /} {^ **} {! %} {\( \)}
  {\[ \]} {< >} {<= >=} {== !=} {& |} {? :} {\" ~}
}

set n 0
foreach o $operators {
  incr n 1
  set f [frame $oper.$n]
  pack $f -side left -fill both -expand 1
  button $f.1 -relief raised -text [lindex $o 0] -width 1 \
    -command "add_name {[lindex $o 0]}"
  button $f.2 -relief raised -text [lindex $o 1] -width 1 \
    -command "add_name {[lindex $o 1]}"
  pack $f.1 $f.2 -side top -fill x  -expand 1
}

#----- command input line

image create bitmap dnimage -foreground $fgColor(button) -data "
\#define down_width 16
\#define down_height 16
static char down_bits[] = {
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xfe, 0x7f, 0xfe, 0x7f,
   0xfc, 0x3f, 0xf8, 0x1f, 0xf0, 0x0f, 0xe0, 0x07, 0xc0, 0x03, 0x80, 0x01,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
"

image create bitmap upimage -foreground $fgColor(button) -data "
\#define up_width 16
\#define up_height 16
static char up_bits[] = {
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, 0x01, 0xc0, 0x03,
   0xe0, 0x07, 0xf0, 0x0f, 0xf8, 0x1f, 0xfc, 0x3f, 0xfe, 0x7f, 0xfe, 0x7f,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
"

set com [frame .com]
pack $com -side top -fill x -padx 2 -pady 2 -anchor n

button $com.but -text "Command" -command do_command
pack $com.but -side left -anchor w

button $com.down -image dnimage -state disabled \
  -command {select_history +1}
pack $com.down -side right -fill y

button $com.up -image upimage -state disabled \
  -command {select_history -1}
pack $com.up -side right -fill y

set ProgData(cmds) [entry $com.e -relief sunken]
pack $com.e -side left -anchor w -fill x -expand 1 -padx 2

bind $com.e <Return> {do_command}
bind $com.e <Escape> {%W delete 0 end}
bind $com.e <Control-Return> add_variable
#focus $com.e

#----- output window

set msg [frame .msg -relief groove -bd 3]
pack $msg -side top -fill both -expand 1 -padx 2 -pady 2

scrollbar $msg.scroll -relief sunken -command "$msg.text yview"
pack $msg.scroll -side right -fill y

text $msg.text -height 10 -yscroll "$msg.scroll set" -cursor {}
pack $msg.text -side left -fill both -expand 1

bind $msg.text <FocusIn> "focus $com.e"

#----- help

proc help_menu {} {
  if {[help_valid cgnscalc]} {
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

help_init {cgnscalc CGNScalc} {cgns CGNS}

#----- procedures

proc command_insert {text} {
  global ProgData
  if {[$ProgData(cmds) selection present]} {
    set first [$ProgData(cmds) index sel.first]
    set last  [$ProgData(cmds) index sel.last]
    $ProgData(cmds) selection clear
    $ProgData(cmds) delete $first $last
    $ProgData(cmds) icursor $first
  }
  $ProgData(cmds) insert insert $text
  focus $ProgData(cmds)
}

proc command_clear {} {
  global ProgData
  if {[$ProgData(cmds) selection present]} {
    set first [$ProgData(cmds) index sel.first]
    set last  [$ProgData(cmds) index sel.last]
    $ProgData(cmds) selection clear
    $ProgData(cmds) delete $first $last
    $ProgData(cmds) icursor $first
  } else {
    $ProgData(cmds) delete 0 end
  }
}

proc command_backspace {} {
  global ProgData
  set n [$ProgData(cmds) index insert]
  if {$n > 0} {
    $ProgData(cmds) delete [expr $n - 1] $n
  }
}

proc writemsg {msg} {
  .msg.text insert end "$msg\n"
  .msg.text yview -pickplace end
  update idletasks
}

proc add_variable {} {
  global ProgData
  set n [$ProgData(vars) curselection]
  if {$n == ""} return
  set name [lindex [split [$ProgData(vars) get $n] "\["] 0]
  if {[lsearch -exact $ProgData(symlist) $name] >= 0} {
    add_name "\\$name"
  } else {
    add_name $name
  }
}

proc add_name {name} {
  set_prefix ""
  set_suffix ""
  command_insert "$name"
}

proc add_function {name} {
  set_prefix ""
  set_suffix ""
  command_insert "$name\("
}

proc add_symbol {sym} {
  if {[string first \[ $sym] != -1} {
    add_name [lindex [split $sym \[] 0]
  } elseif {[string first \( $sym] != -1} {
    if [string match "*(0)" $sym] {
      add_name "[lindex [split $sym \(] 0]\(\)"
    } else {
      add_function [lindex [split $sym \(] 0]
    }
  } else {
    add_name $sym
  }
}

proc set_prefix {p} {
  global Prefix
  if {$p != "" && $p == $Prefix} {set p ""}
  if {$p == ""} {
    .top.func.v.1.2 configure -relief raised
  } else {
    .top.func.v.1.2 configure -relief sunken
  }
  set Prefix $p
}

proc set_suffix {s} {
  global Suffix
  if {$s != "" && $s == $Suffix} {set s ""}
  if {$s == ""} {
    .top.func.v.2.2 configure -relief raised
  } else {
    .top.func.v.2.2 configure -relief sunken
  }
  set Suffix $s
}

proc add_trigfunc {name} {
  global Prefix Suffix
  command_insert "$Prefix$name$Suffix\("
  set_prefix ""
  set_suffix ""
}

proc do_quit {} {
  global ProgData
  CalcDone
  foreach i $ProgData(reg,vals) {
    catch {tclreg_set $ProgData(reg,key) $i $ProgData($i)}
  }
  catch units_write
  catch tclreg_close
  destroy .
  exit 0
}

#----- reset calculator

proc do_reset {} {
  global ProgData
  if [dialog .reset $ProgData(vars) {} Reset \
"This will delete all symbols and clear the output window.
Do you wish to proceed ?" question 0 Yes No] return
  CalcReset
  $ProgData(cmds) delete 0 end
  .msg.text delete 1.0 end
  load_symbols
  clear_history
}

#----- command history

proc update_history {} {
  global Index Count
  if {$Index <= 0} {
    .com.up configure -state disabled
  } else {
    .com.up configure -state normal
  }
  if {$Index >= $Count} {
    .com.down configure -state disabled
  } else {
    .com.down configure -state normal
  }
}

proc add_history {cmd} {
  global Index Count History
  set History($Count) $cmd
  incr Count
  set Index $Count
  if {[winfo exists .histlist]} {
    .histlist.top.list insert end $cmd
    .histlist.top.list selection clear 0 end
    .histlist.top.list see end
  }
  update_history
}

proc clear_history {} {
  global Index Count History
  set Index 0
  set Count 0
  catch {array unset History}
  set History(-1) ""
  if {[winfo exists .histlist]} {
    .histlist.top.list delete 0 end
  }
  update_history
}

proc select_history {dir} {
  global ProgData Index Count History
  set n [expr $Index$dir]
  if {$n >= 0 && $n <= $Count} {
    set Index $n
    $ProgData(cmds) delete 0 end
    if {$Index < $Count} {
      command_insert $History($Index)
    }
    if {[winfo exists .histlist]} {
      .histlist.top.list selection clear 0 end
      if {$Index < $Count} {
        .histlist.top.list selection set $Index
      }
    }
    update_history
  }
}

proc set_history {} {
  global ProgData Index Count History
  set n [.histlist.top.list curselection]
  if {$n != ""} {
    set Index $n
    $ProgData(cmds) delete 0 end
    command_insert $History($Index)
    update_history
  }
}

proc do_history {} {
  global Index Count History
  set w .histlist
  if {[winfo exists $w]} {
    wm deiconify $w
    focus $w
    return
  }
  toplevel $w
  wm title $w "Command History"

  frame $w.top
  pack $w.top -fill both -padx 2 -pady 2 -expand 1
  scrollbar $w.top.ys -command "$w.top.list yview"
  pack $w.top.ys -side right -fill y
  scrollbar $w.top.xs -orient horizontal -command "$w.top.list xview"
  pack $w.top.xs -side bottom -fill x

  listbox $w.top.list -width 40 -height 10 -selectmode single \
    -exportselection 0 -xscroll "$w.top.xs set" -yscroll "$w.top.ys set"
  pack $w.top.list -side left -fill both
  for {set n 0} {$n < $Count} {incr n} {
    $w.top.list insert end $History($n)
  }
  if {$Index >= 0 && $Index < $Count} {
    $w.top.list selection set $Index
    $w.top.list see $Index
  } else {
    $w.top.list see end
  }

  bind $w.top.list <Double-Button-1> set_history

  frame $w.but
  pack $w.but -side bottom -fill x
  button $w.but.select -text Select -width 8 -command set_history
  button $w.but.clear -text Clear -width 8 -command clear_history
  button $w.but.close -text Close -width 8 -command "destroy $w"
  pack $w.but.select $w.but.clear $w.but.close -side left \
    -padx 5 -pady 2 -expand 1

  center_window $w .
  focus $w
}

#----- process command

proc is_equation {cmd} {
  regsub -all "\[ \t\]+" $cmd {} newcmd
  set s [split $newcmd "="]
  if {[llength $s] > 1 && [lindex $s 1] != "" &&
    [string match "\[_a-zA-Z\]*\(\[0-9\]\)" [lindex $s 0]]} {
    return 1
  }
  return 0
}

proc do_command {} {
  global ProgData
  set cmd [string trim [$ProgData(cmds) get]]
  set n [string first "\#" $cmd]
  if {$n >= 0} {
    if {$n == 0} {
      set cmd ""
    } else {
      incr n -1
      set cmd [string trim [string range $cmd 0 $n]]
    }
  }
  if {$cmd == ""} return
  if {[string first @ $cmd] > 0} {
    set l [split $cmd @]
    set cmd [string trim [lindex $l 0]]
    set reg "@[string trim [lindex $l 1]]"
  } else {
    set regnum [$ProgData(regs) curselection]
    if {$regnum != "" && $regnum > 0 && ![is_equation $cmd]} {
      set regname [$ProgData(regs) get $regnum]
      set n [string last "\[" $regname]
      if {$n < 0} {
        set n [string last "\{" $regname]
      }
      set reg "@[string range $regname 0 [expr $n - 1]]"
    } else {
      set reg ""
    }
  }
  writemsg "\nCOMMAND: $cmd$reg"
  add_history $cmd
  if [catch {CalcCommand $cmd$reg} result] {
    foreach line [split $result '\n'] {
      writemsg "ERROR  : $line"
    }
  } else {
    writemsg "RESULT : $result"
    $ProgData(cmds) delete 0 end
    load_symbols
  }
}

#----- load/save command file

proc do_load {filename} {
  global ProgData
  if {$filename == {}} {
    set filename [FileOpen "Load Command File" $ProgData(CMDfile) . \
      {{{CMD Files} {.clc .cmd}} {{All Files} {*}}}]
  }
  if {$filename == ""} return
  update
  if {![file isfile $filename]} {
    errormsg "$filename is not a regular file"
    return
  }
  if [catch {open $filename r} f] {
    errormsg $f
    return
  }
  writemsg "LOADING: $filename"
  set region ""

  while {[gets $f cmd] >= 0} {
    while (1) {
      set n [expr [string length $cmd] - 1]
      if {[string index $cmd $n] != "\\"} break
      set cmd [string trim [string range $cmd 0 [expr $n - 1]]]
      if {[gets $f str] < 0} break
      append cmd [string trim $str]
    }
    set n [string first "\#" $cmd]
    if {$n >= 0} {
      if {$n == 0} {
        set cmd ""
      } else {
        incr n -1
        set cmd [string trim [string range $cmd 0 $n]]
      }
    } else {
      set cmd [string trim $cmd]
    }
    switch [string index $cmd 0] {
      "%" continue
      "?" continue
      "*" continue
      "@" {
        set region [string trim [string range $cmd 1 end]]
        continue
      }
      "&" {
        set cgnsfile [string trim [string range $cmd 1 end]]
        if {$cgnsfile != ""} {
          load_cgns $cgnsfile
        }
      }
      "<" {
        set cmdfile [string trim [string range $cmd 1 end]]
        if {$cmdfile != ""} {
          do_load $cmdfile
        }
      }
      ">" {set cmd [string trim [string range $cmd 1 end]]}
    }
    if {$cmd != ""} {
      add_history $cmd
      if {$region != "" && [string first "@" $cmd] < 0} {
        append cmd "@$region"
      }
      writemsg "\nCOMMAND: $cmd"
      if [catch {CalcCommand $cmd} result] {
        foreach line [split $result '\n'] {
          writemsg "ERROR  : $line"
        }
      } else {
        writemsg "RESULT : $result"
      }
    }
  }
  close $f
  load_symbols
  set ProgData(CMDfile) $filename
}

proc do_save {} {
  global ProgData
  set filename [FileSave "Save Command File" $ProgData(CMDfile) . \
    {{{Command Files} {.clc .cmd .cnv}} {{All Files} {*}}}]
  if {$filename == ""} return
  if [catch {open $filename w+} f] {
    errormsg $f
    return
  }
  writemsg "\nSAVING : $filename"
  set cmd ""
  foreach line [split [.msg.text get 1.0 end] "\n"] {
    if {$line == ""} {
      set cmd ""
      puts $f ""
    } else {
      puts $f "\#$line"
      set prompt [string range $line 0 7]
      if {$prompt == "COMMAND:"} {
        set cmd [string trim [string range $line 8 end]]
      } elseif {$prompt == "RESULT :"} {
        if {$cmd != ""} {
          puts $f $cmd
          set cmd ""
        }
      } else {
        set cmd ""
      }
    }
  }
  close $f
  set ProgData(CMDfile) $filename
}

#----- set base, zone and solution

proc set_base {w base} {
  global ProgData
  foreach f {regs vars} {
    $ProgData($f) delete 0 end
  }
  ComboboxConfig $ProgData(solution) -values {}
  if {[catch {CalcBase [expr $base + 1]} zones]} {
    errormsg $zones
    set zones ""
  } else {
    ComboboxConfig $w -index $base
  }
  ComboboxConfig $ProgData(zone) -values $zones
  if {$zones == ""} {
    load_variables
  } else {
    set_zone $ProgData(zone) 0
  }
  return ""
}

proc set_zone {w zone} {
  global ProgData
  foreach f {regs vars} {
    $ProgData($f) delete 0 end
  }
  if {[catch {CalcZone [expr $zone + 1]} solns]} {
    errormsg $solns
    set solns ""
  } else {
    ComboboxConfig $w -index $zone
  }
  ComboboxConfig $ProgData(solution) -values $solns
  load_regions
  if {$solns == ""} {
    load_variables
  } else {
    set_solution $ProgData(solution) 0
  }
  return ""
}

proc set_solution {w soln} {
  global ProgData
  $ProgData(vars) delete 0 end
  if {[catch {CalcSoln [expr $soln + 1]} vars] || $vars == ""} {
    errormsg $vars
  } else {
    ComboboxConfig $w -index $soln
  }
  load_variables
  return ""
}

#----- load/save CGNS file

proc load_cgns {fname} {
  global ProgData
  if {$fname == ""} {
    set fname [FileOpen "Load CGNS File" $ProgData(CGNSfile) . \
      {{{CGNS Files} {.cgns .cga .cgh .cgx}} {{All Files} *}}]
    if {$fname == ""} return
  }
  if {![file isfile $fname]} {
    errormsg "$fname is not a regular file"
    return
  }
  foreach f {regs vars syms} {
    $ProgData($f) delete 0 end
  }
  foreach i {base zone solution} {
    ComboboxConfig $ProgData($i) -values {}
  }
if 0 {
  if {[file writable $fname]} {
    set ProgData(readonly) 0
    set mode w
  } else {
    set ProgData(readonly) 1
    set mode r
  }
} else {
  set ProgData(readonly) 1
  set mode r
}
  dialog .open -1 -1 "Opening..." \
    "Opening and reading [file tail $fname]" hourglass 0 {}
  .open configure -cursor watch
  . configure -cursor watch
  update
  if {[catch {CalcInit $mode $fname} bases]} {
    set msg $bases
  } else {
    set msg ""
  }
  destroy .open
  . configure -cursor {}
  if {$msg != ""} {
    errormsg $bases
    return
  }
  if {$bases == ""} {
    errormsg "no bases defined in CGNS file"
    return
  }

  if {$ProgData(readonly)} {
    wm title . "CGNScalc : [file tail $fname] (read only)"
  } else {
    wm title . "CGNScalc : [file tail $fname]"
  }
  set ProgData(CGNSfile) $fname
  ComboboxConfig $ProgData(base) -values $bases
  set_base $ProgData(base) 0
}

proc save_cgns {} {
}

wm withdraw .
update idletasks
wm minsize . [winfo reqwidth .] [winfo reqheight .]
catch {
  config_icon . [list cgnscalc cgns] \
    [list $cmd_dir $cmd_dir/images $cmd_dir/../common]
}
focus $ProgData(cmds)
wm deiconify .

if {$argc} {
  set file [lindex $argv [expr $argc - 1]]
  if {[string index $file 0] != "-" && [file exists $file]} {
    load_cgns $file
  }
}

if {[file exists $InitCmds] && [file readable $InitCmds]} {
  do_load $InitCmds
}

