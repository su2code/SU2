#!/bin/sh
# the next line restarts using wish \
exec adfwish -f "$0" "$@"

proc error_exit {msg} {
  wm withdraw .
  tk_dialog .error Error $msg error 0 Exit
  exit 1
}

if {[catch {package require Tk 8.0} msg]} {
  error_exit $msg
}

#----- check for CGNS version

if {$argc} {
  set version [lindex $argv [expr $argc - 1]]
  if {[catch {expr $version < 1.1 || $version > 3.2} ok] || !$ok} {
    error_exit "invalid CGNS version number : $version"
  }
} else {
  set version ""
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

set auto_path "$cmd_dir $cmd_dir/../common $auto_path"
if {[info exists env(TCL_PROC_DIR)]} {
  set auto_path "$env(TCL_PROC_DIR) $auto_path"
}
if [catch {config_defaults 1} msg] {
  error_exit $msg
}

if {$tcl_platform(platform) == "windows"} {
  set vers [join [split $tcl_version .] {}]
  catch {load tclreg$vers registry}
}

#---------- initialize

set file ".cgnstools"
set base "HKEY_CURRENT_USER/Software/CGNS"
catch {tclreg_init -base $base -fname $file}

help_init {cgns CGNS}
cgns_init $version
cgns_tree .

catch {
  config_icon . cgns [list $cmd_dir $cmd_dir/images $cmd_dir/../common]
}

