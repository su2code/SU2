#----- search for specified file

proc findfile:check {type fname} {
  global tcl_platform
  if {$type == "executable" && $tcl_platform(platform) == "windows"} {
    foreach ext {.com .exe .bat} {
      if [file exists $fname$ext] {
        return $fname$ext
      }
    }
  }
  if [file $type $fname] {
    if {$type == "executable" && [file isdirectory $fname]} {
        return {}
    }
    return $fname
  }
  return {}
}

proc findfile:path {str delim} {
  global env
  set path {}
  foreach p [split $str $delim] {
    set dir {}
    foreach d [split $p /] {
      if {[string length $d] > 1 && [string index $d 0] == {$}} {
        if [info exists env([string range $d 1 end])] {
          lappend dir $env([string range $d 1 end])
        } else {
          set dir {}
          break
        }
      } else {
        lappend dir $d
      }
    }
    if [llength $dir] {
      lappend path [join $dir /]
    }
  }
  if [llength $path] {
    return [join $path $delim]
  }
  return ""
}

proc find_file {type fname args} {
  global tcl_platform env

  if {$fname == ""} {
    return {}
  }
  set dr [file dirname $fname]
  set fn [file tail $fname]
  if {$type == ""} {
    set type exists
  }

  # check if path given

  if {$dr != "."} {
    set curdir [pwd]
    catch {cd $dr}
    set filename [pwd]/$fn
    cd $curdir
    return [findfile:check $type $filename]
  }

  # check if running under windows

  if {[info exists tcl_platform(platform)] &&
    $tcl_platform(platform) == "windows"} {
    set delim ";"
  } else {
    set delim ":"
  }

  if {$args == ""} {

    # check current directory

    set filename [findfile:check $type [pwd]/$fn]
    if {$filename != ""} {
      return $filename
    }
    if {$fname != $fn} {
      return {}
    }

    # check PATH environment variable

    if [info exists env(PATH)] {
      foreach dir [split $env(PATH) $delim] {
        set filename [findfile:check $type $dir/$fn]
        if {$filename != ""} {
          return $filename
        }
      }
    }

  } else {

    # check supplied list of paths/directories

    foreach val $args {
      foreach ent [split $val] {
        if {$ent != ""} {
          set path [findfile:path $ent $delim]
          if {$path != ""} {
            foreach dir [split $path $delim] {
              set filename [findfile:check $type $dir/$fn]
              if {$filename != ""} {
                return $filename
              }
            }
          }
        }
      }
    }
  }

  return {}
}
