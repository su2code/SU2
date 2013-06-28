array set _tclReg {
  init 0
  hasreg 0
  usereg 0
  loaded 0
  changed 0
  fname ""
  base ""
}

proc tclreg_init {args} {
  global tcl_platform tcl_version env
  global _tclReg _tclRegData

  catch tclreg_close
  array set _tclReg {
    init 1
    hasreg 0
    usereg 0
    loaded 0
    changed 0
  }
  catch {unset _tclRegData}

  if {$tcl_platform(platform) == "windows"} {
    if {[info commands registry] == ""} {
      set vers [join [split $tcl_version .] {}]
      catch {load tclreg$vers registry}
    }
    if {[info commands registry] != ""} {
      set _tclReg(hasreg) 1
      set _tclReg(usereg) 1
    }
    if {[info exists env(USERPROFILE)]} {
      set home [join [split $env(USERPROFILE) \\] /]
    } elseif {[info exists env(HOME)]} {
      set home [join [split $env(HOME) \\] /]
    } else {
      set home ~
    }
    set _tclReg(fname) "$home/_tclreg"
  } else {
    if {[info exists env(HOME)]} {
      set home $env(HOME)
    } else {
      set home ~
    }
    set _tclReg(fname) "$home/.tclreg"
  }

  foreach {tag value} $args {
    switch -- $tag {
      -usereg {
        set _tclReg(usereg) $value
      }
      -fname {
        if {$tcl_platform(platform) == "windows"} {
          set value [join [split $value \\] /]
        }
        if {[string first / $value] < 0} {
          set _tclReg(fname) "$home/$value"
        } else {
          set _tclReg(fname) $value
        }
      }
      -base {
        set _tclReg(base) $value
      }
    }
  }

  if {$_tclReg(usereg) && !$_tclReg(hasreg)} {return 0}
  return 1
}

proc tclreg_open {} {
  global _tclReg _tclRegData
  if {!$_tclReg(init)} tclreg_init
  if {$_tclReg(usereg)} return
  catch tclreg_close
  catch {unset _tclRegData}
  catch {source $_tclReg(fname)}
  set _tclReg(changed) 0
  set _tclReg(loaded) 1
}

proc tclreg_close {} {
  global _tclReg _tclRegData
  if {$_tclReg(usereg) || !$_tclReg(changed)} return
  catch {unset _tclRegTemp}
  foreach key [array names _tclRegData] {
    set _tclRegTemp($key) $_tclRegData($key)
  }
  catch {source $_tclReg(fname)}
  foreach key [array names _tclRegTemp] {
    set _tclRegData($key) $_tclRegTemp($key)
  }
  set f [open $_tclReg(fname) w+]
  puts $f "array set _tclRegData \{"
  foreach key [array names _tclRegData] {
    puts $f " \{$key\} \{$_tclRegData($key)\}"
  }
  puts $f "\}"
  close $f
  set _tclReg(changed) 0
}

proc tclreg_get {key value} {
  global _tclReg _tclRegData
  if {!$_tclReg(init)} tclreg_init
  if {$_tclReg(usereg)} {
    if {$_tclReg(hasreg)} {
      if {$_tclReg(base) != ""} {
        set key "$_tclReg(base)/$key"
      }
      set key [join [split $key /] \\]
      if {[catch {registry get $key $value} data]} {
        set data ""
      }
    } else {
      set data ""
    }
    return $data
  }
  if {!$_tclReg(loaded)} tclreg_open
  set key [join [split $key \\] /]
  if {[info exists _tclRegData($key,$value)]} {
    return $_tclRegData($key,$value)
  }
  return ""
}

proc tclreg_set {key value data} {
  global _tclReg _tclRegData
  if {!$_tclReg(init)} tclreg_init
  if {$_tclReg(usereg)} {
    if {$_tclReg(hasreg)} {
      if {$_tclReg(base) != ""} {
        set key "$_tclReg(base)/$key"
      }
      set key [join [split $key /] \\]
      registry set $key $value $data
    }
  } else {
    if {!$_tclReg(loaded)} tclreg_open
    set key [join [split $key \\] /]
    set _tclRegData($key,$value) $data
    incr _tclReg(changed)
  }
}

