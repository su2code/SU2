# simple editor procedures

proc editfile:save {w} {
  global _EditFile
  set fd [open $_EditFile($w,name) w+]
  puts $fd [string trimright [$w.edit.text get 1.0 end]]
  close $fd
  set _EditFile($w,dirty) 0
}

proc editfile:save_as {w} {
  global _EditFile tcl_platform
  while {1} {
    set file [FileSave "Save As..." $_EditFile($w,name) $w]
    if {$file == ""} return
    if {[file exists $file] && ![file writable $file]} {
      errormsg "can't save as $file - file is not writeable"
    } elseif {[catch {open $file w+} fd]} {
      errormsg $fd
    } else {
      puts $fd [string trimright [$w.edit.text get 1.0 end]]
      close $fd
      set _EditFile($w,dirty) 0
      if {$tcl_platform(platform) == "windows"} {
        set file [join [split $file /] \\]
      }
      set _EditFile($w,name) $file
      wm title $w $file
      return
    }
  }
}

proc editfile:reload {w} {
  global _EditFile
  if {$_EditFile($w,dirty) && [dialog .reload $w {} "Reload" \
"\"$_EditFile($w,name)\" has changed.
Do you want to discard the changes and reload the file ?" \
     warning 0 Yes No Cancel]} return
  $w.edit.text delete 1.0 end
  if {![catch {open $_EditFile($w,name) r} fd]} {
    $w.edit.text insert end [read $fd]
    close $fd
  }
  set _EditFile($w,dirty) 0
}

proc editfile:close {w} {
  global _EditFile
  if {$_EditFile($w,dirty)} {
    set result [dialog .save $w {} "Save" \
"\"$_EditFile($w,name)\" has changed.
Do you want to save the changes ?" \
     warning 0 Yes No Cancel]
    if {$result == 2} return
    if {$result == 0} {
      if {[file exists $_EditFile($w,name)] && \
         ![file writable $_EditFile($w,name)]} {
        editfile:save_as $w
      } else {
        editfile:save $w
      }
    }
  }
  catch {focus $_EditFile($w,focus)}
  destroy $w
  if {$_EditFile($w,oldGrab) != ""} {
    if {$_EditFile($w,grabStatus) == "global"} {
      grab -global $_EditFile($w,oldGrab)
    } else {
      grab $_EditFile($w,oldGrab)
    }
  }
  destroy $w
}

proc editfile:endhelp {w} {
  global _EditFile
  destroy $w.edit.help $w.edit.sb $w.done
  pack $w.but -side bottom -fill x
  pack $w.edit.ys -side right -fill y
  pack $w.edit.xs -side bottom -fill x
  pack $w.edit.text -side top -fill both -expand 1
  wm title $w $_EditFile($w,name)
  bind $w <Alt-c> {}
}

proc editfile:help {w} {
  set width  [$w.edit.text cget -width]
  set height [$w.edit.text cget -height]
  catch {pack forget $w.edit.text $w.edit.ys $w.edit.xs $w.but}
  button $w.done -text Cancel -underline 0 -pady 0 \
    -highlightthickness 0 -command "editfile:endhelp $w"
  pack $w.done -side bottom -pady 3
  scrollbar $w.edit.sb -command "$w.edit.help yview"
  text $w.edit.help -width $width -height $height -wrap word -cursor {} \
    -yscrollcommand "$w.edit.sb set"
  pack $w.edit.sb -side right -fill y
  pack $w.edit.help -side top -fill both -expand 1

  $w.edit.help insert end {This is a simple editor which supports the\
standard Motif editing characters, in addition to many of the Emacs\
editing characters. The view may be adjusted using the scrollbar or\
by pressing mouse button 2 in the window and dragging. Pressing mouse\
button 1 will set the location of the insertion cursor. Pressing and\
dragging will select a range of characters. Once the button is\
released, the selection may be adjusted by pressing button 1 with\
the shift key down. A double-click will select whole words and a\
triple-click selects lines.

To delete text, select the characters you'd like to delete and type\
Backspace or Delete.  Alternatively, you can type new text, in which\
case it will replace the selected text. To copy the selection into this\
window, select what you want to copy (either here or in another\
application), then click button 2 to copy the selection to the point of\
the mouse cursor.

Backspace and Control-h erase the character to the left of the insertion\
cursor. Delete and Control-d erase the character to the right of the\
insertion cursor. Meta-backspace deletes the word to the left of the\
insertion cursor, and Meta-d deletes the word to the right of the\
insertion cursor. Control-k deletes from the insertion cursor to the\
end of the line, or it deletes the newline character if that is the\
only thing left on the line. Control-o opens a new line by inserting a\
newline character to the right of the insertion cursor. Control-t\
transposes the two characters on either side of the insertion cursor.
}

  $w.edit.help configure -state disabled
  wm title $w "Editor Help"
  bind $w <Alt-c> "editfile:invoke $w.done"
}

proc editfile:invoke {but} {
  if {[winfo viewable $but] && [$but cget -state] == "normal"} {
    $but configure -relief sunken
    update idletasks
    after 250
    $but configure -relief raised
    $but invoke
  }
}

proc edit_file {w fname {parent ""}} {
  global Font _EditFile
  if [winfo exists $w] {
    wm deiconify $w
    raise $w
    focus $w
    return
  }
  toplevel $w
  wm title $w $fname
  wm protocol $w WM_DELETE_WINDOW "editfile:close $w"
  if {$parent != "" && [winfo exists $parent]} {
    wm transient $w [winfo toplevel $parent]
  }

  set _EditFile($w,name) $fname
  set _EditFile($w,dirty) 0
  set _EditFile($w,cursor) "1.0"
  set _EditFile($w,focus) [focus]
  set _EditFile($w,oldGrab) ""

  set f [frame $w.but]
  pack $f -side bottom -fill x
  label $f.lab -textvariable _EditFile($w,cursor) -width 8
  button $f.save -width 8 -text Save -underline 0 -pady 0 \
    -highlightthickness 1 -command "editfile:save $w"
  button $f.saveas -width 8 -text "Save As" -underline 5 -pady 0 \
    -highlightthickness 1 -command "editfile:save_as $w"
  button $f.reload -width 8 -text Reload -underline 0 -pady 0 \
    -highlightthickness 1 -command "editfile:reload $w"
  button $f.close -width 8 -text Close -underline 0 -pady 0 \
    -highlightthickness 1 -command "editfile:close $w"
  button $f.help -width 8 -text Help -underline 0 -pady 0 \
    -highlightthickness 1 -command "editfile:help $w"
  pack $f.lab $f.save $f.saveas $f.reload $f.close $f.help \
    -side left -expand 1 -pady 3

  if {[file exists $fname] && ![file writable $fname]} {
    $f.save configure -state disabled
  }

  set f [frame $w.edit]
  pack $f -side top -fill both -expand 1 -padx 2
  scrollbar $f.ys -command "$f.text yview"
  pack $f.ys -side right -fill y
  scrollbar $f.xs -orient horizontal -command "$f.text xview"
  pack $f.xs -side bottom -fill x
  text $f.text -relief sunken -font $Font(fixed) -width 60 -height 15 \
    -wrap none -xscroll "$f.xs set" -yscroll "$f.ys set"
  pack $f.text -side top -fill both -expand 1

  if {$parent != "" && [winfo exists $parent]} {
    center_window $w $parent
    set _EditFile($w,oldGrab) [grab current $w]
    if {$_EditFile($w,oldGrab) != ""} {
      set _EditFile($w,grabStatus) [grab status $_EditFile($w,oldGrab)]
    }
    catch {grab $w}
  }

  if {[file exists $fname] && [file readable $fname]} {
    set fd [open $fname r]
    $w.edit.text insert end [read $fd]
    close $fd
  }
  $w.edit.text mark set insert 1.0
  focus $w.edit.text

  bind $w.edit.text <KeyPress> {
    if {[string compare %A ""]} {set _EditFile([winfo toplevel %W],dirty) 1}
  }
  bind $w.edit.text <KeyRelease> {
    set _EditFile([winfo toplevel %W],cursor) [%W index insert]
  }
  bind $w.edit.text <ButtonRelease-1> {
    set _EditFile([winfo toplevel %W],cursor) [%W index insert]
  }
  bind $w.edit.text <Control-v> {tk_textPaste %W; break}

  bind $w.edit.text <Alt-s> "editfile:invoke $w.but.save"
  bind $w.edit.text <Alt-a> "editfile:invoke $w.but.saveas"
  bind $w.edit.text <Alt-r> "editfile:invoke $w.but.reload"
  bind $w.edit.text <Alt-c> "editfile:invoke $w.but.close"
  bind $w.edit.text <Alt-h> "editfile:invoke $w.but.help"
}

proc do_edit {{fname ""}} {
  if {$fname == ""} {
    set fname [FileOpen "Edit File..." {}]
  }
  if {$fname != ""} {
    set n 0
    while {1} {
      if ![winfo exists .edit$n] {
        edit_file .edit$n $fname
        return
      }
      incr n
    }
  }
}
