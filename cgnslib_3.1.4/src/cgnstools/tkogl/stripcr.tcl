#
# Rewrites all .c and .h files stripping CR characters
# inserted by DOS/Windows text editors
#
foreach fname [glob *.c *.h] {
    puts $fname
    set f [open $fname]
    set g [open foobar w]
    while {! [eof $f]} {
	gets $f line
	puts $g $line
    }
    close $f
    close $g
    file delete -force $fname
    file rename foobar $fname
}
