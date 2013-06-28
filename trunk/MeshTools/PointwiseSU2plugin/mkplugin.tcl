########################################################################
# Pointwise - Proprietary software product of Pointwise, Inc.
#           Copyright (c) 1995-2009 Pointwise, Inc.
#           All rights reserved.
#
# Plugin SDK, create plugin project tcl-script
########################################################################

namespace eval ::sdk {

    set verbose 0
    set forceOverwrite 0
    set isOverwrite 0
    set pluginName ""
    set pluginType ""
    set notPluginType ""
    set pluginFolder ""
    set tplFolder "src/plugins/templates"
    set idJump 10

    #*******************************************
    proc vmsg { msg } {
        global sdk::verbose
        if { $verbose } {
            puts "V: $msg"
        }
    }

    #*******************************************
    proc showUsage { } {
        global argv0
        global tcl_version
        puts ""
        puts "usage:"
        puts ""
        puts "    $argv0 \[-h\] \[-v\] \[-f\] -str|-uns PluginName"
        puts ""
        puts "where,"
        puts ""
        puts "    -str the plugin supports stuctured grids."
        puts ""
        puts "    -uns the plugin supports unstuctured grids."
        puts ""
        puts "    -h displays this usage."
        puts ""
        puts "    -v displays verbose runtime information."
        puts ""
        puts "    -f forces over-write of an existing plugin."
        puts ""
        puts "    PluginName is name of plugin project to create. The project"
        puts "    folder will be: './src/plugins/PluginName/'"
        puts ""
        puts "FYI..."
        puts "    tcl version '$tcl_version'..."
        puts ""
    }

    #*******************************************
    proc handleSwitch { arg } {
        global sdk::verbose
        global sdk::forceOverwrite
        global sdk::pluginType
        global sdk::notPluginType

        set ret 1

        switch -- $arg {
            -f { # force the over-write of an existing plugin
                vmsg "found switch '$arg'"
                set forceOverwrite 1
            }
            -str { # structured plugin
                vmsg "found switch '$arg'"
                set pluginType    "structured"
                set notPluginType "unstructured"
            }
            -uns { # unstructured plugin
                vmsg "found switch '$arg'"
                set pluginType    "unstructured"
                set notPluginType "structured"
            }
            -v { # force verbose output
                set verbose 1
                vmsg "found switch '$arg'"
            }
            -h { # force usage
                showUsage
                exit 1
            }
            default {
                set ret 0
            }
        }

        return $ret
    }

    #*******************************************
    proc handlePluginName { arg } {
        global sdk::pluginName
        set ret 1
        set arg [ string trim $arg ]
        if { 0 == [ string length $arg ] } {
            puts "ERROR: illegal plugin name defined '$arg'."
            set pluginName $arg
        } elseif { 0 == [ string length $pluginName ] } {
            vmsg "found plugin name '$arg'."
            set pluginName $arg
        } else {
            puts "ERROR: multiple plugin names defined."
            set ret 0
        }

        return $ret
    }

    #*******************************************
    proc postCheck {  } {
        global argv0
        global sdk::pluginName
        global sdk::pluginType
        global sdk::forceOverwrite
        global sdk::isOverwrite
        global sdk::pluginFolder

        set ret 1

        if { 0 == [ string length $pluginName ] } {
            puts "ERROR: plugin name not defined."
            set ret 0
        }

        if { 0 == [ string length $pluginType ] } {
            puts "ERROR: plugin type not defined."
            set ret 0
        }

        set criticalDirs [ list "src/plugins/shared" \
                                "src/plugins/shared/CAEP" \
                                "src/plugins/shared/PWGM" \
                                "src/plugins/shared/PWP" \
                                "src/plugins/templates/PWP" \
                                "src/plugins/templates/CAEP" ]
        foreach dd $criticalDirs {
            if { [ file isdirectory $dd ] } {
                vmsg "verified folder: '$dd'"
            } else {
                puts "ERROR: critical folder '$dd' not found."
                set ret 0
            }
        }

        if { 0 != [ string length $pluginName ] } {
            set pluginFolder "src/plugins/$pluginName"
            if { ![ file exists $pluginFolder ] } {
                vmsg "plugin folder: $pluginFolder (NEW)"
            } elseif { $forceOverwrite } {
                vmsg "plugin folder: $pluginFolder (OVER-WRITE)"
                set isOverwrite 1
            } else {
                puts "ERROR: The plugin folder '$pluginFolder' already exists."
                puts "       To force over-write use:"
                puts ""
                puts "         $argv0 -[ string range $pluginType 0 2 ] -f $pluginName"
                puts ""
                set ret 0
            }
        }

        return $ret
    }

    #*******************************************
    proc processCmdLine { theArgs } {
        set ret 1

        foreach arg $theArgs {
            set isSwitch [regexp {^-.+} $arg]
            if { !$isSwitch && ![handlePluginName $arg] } {
                set ret 0
            } elseif { $isSwitch && ![handleSwitch $arg] } {
                puts "ERROR: bad switch '$arg'."
                set ret 0
            }
        }

        if { ![postCheck] } {
            set ret 0
        }
        return $ret
    }

    #*******************************************
    proc openFile { upfp fname {access r} } {
        upvar $upfp fp
        set fp [open $fname $access]
        if { "" != $fp } {
            vmsg "opened '$fname'"
            # treat files with ".ext" as unix files
            set unix [ list ".mk" ".other" ]
            set ext [ file extension $fname ]
            set ndx [lsearch -exact $unix $ext ]
            if { -1 != $ndx } {
                # force unix eol chars
                fconfigure $fp -translation lf
                vmsg "forcing unix EOL for '$fname'"
            }
        } else {
            puts "ERROR: could not open '$fname'"
        }
        return [ expr { "" != $fp } ]
    }

    #*******************************************
    proc filterString { haystack filters } {
        foreach {needle newtxt} $filters {
            vmsg "  filter: '$needle' --> '$newtxt'"
            regsub -all $needle $haystack $newtxt haystack
        }
        return $haystack
    }

    #*******************************************
    proc filterFile { fileName { filters "" } } {
        set ret 0
        vmsg "filter '$fileName'..."
        if { [ openFile fp $fileName r+ ] } {
            # read whole file
            set contents [read $fp]
            # rewind file
            seek $fp 0
            # replace file contents with filtered data
            puts -nonewline $fp [ filterString $contents $filters ]
            set ret 1
            close $fp
        }
        return $ret
    }

    #*******************************************
    proc copyFile { srcFileName destFileName { filters "" } } {
        set ret 0
        vmsg "copy '$srcFileName' to '$destFileName'..."
        if { [ openFile sf $srcFileName r ] } {
            if { [ openFile df $destFileName w ] } {
                puts -nonewline $df [ filterString [read $sf] $filters ]
                set ret 1
                close $df
            }
            close $sf
        }
        return $ret
    }

    #*******************************************
    proc copyTemplateFile { srcFileName { destFileName "" } { filters "" } } {
        global sdk::pluginFolder
        global sdk::tplFolder

        if { [ string equal $destFileName "" ] } {
            set destFileName [ file tail $srcFileName ]
            vmsg "default dest template filename '$destFileName'"
        }
        set destFileName [ file join $pluginFolder $destFileName ]
        return [ copyFile "$tplFolder/$srcFileName" $destFileName $filters ]
    }

    #*******************************************
    proc mkDir { dir } {
        set ret 0
        if { [ file exists $dir ] } {
            set ret 1
        } else {
            file mkdir $dir
            set ret [ file exists $dir ]
        }
        return $ret
    }

    #*******************************************
    proc random {{range 100}} {
        return [ expr { int( rand() * $range ) } ]
    }

    #*******************************************
    proc nextId { matches } {
        global sdk::idJump
        set ret $idJump
        set len [llength $matches]
        if { 0 != $len } {
            incr len -1
            set idList ""
            while { $len > 0 } {
                lappend idList [lindex $matches $len]
                incr len -2
            }
            set idList [ lsort -integer -decreasing $idList ]
            set ret [expr [lindex $idList 0] + $idJump ]
        }
        return $ret
    }

    #*******************************************
    proc getPluginId { } {
        global sdk::pluginName
        global sdk::idJump
        set ret $idJump
        set idFile "src/plugins/pluginRegistry.h"

        if { [ file exists $idFile ] } {
            vmsg "looking for plugin ID in '$idFile'..."
            if { [ openFile fp $idFile r ] } {
                set contents [read $fp]
                close $fp
                set pattern "#define ID_$pluginName (\\d+)"
                # 'ret' is not changed if match fails
                if { 0 == [regexp -- "$pattern" "$contents" -> ret] } {
                    set pattern "#define ID_\\w+ (\\d+)"
                    set matches [regexp -all -inline -- "$pattern" "$contents"]
                    set ret [ nextId $matches ]
                }
            } else {
                set ret 0
            }
        }
        return $ret
    }

    #*******************************************
    proc createPlugin {  } {
        global sdk::isOverwrite
        global sdk::pluginName
        global sdk::pluginFolder
        global sdk::pluginType
        global sdk::notPluginType

        # assume all OK
        set ret 0
        set pluginId [ getPluginId ]

        if { 0 == $pluginId } {
            puts "ERROR: could not create plugin ID."
            return -2
        } elseif { $isOverwrite } {
            puts "over-writing existing plugin '$pluginName' with ID=$pluginId..."
        } else {
            puts "creating new $pluginType plugin '$pluginName' with ID=$pluginId..."
        }
        # Create folder under src/plugins/$pluginName
        if { ![ mkDir $pluginFolder ] } {
            puts "ERROR: could not create folder '$pluginFolder'."
            set ret -3
        } else {
            vmsg "created folder '$pluginFolder'..."
            # Copy src/plugins/templates/CAEP/*
            # to   src/plugins/$pluginName/*
            copyTemplateFile "CAEP/module.mk"
            copyTemplateFile "CAEP/rtCaepInstanceData.h"
            copyTemplateFile "CAEP/runtimeWrite.c"

            # copy src/plugins/templates/CAEP/CaeTemplate.vcproj
            # to   src/plugins/$pluginName/$pluginName.vcproj
            # with filter
            set filters [ list  "CaeTemplate" "$pluginName" \
                                "plugintype" "$pluginType" ]
            copyTemplateFile "CAEP/CaeTemplate.vcproj" "$pluginName.vcproj" $filters

            # copy src/plugins/templates/CAEP/module.mk
            # to   src/plugins/$pluginName/module.mk
            # with filter
            set filters [   list \
                            "CaeXxxxx"       "$pluginName" \
                            "__ID__"         "$pluginId" \
                            "CAEXXXXX"       [string toupper $pluginName] \
                            "NOTPLUGINTYPE"  [string toupper $notPluginType] ]
            copyTemplateFile "CAEP/module.mk" "module.mk" $filters
            copyTemplateFile "CAEP/modulelocal-sample.mk" "modulelocal-sample.mk" $filters
            copyTemplateFile "CAEP/rtCaepInitItems.h" "" $filters
            copyTemplateFile "CAEP/rtCaepSupportData.h" "" $filters

            # Copy src/plugins/templates/PWP/*
            # to   src/plugins/$pluginName/*
            copyTemplateFile "PWP/rtPwpInitItems.h"
            copyTemplateFile "PWP/rtPwpPluginInfo.h"
            copyTemplateFile "PWP/rtPwpVersions.h"

            # update src/plugins/modulelocal.mk
            # with filter (inserts the plugin target)
            #  1st filter removes any existing target (needed for -f switch)
            #  2nd filter (re)adds new target
            set filters [ list \
                    "\t$pluginName \\\\\n" "" \
                    "\\\$\\(NULL\\)"       "$pluginName \\\n\t\$(NULL)" ]
            filterFile "src/plugins/modulelocal.mk" $filters

            # update src/plugins/pluginRegistry.h
            # with filter (inserts the plugin target)
            #  1st filter removes any existing target (needed for -f switch)
            #  2nd filter (re)adds new target
            set filters [ list \
                "#define ID_$pluginName \[0-9\]+\n"  "" \
                "// add new plugin ID above this line"  "#define ID_$pluginName $pluginId\n// add new plugin ID above this line" ]
            filterFile "src/plugins/pluginRegistry.h" $filters
        }
        return $ret
    }

    #*******************************************
    #***                  MAIN               ***
    #*******************************************
    set tcl_ret -1
    if { [ processCmdLine $argv ] } {
        vmsg "using tcl version '$tcl_version'..."
        set tcl_ret [ createPlugin ]
    } else {
        puts "aborted!"
    }
}
if { 0 != $sdk::tcl_ret } {
    sdk::showUsage
}
exit $sdk::tcl_ret;
