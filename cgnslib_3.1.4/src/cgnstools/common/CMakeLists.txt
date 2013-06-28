##########
# common #
##########

if (WIN32)
  install(FILES
	balloon.tcl
	color.tcl
	combobox.tcl
	config.tcl
	dialog.tcl
	editfile.tcl
	filesel.tcl
	findfile.tcl
	frame.tcl
	help.tcl
	menubar.tcl
	tclIndex
	tclreg.tcl
	tkdir.tcl
	tree.tcl
	units.tcl
	cgns.ico
	cgns.png
	cgns-icon.xbm
	cgns-mask.xbm
	DESTINATION share)
else (WIN32)
  install(FILES
	balloon.tcl
	color.tcl
	combobox.tcl
	config.tcl
	dialog.tcl
	editfile.tcl
	filesel.tcl
	findfile.tcl
	frame.tcl
	help.tcl
	menubar.tcl
	tclIndex
	tclreg.tcl
	tkdir.tcl
	tree.tcl
	units.tcl
	cgns.ico
	cgns.png
	cgns-icon.xbm
	cgns-mask.xbm
	DESTINATION share/cgnstools)
endif (WIN32)

