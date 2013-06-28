/****************************************************************************
 *
 * Pointwise Plugin utility functions
 *
 * Proprietary software product of Pointwise, Inc.
 * Copyright (c) 1995-2012 Pointwise, Inc.
 * All rights reserved.
 *
 ***************************************************************************/

********************
*** Introduction ***
********************

This document covers the process of creating a "starting" plugin project. 
You can use this project to implement your plugins export logic.

It is suggested you verify your Pointwise CAE plugin SDK installation before 
making any of the changes detailed below.


Verify on windows (Visual Studio 2008):
---------------------------------------
Open the VisualStudio solution file .../PluginSDK/PluginSDK.sln (created in VS2008).
Select the "Build/Batch Build..." menu.
Press the "Select All" button.
Press the "Rebuild" button.

The CaeUnsXML and CaeStrXML sample plugins should build for all platforms 
(win32, win64) and build types (debug, release). The resulting DLL plugins are 
located under the .../PluginSDK/dist/<platform> folders.

NOTE: These instructions assume your installation is VisualStudio2008 and that 
it supports building both win32 and the x64 platform types. Using other 
versions of Visual Studio is untested and unsupported.

Verify on unix and macOSX:
---------------------------------------
Open a command shell.
cd to .../PluginSDK/.
Run "make plugins-dr"

The CaeUnsXML and CaeStrXML sample plugins should build for the platform
(linux, linux_x86_64, macosx) and all build types (debug, release). The 
resulting plugins (.so or .dynlib) are located under the 
.../PluginSDK/dist/<platform> folders.


*************************
*** tclsh Information ***
*************************

If you have access to a tcl shell (tclsh) on your development machine, you can
use the included scripts (mkplugin and mkplugin.bat) to easily create a new CAE
plugin project. For details, see the "Create a New CAE Plugin Project with the
mkplugin Script" section below.

The script manages the proper integration of the new plugin project into the
SDK framework. For example;

    * Each plugin is assigned a unique id.
    * pluginRegistry.h is updated.
    * A make target is added.
    * Some plugin runtime data is initialized.

It is recommended that you use the script, but if you do not have access to a
tcl shell, see the "Manual Steps to Create a New CAE Plugin Project" section
below.

A tcl shell is included with your Pointwise installation. The pointwise tcl shell
is invoked on unix using the script:

    /your/pointwise/installation/pointwise -b

and on winOS using:

    bin/tclsh

If Pointwise is not installed on your development machine, you will need to
download and install a tcl shell on your machine. There are free installations
available on the web. Try http://www.activestate.com/activetcl/ for a good tcl
shell.

The scripts assume that the tcl shell (tclsh) is can be launched by name only.
That is, typing "tclsh" on the command line will start the tcl shell without
errors. If not already working, you may want to add the tclsh folder to the
shell's path or create an alias (unix) or a tclsh.bat file (winOS).

If you do not want to make shell changes as outlined above, you can invoke the
underlying tcl file using the following command line:

/path/to/your/tclsh mkplugin.tcl <command-line-arguments>

The <command-line-arguments> are detailed in the script section.


**************************************
*** A Few Comments Before We Begin ***
**************************************

A Pointwise CAE Plugin must be declared as exporting "Structured" or 
"Unstructured" grids. The CAE solver being targeted determines the type. If a 
solver supports both structured and unstructured grids, 2 separate plugins 
must be created.

The set of PluginSDK API calls available to a given plugin will be determined 
by the plugin type. If you get an "undefined" function compile error, you may 
be using an API call unsupported for your plugin's type.

Though not technically required, it is suggested that the plugin name start
with a prefix designating the type of plugin. The following are suggested
prefixs:

	API           Prefix  Example Plugin Name(s)
	------------  ------  ------------------------
	Export-CAE    Cae     CaeOpenFOAM, CaeXml
	
In the sections below, we will be creating the new plugin named
"CaeMyExporter". When complete, the cross-platform build will produce binaries
similar to the following examples:

for BUILD=Release
    "/path/to/your/PluginSDK/dist/win32/plugins/CaeMyExporter.dll"
    "/path/to/your/PluginSDK/dist/linux/plugins/libCaeMyExporter.so"

for BUILD=Debug
    "/path/to/your/PluginSDK/dist/win32/plugins/debug/CaeMyExporterd.dll"
    "/path/to/your/PluginSDK/dist/linux/plugins/debug/libCaeMyExporterd.so"


****************************************************************
*** Create a New CAE Plugin Project with the mkplugin Script ***
****************************************************************

Use the following steps to create a new CAE Plugin project:

Open a shell on your development platform.

Change to the "/path/to/your/PluginSDK" folder.

On the command line, enter the following...

for unix and macOSX:
    ./mkplugin -str CaeMyExporter
    ./mkplugin -uns CaeMyExporter

for windows:
    mkplugin -str CaeMyExporter
    mkplugin -uns CaeMyExporter

    The -str/-uns switch indicates we are creating a structured/unstructured 
    plugin project.

If all goes well, you should see a message similar to one of the following:

    "creating new unstructured plugin 'CaeMyExporter' with ID=20..."
    "creating new structured plugin 'CaeMyExporter' with ID=20..."

To see the full, command line usage, use the -h switch as:
    ./mkplugin -h    (unix and macOSX)
    mkplugin -h      (windows)


*******************************************************
*** Manual Steps to Create a New CAE Plugin Project ***
*******************************************************

Create the project folder as:
    "/path/to/your/PluginSDK/src/plugins/CaeMyExporter/".

Copy the the following project files over to your new project folder
(do NOT copy the folders, just the files):

	.../src/plugins/templates/CAEP/*.*
	.../src/plugins/templates/PWP/*.*

--------------------

Make the following changes to the files in your new project folder. For the
steps below, I will use the folder name "CaeMyExporter".

--------------------

Rename the file CaeMyExporter/CaeTemplate.vcproj to
CaeMyExporter/CaeMyExporter.vcproj.

--------------------

Then, with a text editor, open the CaeMyExporter.vcproj file and replace all
string occurrences according to the following table.

Replace all "CaeTemplate" substrings with the project folder name. The case must
match *exactly*. In this example, replace "CaeTemplate" with "CaeMyExporter".

For a STRUCTURED plugin, replace all "plugintype" substrings with the string 
"structured". The case must match *exactly*.

For an UNSTRUCTURED plugin, replace all "plugintype" substrings with the 
string "unstructured". The case must match *exactly*.

Save CaeMyExporter.vcproj.

--------------------

Next, with a text editor, open the CaeMyExporter/module.mk file.

Replace all "CaeXxxxx" substrings with the project folder name. The case must
match *exactly*. In this example, replace "CaeXxxxx" with "CaeMyExporter".

Replace all "CAEXXXXX" substrings with the project folder name in UPPERCASE.
In this example, replace "CAEXXXXX" with "CAEMYEXPORTER".

Save module.mk.

--------------------

Next, with a text editor, open the CaeMyExporter/rtCaepInitItems.h file.

Replace all "CaeXxxxx" substrings with the project folder name. The case must
match *exactly*. In this example, replace "CaeXxxxx" with "CaeMyExporter".

Replace the "__ID__" substring with the plugin's unique ID number. This ID 
*must* be an integer value that is unique for all plugins in your build 
environment. This ID will also be added to the pluginRegistry.h file below. 
In this example, replace "__ID__" with "20".

Save rtCaepInitItems.h.

NOTE: You will need to make additional changes to this file to properly 
implement your plugin. Those changes are *not* covered here.

--------------------

Next, with a text editor, open the CaeMyExporter/rtCaepSupportData.h file.

Replace all "CaeXxxxx" substrings with the project folder name. The case must
match *exactly*. In this example, replace "CaeXxxxx" with "CaeMyExporter".

Save rtCaepSupportData.h.

NOTE: You will need to make additional changes to this file to properly 
implement your plugin. Those changes are *not* covered here.

--------------------

Next, with a text editor, open the .../plugins/modulelocal.mk file.

Add the CaeMyExporter target to the target-list defined in ALL_PLUGIN_MODULES.

If you do *not* want to build the example XML plugins shipped with the SDK, 
you may remove them from the target-list.

Save modulelocal.mk.

--------------------

Next, with a text editor, open the .../plugins/pluginRegistry.h file.

Add a macro defining the plugin's unique ID value. This ID *must* be unique 
for all plugins in your build environment. Be sure to place the new macro at 
the end of the list *above* the appropriate comment. For example,

    #define ID_CaeStrXML  9
    #define ID_CaeUnsXML 10
    #define ID_CaeMyExporter 20
    // add new plugin ID above this line

Save pluginRegistry.h.

--------------------

The End (for now)
