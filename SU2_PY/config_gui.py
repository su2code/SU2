#!/usr/bin/env python

## \file config_gui.py
#  \brief _____________.
#  \author A. Aranake
#  \version 7.0.4 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
# 
# The SU2 Project is maintained by the SU2 Foundation 
# (http://su2foundation.org)
#
# Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# make print(*args) function available in PY2.6+, does'nt work on PY < 2.6
from __future__ import print_function

from parse_config import *
import wx, sys
import wx.lib.scrolledpanel as sp


# ctrldict is a dictionary
# keys of ctrldict are option categories
# Values of ctrldict are lists (optlist) of lists (option_data)
# option_data takes the form [option_name, option_value, StaticText,Control]

class YesNoBox(wx.CheckBox):
  '''The regular checkbox returns True or False, this one returns YES or NO to match the SU2 config format.'''
  
  def __init__(self, parent, *args, **kwargs):
    wx.CheckBox.__init__(self, parent, *args, **kwargs)

  def GetValue(self):
    if self.Value:
      return "YES"
    else:
      return "NO"

  def SetValue(self,stringval):
    if stringval=="YES":
      self.Value=True
    else:
      self.Value=False

class LabeledComboBox():
  '''Wrap a StaticText and TextCtrl into a single object'''

  def __init__(self,parent,option_name,txtlabel,option_default,option_values,option_type,option_description):
    self.label       = wx.StaticText(parent,label=txtlabel+", "+option_type)
    self.option_name = option_name
    self.control     = wx.ComboBox(parent,value=option_default,choices=option_values)
    self.control.SetSize((400,20))
    self.sizer       = wx.BoxSizer(wx.HORIZONTAL)
    self.sizer.Add(self.label,wx.EXPAND)
    self.sizer.AddSpacer(20)
    self.sizer.Add(self.control,wx.EXPAND)
    self.sizer.SetMinSize((400,20))

    self.option_default     = option_default
    self.option_type        = option_type
    self.option_description = option_description

  def GetSizer(self):
    return self.sizer

  def GetValue(self):
    return self.control.GetValue()

  def SetValue(self,val):
    self.control.SetValue(val)

  def SetDefaultValue(self):
    self.control.SetValue(self.option_default)

  def GetCtrl(self):
    return self.control

class LabeledTextCtrl():
  '''Wrap a StaticText and TextCtrl into a single object'''

  def __init__(self,parent,option_name,txtlabel,option_default,option_type):
    self.label       = wx.StaticText(parent,label=txtlabel+", "+option_type)
    self.option_name = option_name
    self.control     = wx.TextCtrl(parent,value=option_default)
    self.control.SetSize((400,20))
    self.sizer       = wx.BoxSizer(wx.HORIZONTAL)
    self.sizer.Add(self.label,wx.EXPAND)
    self.sizer.AddSpacer(20)
    self.sizer.Add(self.control,wx.EXPAND)
    self.sizer.SetMinSize((400,20))

    self.option_default = option_default
    self.option_type    = option_type

  def GetSizer(self):
    return self.sizer

  def GetValue(self):
    return self.control.GetValue()

  def SetValue(self,val):
    self.control.SetValue(val)

  def SetDefaultValue(self):
    self.control.SetValue(self.option_default)

  def GetCtrl(self):
    return self.control

class config_gui(wx.Frame):
  def __init__(self,option_data):

    # wx Initialization
    wx.Frame.__init__(self, None, title="SU2 config file editor")

    # Define the main sizers
    self.frame_sizer  = wx.BoxSizer(wx.HORIZONTAL)
    self.main_sizer   = wx.BoxSizer(wx.HORIZONTAL)
    self.left_sizer   = wx.BoxSizer(wx.VERTICAL)
    self.right_sizer  = wx.BoxSizer(wx.VERTICAL)

    # Use a scrolled panel on the right side
    self.main_panel   = wx.Panel(self)
    self.scroll_sizer = wx.BoxSizer(wx.VERTICAL)
    self.right_panel  = sp.ScrolledPanel(self.main_panel,size=(500,500))
    self.right_panel.SetupScrolling()

    # Left side - list of option categories
    self.list_ctrl = wx.ListCtrl(self.main_panel,style=wx.LC_REPORT|wx.BORDER_SUNKEN,size=(300,600))
    self.list_ctrl.InsertColumn(0, 'Option Category')

    bigfont = wx.Font(20,wx.MODERN,wx.NORMAL,wx.BOLD)

    # Read the option_data and build controls
    self.ctrldict  = {}
    self.optlabels = {}
    for j,category in enumerate(option_data):

      self.list_ctrl.InsertStringItem(j, category)  # Add category to left size list

      self.optlabels[category] = wx.StaticText(self.right_panel,label=category)
      self.optlabels[category].SetFont(bigfont)

      if j>0:
        self.scroll_sizer.AddSpacer(20)
      self.scroll_sizer.Add(self.optlabels[category],wx.EXPAND)

      self.ctrldict[category] = []
      yctr = 0
      for j,opt in enumerate(option_data[category]):
        if opt.option_type in ["EnumOption","MathProblem","SpecialOption","ConvectOption"]:
          self.ctrldict[category].append(LabeledComboBox(self.right_panel,opt.option_name,opt.option_name,opt.option_default,opt.option_values,opt.option_type,opt.option_description))
        else:
          self.ctrldict[category].append(LabeledTextCtrl(self.right_panel,opt.option_name,opt.option_name,opt.option_default,opt.option_type))

      for control in self.ctrldict[category]:
        self.scroll_sizer.Add(control.GetSizer(),wx.EXPAND)     # Add each control to the sizer
        self.lastctrl = control.GetCtrl()

    # Set right_panel to scroll vertically
    self.right_panel.SetSizer(self.scroll_sizer)

    # Set up menu
    menuBar = wx.MenuBar()
    m_file = wx.Menu()
    m_save = m_file.Append(wx.ID_SAVE, "&Save", "Save an SU2 .cfg file")
    m_open = m_file.Append(wx.ID_OPEN, "&Open", "Load an SU2 .cfg file")
    m_exit = m_file.Append(wx.ID_EXIT, "E&xit", "Close window and exit program.")

    menuBar.Append(m_file, "&File")
    self.SetMenuBar(menuBar)
    self.CreateStatusBar()
  
    # Specify which functions to call when stuff is changed
    self.Bind(wx.EVT_LIST_ITEM_SELECTED, self.list_click, self.list_ctrl)
    self.Bind(wx.EVT_MENU, self.OnSave, m_save)
    self.Bind(wx.EVT_MENU, self.OnOpen, m_open)
    self.Bind(wx.EVT_SIZE, self.OnResize)
    self.right_panel.SetAutoLayout(1)

    # Add it all to the panel and draw
    self.left_sizer.SetMinSize((300,600))
    self.right_sizer.SetMinSize((300,600))
    self.list_ctrl.SetColumnWidth(0,500)

    self.left_sizer.Add(self.list_ctrl,0,wx.EXPAND)
    self.right_sizer.Add(self.right_panel,0,wx.EXPAND)
    self.main_sizer.Add(self.left_sizer,0,wx.EXPAND)
    self.main_sizer.Add(self.right_sizer,0,wx.EXPAND)
    self.frame_sizer.Add(self.main_sizer,0,wx.EXPAND)

    self.main_panel.SetSizer(self.main_sizer)
    self.SetSizer(self.frame_sizer)
    self.SetInitialSize()

  def OnCheck(self,event):
    print(event)
    print(dir(event))
    print(event.GetEventObject())
    print(event.GetEventObject().GetValue())

  def OnResize(self,event):
  # There is surely a better way to do this....
    framesize = self.GetSize()
    self.main_panel.SetSize(framesize)
    self.list_ctrl.SetSize((300,framesize[1]-50))
    self.right_panel.SetSize((framesize[0]-300,framesize[1]-50))


  def OnSave(self,event):

    # Dialog to set output filename
    SaveDialog = wx.FileDialog(self,style=wx.FD_SAVE,wildcard='*.cfg')
    SaveDialog.ShowModal()
    outfile    = SaveDialog.GetPath()

    f = open(outfile,'w')

    for category,optlist in self.ctrldict.items():
      f.write('%% %s \n'%category) 
      for option in optlist:
        value = option.GetValue()
        if not value=="":
          f.write('%s=%s\n'%(option.option_name,option.GetValue()))
    f.close()

  def OnOpen(self,event):

    # Dialog to select file
    OpenDialog = wx.FileDialog(self,wildcard='*.cfg')
    OpenDialog.ShowModal()
    infile     = OpenDialog.GetPath()

    # Open file, Load values into a dictionary
    cfgfile    = open(infile,'r')
    infiledict = {}
    lines      = cfgfile.readlines()
    for line in lines: 
      if (line[0]=='%' or line.find('=')==-1): # Skip lines w/o = sign
        continue
      else:    # Get value
        key,val=[text.strip() for text in line.split('=')]
        key = key.upper()    # AoA should work as well as AOA
        infiledict[key] = val

    # Loop through controls and set them to the new values
    for category,optlist in self.ctrldict.items():
      for option in optlist:
        if option.option_name in infiledict:
          option.SetValue(infiledict[option.option_name])
        else:
          option.SetValue("")

  def OnClose(self,event):
    sys.exit(1)

  def OnAbout(self,event):
    print("OnAbout")

  def list_click(self, event):

    category = event.GetText()

    self.right_panel.ScrollChildIntoView(self.lastctrl)
    self.right_panel.ScrollChildIntoView(self.optlabels[category])

     
def prepare_data():
  """ Method to get configuration data from source files. 
      Outputs a dictionary of categories as keys and lists of config_options as values
  """

  # These variables should point to the configuration files
  su2_basedir = os.environ['SU2_HOME']
  config_cpp = os.path.join(su2_basedir,'Common/src/config_structure.cpp')
  config_hpp = os.path.join(su2_basedir,'Common/include/option_structure.hpp')

  # Check that files exist
  if not os.path.isfile(config_cpp):
    sys.exit('Could not find cpp file, please check that su2_basedir is set correctly in config_gui.py')
  if not os.path.isfile(config_hpp):
    sys.exit('Could not find hpp file, please check that su2_basedir is set correctly in config_gui.py')
 
  # Run the parser
  option_list = parse_config(config_cpp, config_hpp)

  # Organize data into dictionary with categories as keys
  option_data = {}
  for opt in option_list:
    if not opt.option_category in option_data:
      option_data[opt.option_category] = []
    option_data[opt.option_category].append(opt)

  return option_data

if __name__=="__main__":

  # Read source files, get option data
  option_data = prepare_data()     
  
  # Launch GUI
  app = wx.App(None)       
  frame = config_gui(option_data)  
  frame.Show()
  app.MainLoop()
