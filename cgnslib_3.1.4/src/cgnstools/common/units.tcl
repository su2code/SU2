
#----- conversions to SI units

# mass -> kg

proc convert_mass {val unit {dir from}} {
  switch -- $unit {
    kg      {set scl 1.0}
    lbm     {set scl 0.45359237}
    slug    {set scl 14.593881}
    g       {set scl 0.001}
    oz      {set scl 0.028349523}
    ton     {set scl 907.18474}
    default {return -code error "unknown mass conversion - $unit"}
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

# length -> m

proc convert_length {val unit {dir from}} {
  switch -- $unit {
    m       {set scl 1.0}
    ft      {set scl 0.3048}
    cm      {set scl 0.01}
    mm      {set scl 0.001}
    in      {set scl 0.0254}
    yd      {set scl 0.9144}
    km      {set scl 1000.0}
    mi      {set scl 1609.344}
    micron  {set scl 1.0e-06}
    thou    -
    mil     {set scl 2.54e-05}
    default {return -code error "unknown length conversion - $unit"}
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

# time -> s

proc convert_time {val unit {dir from}} {
  switch -- $unit {
    sec     -
    s       {set scl 1.0}
    min     {set scl 60.0}
    hr      {set scl 3600.0}
    default {return -code error "unknown time conversion - $unit"}
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

# temperature -> K

proc convert_temperature {val unit {dir from}} {
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    switch -- $unit {
      K {return [expr $exp]}
      C {return [expr ($exp) + 273.16]}
      R {return [expr ($exp) / 1.8]}
      F {return [expr (($exp) - 32.0) / 1.8 + 273.16]}
    }
  } else {
    switch -- $unit {
      K {return [expr $exp]}
      C {return [expr ($exp) - 273.16]}
      R {return [expr ($exp) * 1.8]}
      F {return [expr (($exp) - 273.16) * 1.8 + 32.0]}
    }
  }
  return -code error "unknown temperature conversion - $unit"
}

# angle -> rad

proc convert_angle {val unit {dir from}} {
  switch -- $unit {
    rad     {set scl 1.0}
    deg     {set scl 0.017453293}
    rev     {set scl 6.2831853}
    grad    {set scl 0.01570796}
    default {return -code error "unknown angle conversion - $unit"}
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

# area -> m^2

proc convert_area {val unit {dir from}} {
  switch -- $unit {
    m^2     {set scl 1.0}
    ft^2    {set scl 0.09290304}
    cm^2    {set scl 1.0e-04}
    mm^2    {set scl 1.0e-06}
    in^2    {set scl 6.4516e-04}
    yd^2    {set scl 0.83612736}
    default {return -code error "unknown area conversion - $unit"}
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

# volume -> m^3

proc convert_volume {val unit {dir from}} {
  switch -- $unit {
    m^3     {set scl 1.0}
    ft^3    {set scl 0.028316847}
    cc      -
    ml      -
    cm^3    {set scl 1.0e-06}
    mm^3    {set scl 1.0e-09}
    in^3    {set scl 1.6387064e-05}
    yd^3    {set scl 0.76455486}
    liter   -
    l       {set scl 0.001}
    gal     {set scl 0.003785412}
    cup     {set scl 2.3658825e-04}
    qt      {set scl 9.46353e-04}
    default {return -code error "unknown volume conversion - $unit"}
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

# velocity -> m/s

proc convert_velocity {val unit {dir from}} {
  switch -- $unit {
    m/s     {set scl 1.0}
    ft/s    {set scl 0.3048}
    cm/s    {set scl 0.01}
    mm/s    {set scl 0.001}
    in/s    {set scl 0.0254}
    km/h    {set scl 0.27777778}
    mi/h    {set scl 0.44704}
    mach    {set scl 340.29}
    default {return -code error "unknown velocity conversion - $unit"}
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

# rotation rate -> rad/s

proc convert_rotation {val unit {dir from}} {
  switch -- $unit {
    rad/s   {set scl 1.0}
    deg/s   {set scl 0.017453293}
    rev/s   {set scl 6.2831853}
    rad/min {set scl 0.016666667}
    deg/min {set scl 2.9088821e-04}
    rev/min -
    RPM     {set scl 0.10471976}
    default {return -code error "unknown rotation rate conversion - $unit"}
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

# density -> kg/m^3

proc convert_density {val unit {dir from}} {
  switch -- $unit {
    kg/m^3    {set scl 1.0}
    lbm/ft^3  {set scl 16.018463}
    slug/ft^3 {set scl 515.37804}
    g/cm^3    {set scl 1000.0}
    default   {return -code error "unknown density conversion - $unit"}
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

# pressure -> N/m^2 (Pa)

proc convert_pressure {val unit {dir from}} {
  switch -- $unit {
    N/m^2    -
    Pa       {set scl 1.0}
    MPa      {set scl 1.0e+06}
    lbf/in^2 -
    psi      {set scl 6894.7572}
    lbf/ft^2 -
    psf      {set scl 47.880258}
    bar      {set scl 1.0e+05}
    atm      {set scl 101325.0}
    dyn/cm^2 {set scl 0.1}
    inHg     {set scl 3386.388}
    mmHg     -
    torr     {set scl 133.3224}
    default  {return -code error "unknown pressure conversion - $unit"}
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

# enthalpy -> J/kg

proc convert_enthalpy {val unit {dir from}} {
  switch -- $unit {
    J/kg        {set scl 1.0}
    kJ/kg       {set scl 1000.0}
    BTU/lbm     {set scl 2325.9723}
    BTU/slug    {set scl 72.293539}
    ft-lbf/lbm  {set scl 2.9890669}
    ft-lbf/slug {set scl 0.09290318}
    default     {return -code error "unknown enthalpy conversion - $unit"}
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

# force -> N

proc convert_force {val unit {dir from}} {
  switch -- $unit {
    N       {set scl 1.0}
    lbf     {set scl 4.4482216}
    pdl     {set scl 0.13825516}
    dyn     {set scl 1.0e-05}
    default {return -code error "unknown force conversion - $unit"}
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

# energy/work/torque -> J

proc convert_energy {val unit {dir from}} {
  switch -- $unit {
    N-m     -
    m-N     -
    J       {set scl 1.0}
    kJ      {set scl 1000.0}
    ft-lbf  {set scl 1.3558179}
    ft-pdl  {set scl 0.042140172}
    BTU     {set scl 1055.0433}
    cal     {set scl 4.184}
    kcal    {set scl 4184.0}
    erg     {set scl 1.0e-07}
    default {return -code error "unknown energy conversion - $unit"}
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

proc convert_work {val unit {dir from}} {
  return [convert_energy $val $unit $dir]
}

proc convert_torque {val unit {dir from}} {
  return [convert_energy $val $unit $dir]
}

# power -> W

proc convert_power {val unit {dir from}} {
  switch -- $unit {
    J/s      -
    W        {set scl 1.0}
    kW       {set scl 1000.0}
    ft-lbf/s {set scl 1.3558179}
    ft-pdl/s {set scl 0.042140172}
    BTU/s    {set scl 1055.0433}
    BTU/hr   {set scl 0.29306758}
    hp       {set scl 745.69987}
    default  {return -code error "unknown power conversion - $unit"}
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

# mass flow rate -> kg/s

proc convert_massflow {val unit {dir from}} {
  switch -- $unit {
    kg/s    {set scl 1.0}
    lbm/s   {set scl 0.45359237}
    slug/s  {set scl 14.593881}
    g/s     {set scl 0.001}
    kg/h    {set scl 2.7777778e-04}
    lbm/h   {set scl 1.2599788e-04}
    slug/h  {set scl 0.0040538558}
    default {return -code error "unknown mass flow rate conversion - $unit"}
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

# volume flow -> m^3/s

proc convert_volflow {val unit {dir from}} {
  switch -- $unit {
    m^3/s    {set scl 1.0}
    ft^3/s   {set scl 0.028316847}
    gal/s    {set scl 0.003785412}
    l/s      {set scl 0.001}
    ml/s     -
    cc/s     -
    cm^3/s   {set scl 1.0e-06}
    m^3/min  {set scl 0.016666667}
    ft^3/min {set scl 4.7194744e-04}
    gal/min  -
    gpm      {set scl 6.30902e-05}
    default  {return -code error "unknown volume flow conversion - $unit"}
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

# dynamic viscosity -> N-s/m^2

proc convert_dynvisc {val unit {dir from}} {
  switch -- $unit {
    kg/m-s     -
    N-s/m^2    {set scl 1.0}
    slug/ft-s  -
    lbf-s/ft^2 {set scl 47.880187}
    lbm/ft-s   {set scl 1.4881639}
    dyn-s/cm^2 -
    g/cm-s     -
    Poise      {set scl 0.1}
    default    {
      return -code error "unknown dynamic viscosity conversion - $unit"
    }
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

# kinematic viscosity -> m^2/s

proc convert_kinvisc {val unit {dir from}} {
  switch -- $unit {
    Stoke   -
    m^2/s   {set scl 1.0}
    ft^2/s  {set scl 0.09290304}
    default {
      return -code error "unknown kinematic viscosity conversion - $unit"
    }
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

# specific heat -> J/kg-K

proc convert_spheat {val unit {dir from}} {
  switch -- $unit {
    J/kg-K        {set scl 1.0}
    kJ/kg-K       {set scl 1000.0}
    erg/g-K       {set scl 0.0001}
    ft-lbf/lbm-R  {set scl 5.3803205}
    ft-lbf/slug-R {set scl 0.16722572}
    BTU/lbm-R     {set scl 4186.7502}
    BTU/slug-R    {set scl 130.12837}
    default       {
      return -code error "unknown specific heat conversion - $unit"
    }
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

# conductivity -> N/(s*K)

proc convert_conduct {val unit {dir from}} {
  switch -- $unit {
    J/m-s-K    -
    N/s-K      {set scl 1.0}
    lbf/s-R    {set scl 8.0067989}
    BTU/ft-s-R {set scl 6230.5706}
    cal/m-s-K  {set scl 4.184}
    dyn/s-K    {set scl 1.0e-05}
    default    {return -code error "unknown conductivity conversion - $unit"}
  }
  regsub -all / $val {*1.0/} exp
  if {$dir == "from"} {
    return [expr ($exp) * $scl]
  }
  return [expr ($exp) / $scl]
}

# SI units

proc unitsSI {type} {
  switch -- $type {
    mass        {return [list kg 1.0]}
    length      {return [list m 1.0]}
    time        {return [list s 1.0]}
    temperature {return [list K 1.0]}
    angle       {return [list deg 57.2958]}
    area        {return [list m^2 1.0]}
    volume      {return [list m^3 1.0]}
    velocity    {return [list m/s 1.0]}
    rotation    {return [list RPM 9.54930]}
    density     {return [list kg/m^3 1.0]}
    pressure    {return [list MPa 1.0e-06]}
    enthalpy    {return [list kJ/kg 0.001]}
    force       {return [list N 1.0]}
    energy      {return [list J 1.0]}
    work        {return [list J 1.0]}
    power       {return [list W 1.0]}
    torque      {return [list m-N 1.0]}
    massflow    {return [list kg/s 1.0]}
    volflow     {return [list m^3/s 1.0]}
    dynvisc     {return [list N-s/m^2 1.0]}
    kinvisc     {return [list m^2/s 1.0]}
    spheat      {return [list kJ/kg-K 0.001]}
    conduct     {return [list N/s-K 1.0]}
  }
  return [list $type 1.0]
}

# SI -> EE units

proc unitsEE {type} {
  switch -- $type {
    mass        {return [list lbm 2.20462]}
    length      {return [list ft 3.28084]}
    time        {return [list s 1.0]}
    temperature {return [list R 1.8]}
    angle       {return [list deg 57.2958]}
    area        {return [list ft^2 10.7639]}
    volume      {return [list ft^3 35.3147]}
    velocity    {return [list ft/s 3.28084]}
    rotation    {return [list RPM 9.54930]}
    density     {return [list lbm/ft^3 0.062428]}
    pressure    {return [list psi 0.000145038]}
    enthalpy    {return [list BTU/lbm 0.000429928]}
    force       {return [list lbf 0.224809]}
    energy      {return [list BTU 0.000947828]}
    work        {return [list ft-lbf 0.737562]}
    power       {return [list hp 0.00134102]}
    torque      {return [list ft-lbf 0.737562]}
    massflow    {return [list lbm/s 2.20462]}
    volflow     {return [list ft^3/s 35.3147]}
    dynvisc     {return [list lbf-s/ft^2 0.0208855]}
    kinvisc     {return [list ft^2/s 10.7639]}
    spheat      {return [list BTU/lbm-R 0.000238849]}
    conduct     {return [list lbf/s-R 0.124894]}
  }
  return [list $type 1.0]
}

# SI -> EG units

proc unitsEG {type} {
  switch -- $type {
    mass        {return [list slug 0.0685219]}
    length      {return [list ft 3.28084]}
    time        {return [list s 1.0]}
    temperature {return [list R 1.8]}
    angle       {return [list deg 57.2958]}
    area        {return [list ft^2 10.7639]}
    volume      {return [list ft^3 35.3147]}
    velocity    {return [list ft/s 3.28084]}
    rotation    {return [list RPM 9.54930]}
    density     {return [list slug/ft^3 0.00194032]}
    pressure    {return [list lbf/ft^2 0.0208854]}
    enthalpy    {return [list ft-lbf/slug 10.7639]}
    force       {return [list lbf 0.224809]}
    energy      {return [list ft-lbf 0.737562]}
    work        {return [list ft-lbf 0.737562]}
    power       {return [list ft-lbf/s 0.737562]}
    torque      {return [list ft-lbf 0.737562]}
    massflow    {return [list slug/s 0.0685219]}
    volflow     {return [list ft^3/s 35.3147]}
    dynvisc     {return [list lbf-s/ft^2 0.0208855]}
    kinvisc     {return [list ft^2/s 10.7639]}
    spheat      {return [list ft-lbf/slug-R 5.97994]}
    conduct     {return [list lbf/s-R 0.124894]}
  }
  return [list $type 1.0]
}

#---------- unit conversion panel ----------

set _UnitTypes {mass length time temperature angle area volume \
  velocity rotation density pressure enthalpy force energy \
  power massflow volflow dynvisc kinvisc spheat conduct}

array set _UnitList {
  mass        {kg lbm slug g oz ton}
  length      {m ft cm mm in yd km mi micron mil}
  time        {s min hr}
  temperature {K R C F}
  angle       {rad deg rev grad}
  area        {m^2 ft^2 cm^2 mm^2 in^2 yd^2}
  volume      {m^3 ft^3 cm^3 mm^3 in^3 yd^3 gal l ml cup qt}
  velocity    {m/s ft/s cm/s mm/s in/s km/h mi/h mach}
  rotation    {rad/s RPM deg/s rev/s rad/min deg/min}
  density     {kg/m^3 lbm/ft^3 slug/ft^3 g/cm^3}
  pressure    {Pa MPa psi psf bar atm lbf/in^2 lbf/ft^2 dyn/cm^2 \
               inHg mmHg torr}
  enthalpy    {J/kg kJ/kg BTU/lbm BTU/slug ft-lbf/lbm ft-lbf/slug}
  force       {N lbf pdl dyn}
  energy      {J kJ ft-lbf ft-pdl BTU m-N cal kcal erg}
  power       {W kW J/s ft-lbf/s ft-pdl/s BTU/s BTU/hr hp}
  massflow    {kg/s lbm/s slug/s g/s kg/h lbm/h slug/h}
  volflow     {m^3/s ft^3/s cm^3/s gal/s gpm l/s m^3/min ft^3/min}
  dynvisc     {N-s/m^2 lbf-s/ft^2 kg/m-s lbm/ft-s slug/ft-s Poise \
               dyn-s/cm^2 g/cm-s}
  kinvisc     {m^2/s ft^2/s Stoke}
  spheat      {J/kg-K kJ/kg-K ft-lbf/lbm-R ft-lbf/slug-R \
               BTU/lbm-R BTU/slug-R erg/g-K}
  conduct     {N/s-K lbf/s-R BTU/ft-s-R J/m-s-K cal/m-s-K dyn/s-K}
}

array set _UnitData {
  mass:label        Mass
  length:label      Length
  time:label        Time
  temperature:label Temperature
  angle:label       Angle
  area:label        Area
  volume:label      Volume
  velocity:label    Velocity
  rotation:label    Rotation
  density:label     Density
  pressure:label    Pressure
  enthalpy:label    Enthalpy
  force:label       Force
  energy:label      Energy/Work/Torque
  power:label       Power
  massflow:label    "Mass Flow Rate"
  volflow:label     "Volume Flow Rate"
  dynvisc:label     "Dynamic Viscosity"
  kinvisc:label     "Kinematic Viscosity"
  spheat:label      "Specific Heat"
  conduct:label     Conductivity
  mass:last         {0 1}
  length:last       {0 1}
  time:last         {0 0}
  temperature:last  {0 1}
  angle:last        {0 1}
  area:last         {0 1}
  volume:last       {0 1}
  velocity:last     {0 1}
  rotation:last     {0 1}
  density:last      {0 1}
  pressure:last     {0 2}
  enthalpy:last     {0 2}
  force:last        {0 1}
  energy:last       {0 2}
  power:last        {0 3}
  massflow:last     {0 1}
  volflow:last      {0 1}
  dynvisc:last      {0 1}
  kinvisc:last      {0 1}
  spheat:last       {0 2}
  conduct:last      {0 1}
  type              mass
  input:value       1
  input:unit        kg
  output:value      1
  output:unit       kg
  saved             ""
  regkey            "UnitConversions"
}

proc units_convert {{loc .}} {
  global _UnitData _UnitTypes
  set w .unit_conversions
  if {[winfo exists $w]} {
    wm deiconify $w
    raise $w
    focus $w.values.input.ent
    $w.values.input.ent selection range 0 end
    return
  }
  toplevel $w
  units:create $w
  catch {center_window $w $loc}
}

proc units_read {{key ""}} {
  global _UnitTypes _UnitData
  if {[info procs tclreg_get] == ""} {return}
  if {$key == ""} {
    set key $_UnitData(regkey)
  }
  foreach i $_UnitTypes {
    if {![catch {tclreg_get $key $i} last] &&
      [llength $last] == 2} {
      set _UnitData($i:last) $last
    }
  }
  foreach i {type saved} {
    if {![catch {tclreg_get $key $i} val] && $val != ""} {
      set _UnitData($i) $val
    }
  }
  foreach i {input output} {
    if {![catch {tclreg_get $key $i} data] &&
      [llength $data] == 2} {
      set _UnitData($i:value) [lindex $data 0]
      set _UnitData($i,unit) [lindex $data 1]
    }
  }
}

proc units_write {{key ""}} {
  global _UnitTypes _UnitData
  if {[info procs tclreg_set] == ""} {return}
  if {$key == ""} {
    set key $_UnitData(regkey)
  }
  foreach i $_UnitTypes {
    catch {tclreg_set $key $i $_UnitData($i:last)}
  }
  foreach i {type saved} {
    catch {tclreg_set $key $i $_UnitData($i)}
  }
  foreach i {input output} {
    catch {tclreg_set $key $i "$_UnitData($i:value) $_UnitData($i:unit)"}
  }
}

proc units:create {top} {
  global _UnitData _UnitTypes

  wm title $top "Unit Conversions"
  wm protocol $top WM_DELETE_WINDOW "units:destroy $top"
  wm resizable $top 0 0

  if {$top == "."} {
    set w ""
  } else {
    set w $top
  }

  set f [frame $w.type]
  pack $f -side top -anchor w -padx 5 -pady 5

  label $f.lab -text "Unit Type"
  pack $f.lab -side left
  menubutton $f.but -indicatoron 1 -width 20 -menu $f.but.menu \
    -relief raised -text Mass -padx 3 -pady 3 \
    -highlightthickness 1 -takefocus 1
  pack $f.but -side left -padx 5

  menu $f.but.menu -tearoff 0
  foreach i $_UnitTypes {
    $f.but.menu add command -label $_UnitData($i:label) \
      -command "units:change {$w} $i"
  }

  frame $w.units
  pack $w.units -side top -fill x -expand 1

  set f [frame $w.units.input]
  pack $f -side left -padx 5 -fill x -expand 1

  label $f.lab -text "From Units Of"
  pack $f.lab -side top -anchor w
  scrollbar $f.scroll -command "$f.list yview"
  pack $f.scroll -side right -fill y
  listbox $f.list -width 15 -height 5 -yscroll "$f.scroll set" \
    -selectmode browse -exportselection 0
  pack $f.list -side left -fill both -expand 1

  bind $f.list <<ListboxSelect>> "units:compute {$w}"

  set f [frame $w.units.output]
  pack $f -side right -padx 5 -fill x -expand 1

  label $f.lab -text "To Units Of"
  pack $f.lab -side top -anchor w
  scrollbar $f.scroll -command "$f.list yview"
  pack $f.scroll -side right -fill y
  listbox $f.list -width 15 -height 5 -yscroll "$f.scroll set" \
    -selectmode browse -exportselection 0
  pack $f.list -side top -fill both -expand 1

  bind $f.list <<ListboxSelect>> "units:compute {$w}"

  frame $w.values
  pack $w.values -side top -padx 5 -pady 5 -fill x -expand 1

  set f [frame $w.values.input]
  pack $f -side top -fill x -expand 1
  label $f.lab -text Input -width 6 -anchor w
  pack $f.lab -side left
  entry $f.ent -width 20 -textvariable _UnitData(input:value)
  pack $f.ent -side left -fill x -expand 1
  label $f.unit -width 9 -textvariable _UnitData(input:unit) -anchor w
  pack $f.unit -side left

  set f [frame $w.values.output]
  pack $f -side top -fill x -expand 1
  label $f.lab -text Output -width 6 -anchor w
  pack $f.lab -side left
  entry $f.ent -width 20 -textvariable _UnitData(output:value) \
    -state disabled -bg [$f.lab cget -bg]
  pack $f.ent -side left -fill x -expand 1
  label $f.unit -width 9 -textvariable _UnitData(output:unit) -anchor w
  pack $f.unit -side left

  set f [frame $w.but]
  pack $f -side top -padx 5 -pady 5 -fill x -expand 1

  button $f.sto -text STO -width 5 -highlightthickness 1 -pady 0 \
    -underline 0 -command "units:store {$w}"
  button $f.rcl -text RCL -width 5 -highlightthickness 1 -pady 0 \
    -underline 0 -command "units:recall {$w}"
  button $f.copy -text Copy -width 5 -highlightthickness 1 -pady 0 \
    -underline 0 -command "units:copy {$w}"
  button $f.help -text Help -width 5 -highlightthickness 1 -pady 0 \
    -underline 0 -command units_help
  if {$top == "."} {
    button $f.exit -text Exit -width 5 -highlightthickness 1 -pady 0 \
      -underline 0 -command "units:destroy $top"
  } else {
    button $f.exit -text Close -width 5 -highlightthickness 1 -pady 0 \
      -underline 4 -command "units:destroy $top"
  }
  pack $f.sto $f.rcl $f.copy $f.help $f.exit -side left -expand 1

  if {$_UnitData(saved) == ""} {$f.rcl configure -state disabled}
  if {[info procs units_help] == ""} {$f.help configure -state disabled}

  bind $top <Alt-s> "units:invoke $f.sto"
  bind $top <Alt-r> "units:invoke $f.rcl"
  bind $top <Alt-c> "units:invoke $f.copy"
  bind $top <Alt-h> "units:invoke $f.help"
  bind $top <Alt-e> "units:invoke $f.exit"

  units:set $w

  focus $w.values.input.ent
  $w.values.input.ent selection range 0 end

  trace variable _UnitData(input:value) w "units:update {$w}"
}

proc units:invoke {but} {
  if {[$but cget -state] == "normal"} {
    $but configure -relief sunken
    update idletasks
    after 250
    $but configure -relief raised
    $but invoke
  }
}

proc units:destroy {top} {
  if {$top == "."} {
    set w ""
  } else {
    set w $top
  }
  trace vdelete _UnitData(input:value) w "units:update {$w}"
  units:last $w
  destroy $top
}

proc units:last {w} {
  global _UnitData
  set input [$w.units.input.list curselection]
  set output [$w.units.output.list curselection]
  set _UnitData($_UnitData(type):last) [list $input $output]
}

proc units:set {w} {
  global _UnitData _UnitList
  set unit $_UnitData(type)
  $w.type.but configure -text $_UnitData($unit:label)
  set n 0
  foreach i {input output} {
    $w.units.$i.list delete 0 end
    foreach j $_UnitList($unit) {
      $w.units.$i.list insert end $j
    }
    set sel [lindex $_UnitData($unit:last) $n]
    if {$sel != ""} {
      $w.units.$i.list selection set $sel $sel
      $w.units.$i.list see $sel
    }
    incr n
  }
  units:compute $w
}

proc units:change {w unit} {
  global _UnitData
  units:last $w
  set _UnitData(type) $unit
  units:set $w
}

proc units:update {w name1 name2 op} {
  units:compute $w
}

proc units:compute {w} {
  global _UnitData
  foreach i {input output} {
    set sel [$w.units.$i.list curselection]
    if {$sel == {}} {
      set _UnitData($i:unit) ""
    } else {
      set _UnitData($i:unit) [$w.units.$i.list get $sel $sel]
    }
  }
  set conv convert_$_UnitData(type)
  if {[info commands $conv] == "" ||
    $_UnitData(input:unit) == "" ||
    $_UnitData(output:unit) == ""} {
    set _UnitData(output:value) ""
  } elseif {$_UnitData(input:unit) == $_UnitData(output:unit)} {
    regsub -all / $_UnitData(input:value) {*1.0/} exp
    if {[catch {expr double($exp)} output]} {
      set _UnitData(output:value) ""
    } else {
      set _UnitData(output:value) [format "%g" $output]
    }
  } else {
    if {[catch {$conv $_UnitData(input:value) $_UnitData(input:unit)} val] ||
      [catch {$conv $val $_UnitData(output:unit) to} output]} {
      set _UnitData(output:value) ""
    } else {
      set _UnitData(output:value) [format "%g" $output]
    }
  }
  if {$_UnitData(output:value) == ""} {
    $w.but.sto configure -state disabled
    $w.but.copy configure -state disabled
  } else {
    $w.but.sto configure -state normal
    $w.but.copy configure -state normal
  }
}

proc units:store {w} {
  global _UnitData
  set _UnitData(saved) $_UnitData(output:value)
  $w.but.rcl configure -state normal
}

proc units:recall {w} {
  global _UnitData
  set _UnitData(input:value) $_UnitData(saved)
  $w.values.input.ent selection clear
  $w.values.input.ent icursor end
}

proc units:copy {w} {
  global _UnitData
  clipboard clear
  clipboard append $_UnitData(output:value)
}

