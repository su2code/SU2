---
title: Adding a Config Option
permalink: /docs_v7/Adding-a-Config-Option/
---

This page shows how to add an option to the SU2 configuration file, how to document it, and how to rename or remove
one. The examples are taken from the code.

The class `CConfig` reads the configuration file. An option touches these files:

| File | What goes there |
|---|---|
| `Common/include/CConfig.hpp` | The member that stores the value, and its getter |
| `Common/src/CConfig.cpp` | The registration in `CConfig::SetConfig_Options`, its comment, and the checks in `CConfig::SetPostprocessing` |
| `Common/include/option_structure.hpp` | The enum and its map, for options that take a keyword |
| `config_template.cfg` | An example line with a short explanation |

The [Option Reference](../Option-Reference/) of this website is generated from `Common/src/CConfig.cpp`,
`Common/include/option_structure.hpp` and `config_template.cfg` by `Docs/website/tools/build_site.py`. What you
write in the code is what users read on the website.

`CConfig.hpp`, `option_structure.hpp` and `CConfig.cpp` are excluded from clang-format (see
`.pre-commit-config.yaml`), so follow the layout of the code around your change.

## 1. Add the member and the getter

Declare the member in the private part of `CConfig` in `Common/include/CConfig.hpp`, with a short Doxygen comment:

```cpp
  su2double hRoughness;             /*!< \brief RMS roughness for Transition model. */
```

Add a `const` getter in the public part:

```cpp
  /*!
   * \brief Get RMS roughness for Transtion model from config
   * \return Value of roughness.
   */
  su2double GethRoughness(void) const { return hRoughness; }
```

The rest of the code reads the value through the getter, e.g. `config->GethRoughness()`.
See the [Code Reference](../Code-Reference/) for the existing members of `CConfig`.

## 2. Register the option

Add one `add...Option` call in `CConfig::SetConfig_Options` in `Common/src/CConfig.cpp`, in the group where the
option belongs. The first argument is the name used in the configuration file, in capital letters. The second is
the member. Most calls take a default value as the last argument:

```cpp
  /*!\brief HROUGHNESS \n DESCRIPTION: Value of RMS roughness for transition model \n DEFAULT: 1E-6 \ingroup Config*/
  addDoubleOption("HROUGHNESS", hRoughness, 1e-6);
```

The name must be unique. Each `add...Option` function checks this with an `assert`, which is active only in a debug
build (`--buildtype=debug`).

### Kinds of options

The functions are declared in `Common/include/CConfig.hpp`. The common ones are:

| Function | Arguments after the name | Example in `CConfig.cpp` |
|---|---|---|
| `addDoubleOption` | member (`su2double`), default | `addDoubleOption("HROUGHNESS", hRoughness, 1e-6);` |
| `addBoolOption` | member (`bool`), default | `addBoolOption("WRT_RESTART_COMPACT", Wrt_Restart_Compact, true);` |
| `addIntegerOption`, `addLongOption`, `addUnsignedLongOption`, `addUnsignedShortOption` | member, default | `addUnsignedLongOption("SCREEN_WRT_FREQ_INNER", ScreenWrtFreq[2], 1);` |
| `addStringOption` | member (`string`), default | `addStringOption("FLUID_NAME", FluidName, string("nitrogen"));` |
| `addEnumOption` | member, map, default | `addEnumOption("KIND_TURB_MODEL", Kind_Turb_Model, Turb_Model_Map, TURB_MODEL::NONE);` |
| `addEnumListOption` | size, pointer, map | `addEnumListOption("SST_OPTIONS", nSST_Options, SST_Options, SST_Options_Map);` |
| `addDoubleListOption`, `addShortListOption`, `addUShortListOption`, `addULongListOption`, `addStringListOption` | size, pointer | `addDoubleListOption("PRANDTL_LAM", nPrandtl_Lam , Prandtl_Lam);` |
| `addDoubleArrayOption`, `addUShortArrayOption`, `addULongArrayOption` | fixed size, `allow_fewer`, array | `addDoubleArrayOption("BODY_FORCE_VECTOR", 3, false, body_force);` |
| `addPythonOption` | none | `addPythonOption("OPT_OBJECTIVE");` |

In the configuration file, the values of lists and arrays are separated by spaces or commas and may be put in
parentheses, e.g. `SST_OPTIONS= (V2003m, VORTICITY)`.

Notes on the kinds:

- **Bool** options accept `YES` and `NO`.
- **Enum** options accept one keyword of the map. The member can be the enum type. See [Add an enum](#4-add-an-enum).
- **Lists** have any length. The function sets the size member to the number of values and allocates the array; the
  option object frees it. `NONE` gives an empty list. **Lists have no default**: if the option is not in the file,
  the size is 0. Set a default in `SetPostprocessing` if you need one (see below).
- **Arrays** have a fixed size. With `allow_fewer = true`, fewer values are accepted. **Arrays have no default
  argument**: fill the array before the call, as the code does for `BODY_FORCE_VECTOR`:

  ```cpp
  body_force[0] = 0.0; body_force[1] = 0.0; body_force[2] = 0.0;
  /* DESCRIPTION: Vector of body force values (BodyForce_X, BodyForce_Y, BodyForce_Z) */
  addDoubleArrayOption("BODY_FORCE_VECTOR", 3, false, body_force);
  ```

- **Python** options are read only by the Python scripts in `SU2_PY`. The C++ code accepts them without checking
  them, and has no member for them.

Options with a special format (boundary conditions, design variables, FFD boxes, ...) have their own functions,
such as `addInletOption`, `addRiemannOption`, `addPeriodicOption` and `addDVParamOption`. Their parsers are the
`COption...` classes in `Common/include/option_structure.inl`. Look at an existing option of the same kind before
you write a new parser.

### Multizone problems

In a multizone problem each zone has its own configuration file. An option that is not set in the file of a zone
takes the value of the main configuration file (`CConfig::SetDefaultFromConfig`), except `SCREEN_OUTPUT` and
`HISTORY_OUTPUT`. Options that are not set anywhere take their default value.

## 3. Write the comment

Put a Doxygen comment on the line above the `add...Option` call, in this form:

```cpp
/*!\brief NAME \n DESCRIPTION: what the option does \n Options: allowed values \n DEFAULT: value \ingroup Config*/
```

Real examples:

```cpp
  /*!\brief REYNOLDS_NUMBER \n DESCRIPTION: Reynolds number (non-dimensional, based on the free-stream values). Needed for viscous solvers. For incompressible solvers the Reynolds length will always be 1.0 \n DEFAULT: 0.0 \ingroup Config */
  addDoubleOption("REYNOLDS_NUMBER", Reynolds, 0.0);

  /*!\brief KIND_TURB_MODEL \n DESCRIPTION: Specify turbulence model \n Options: see \link Turb_Model_Map \endlink \n DEFAULT: NONE \ingroup Config*/
  addEnumOption("KIND_TURB_MODEL", Kind_Turb_Model, Turb_Model_Map, TURB_MODEL::NONE);
```

`build_site.py` reads the comment as follows. Follow these rules, or the Option Reference will be incomplete:

- The comment must start with `/*!`. A comment that starts with `/*` (such as `/* DESCRIPTION: ... */`) is not read,
  and the option gets "No description in the code yet".
- The comment must be the last `/*!` comment before the call.
- The parts are separated by `\n`. The keys are `DESCRIPTION:`, `Options:` (or `OPTIONS:`) and `DEFAULT:`, each with
  a colon. `\link Name \endlink` is shown as `Name`.
- The description should say what the option does, the unit if there is one, and which solvers use it.

The website also adds, without any comment:

- the type, from the `add...Option` function;
- for enum options, all the keywords of the map;
- the default, from the last argument of the call (for an enum, the keyword of the map that has this value);
- the line of the option in `config_template.cfg`, or "Not in config_template.cfg";
- a link to the line of the call in `CConfig.cpp`.

### Categories

The options are grouped into pages by the category comments in `SetConfig_Options`, for example:

```cpp
  /*!\par CONFIG_CATEGORY: Freestream Conditions \ingroup Config*/
```

Every option after such a comment belongs to that category, until the next one. Put a new option after the category
where it belongs. Add a new category only for a new group of options.

## 4. Add an enum

An option that takes a keyword needs an enum and a map in `Common/include/option_structure.hpp`. Use an
`enum class`. This is the enum of `VOLUME_OUTPUT_PRECISION`:

```cpp
/*!
 * \brief Floating-point type of the fields of the volume and surface files (CGNS and Paraview XML).
 */
enum class VOLUME_OUTPUT_PRECISION {
  SINGLE,             /*!< \brief Single precision (float). */
  DOUBLE              /*!< \brief Double precision (double). */
};
static const MapType<std::string, VOLUME_OUTPUT_PRECISION> VolumeOutputPrecision_Map = {
  MakePair("SINGLE", VOLUME_OUTPUT_PRECISION::SINGLE)
  MakePair("DOUBLE", VOLUME_OUTPUT_PRECISION::DOUBLE)
};
```

The member in `CConfig.hpp`, the registration and the getter:

```cpp
  VOLUME_OUTPUT_PRECISION Volume_Output_Precision; /*!< \brief Floating-point type of the volume and surface fields. */
```

```cpp
  addEnumOption("VOLUME_OUTPUT_PRECISION", Volume_Output_Precision, VolumeOutputPrecision_Map,
                VOLUME_OUTPUT_PRECISION::SINGLE);
```

```cpp
  bool GetVolume_Output_Double_Precision(void) const {
    return Volume_Output_Precision == VOLUME_OUTPUT_PRECISION::DOUBLE;
  }
```

Many getters of enum options return the enum itself, e.g. `GetKind_Turb_Model()`.

Things to know about the maps:

- `MakePair` is a macro that adds the comma itself. Do not put commas between the `MakePair` lines.
- The maps are real `std::map` objects only where `ENABLE_MAPS` is defined before `CConfig.hpp` is included
  (`Common/src/CConfig.cpp` and `SU2_DEF/src/drivers/CDiscAdjDeformationDriver.cpp`). Everywhere else
  `MapType` is the empty class `CEmptyMap`, to make compiling faster. Use the enum in the rest of the code, not the map.
- `build_site.py` finds the keywords with the pattern `static const MapType<std::string, TYPE> Name = { MakePair(...) ... };`
  in `Common/include`. Keep this form. Write the default of `addEnumOption` exactly as in the map
  (`VOLUME_OUTPUT_PRECISION::SINGLE`), so the website can show it as a keyword.

## 5. Check the value and set derived values

The constructors of `CConfig` call, among other steps and in this order: `SetConfig_Parsing` (reads the values),
`SetDefault` (sets the default of the options that were not in the file), `SetPostprocessing` (checks and derived
values) and `SetMarkers` (boundary markers). Most constructors then call `SetOutput`, which prints the summary of
the settings at start-up.

Put the checks of the new option in `CConfig::SetPostprocessing`. Stop with a clear message that names the options
and the allowed values:

```cpp
  if (Kind_Turb_Model == TURB_MODEL::NONE && Kind_Trans_Model != TURB_TRANS_MODEL::NONE) {
    SU2_MPI::Error("KIND_TURB_MODEL cannot be NONE to use a transition model", CURRENT_FUNCTION);
  }
```

`OptionIsSet("NAME")` tells if the option was in the configuration file. Use it when the default depends on other
options, or to change a value silently only when the user did not ask for it:

```cpp
  if (MUSCL_Flow && (Kind_ConvNumScheme_Flow == SPACE_CENTERED)) {
    if (OptionIsSet("MUSCL_FLOW")) {
      SU2_MPI::Error("Centered schemes do not use MUSCL reconstruction (use MUSCL_FLOW= NO).", CURRENT_FUNCTION);
    } else {
      MUSCL_Flow = false;
    }
  }
```

Lists get their defaults here. For example, the lists of fluid properties use a small helper:

```cpp
  auto SetDefaultIfEmpty = [](su2double*& array, unsigned short& size, const su2double& default_val) {
    if (array == nullptr) {
      array = new su2double[1];
      array[0] = default_val;
      size = 1;
    }
  };
```

```cpp
  SetDefaultIfEmpty(Prandtl_Lam, nPrandtl_Lam, Prandtl_Lam_Default);
```

Options that name boundary markers are also processed in `CConfig::SetMarkers`. Follow an existing marker option,
such as `MARKER_PLOTTING`.

## 6. Add the option to config_template.cfg

Add the option to `config_template.cfg` in the section where it belongs. Write one or more comment lines with what
it does and the allowed values, then the option with its default value:

```
% Value of RMS roughness for transition model
HROUGHNESS= 1.0e-6
```

The Option Reference shows this line, and counts the options that are missing from the template. The checklist of
the pull request template also asks for it.

## Rename or remove an option

SU2 stops when the configuration file has an option that does not exist. The message names the line and the option
and suggests the closest existing name:

```
Line 12 NAME: invalid option name. Check current SU2 options in config_template.cfg.
Did you mean ...?
```

For a renamed or removed option, `CConfig::SetConfig_Parsing` can print a better message. It has a chain of `if`
blocks with the old names. Add the old name there, with the new name if there is one:

```cpp
          else if (!option_name.compare("DYN_TIMESTEP"))
            newString.append("DYN_TIMESTEP is deprecated. Use TIME_STEP instead.\n\n");
          else if (!option_name.compare("DYN_TIME"))
            newString.append("DYN_TIME is deprecated. Use MAX_TIME instead.\n\n");
```

For a removed keyword of an enum option, there is a second chain after the call to `SetValue`. It runs only when the
value is not valid, e.g. for `KIND_TURB_MODEL`:

```cpp
        if (!option_name.compare("KIND_TURB_MODEL")) {
          if (option_value[0] == "SST_SUST")
            errorString.append("Option KIND_TURB_MODEL=SST_SUST is deprecated. Use KIND_TURB_MODEL=SST, SST_OPTIONS=SUSTAINING instead.\n");
```

These are the only mechanisms. There is no automatic translation of old names: the run stops with the message.

When you rename or remove an option, also:

1. Update `config_template.cfg`.
2. Update every configuration file in `TestCases/` that uses it (search with `git grep`). The TestCases and
   Tutorials repositories can have files that use it too.
3. Update the website pages and tutorials that name it. The [Docs Check](../Docs-Check/) page lists the options named
   in the docs that do not exist in the code.

## Python scripts

The C++ option needs no change in Python. `read_config` in `SU2_PY/SU2/io/config.py` stores every option it does not
know as a string. Add a case there only if a Python script must read the value as a number, a list or a dictionary.
The Python wrapper (`SU2_PY/pySU2/pySU2.i`) wraps the drivers, not `CConfig`.

## Test the option

Add a regression test that uses the new option, or add the option to the configuration file of an existing test.
See [Running Regression Tests](../Running-Regression-Tests/).

## Checklist

1. Member and getter in `Common/include/CConfig.hpp`.
2. `add...Option` call with a `/*!\brief ... \ingroup Config*/` comment in `CConfig::SetConfig_Options`.
3. Enum and map in `Common/include/option_structure.hpp`, if the option takes a keyword.
4. Checks and defaults in `CConfig::SetPostprocessing`.
5. Example in `config_template.cfg`.
6. A regression test that uses the option.

Last checked against SU2 develop `11c177b904`.
