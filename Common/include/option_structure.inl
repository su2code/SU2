/*!
 * \file option_structure.inl
 * \brief Template derived classes from COption, defined here as we
 *        only include them where needed to reduce compilation time.
 * \author J. Hicken, B. Tracey
 * \version 7.0.3 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

template <class Tenum>
class COptionEnum : public COptionBase {

  map<string, Tenum> m;
  unsigned short & field; // Reference to the feildname
  Tenum def; // Default value
  string name; // identifier for the option

public:
  COptionEnum(string option_field_name, const map<string, Tenum> m, unsigned short & option_field, Tenum default_value) : field(option_field) {
    this->m = m;
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionEnum() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // Check if there is more than one string
    string out = optionCheckMultipleValues(option_value, "enum", this->name);
    if (out.compare("") != 0) {
      return out;
    }

    // Check to see if the enum value is in the map
    if (this->m.find(option_value[0]) == m.end()) {
      string str;
      str.append(this->name);
      str.append(": invalid option value ");
      str.append(option_value[0]);
      str.append(". Check current SU2 options in config_template.cfg.");
      return str;
    }
    // If it is there, set the option value
    Tenum val = this->m[option_value[0]];
    this->field = val;
    return "";
  }

  void SetDefault() {
    this->field = this->def;
  }
};

class COptionDouble : public COptionBase {
  su2double & field; // Reference to the fieldname
  su2double def; // Default value
  string name; // identifier for the option

public:
  COptionDouble(string option_field_name, su2double & option_field, su2double default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionDouble() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // check if there is more than one value
    string out = optionCheckMultipleValues(option_value, "su2double", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    istringstream is(option_value[0]);
    su2double val;
    if (is >> val) {
      this->field = val;
      return "";
    }
    return badValue(option_value, "su2double", this->name);
  }
  void SetDefault() {
    this->field = this->def;
  }
};

class COptionString : public COptionBase {
  string & field; // Reference to the fieldname
  string def; // Default value
  string name; // identifier for the option

public:
  COptionString(string option_field_name, string & option_field, string default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionString() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // check if there is more than one value
    string out = optionCheckMultipleValues(option_value, "su2double", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    this->field.assign(option_value[0]);
    return "";
  }
  void SetDefault() {
    this->field = this->def;
  }
};

class COptionInt : public COptionBase {
  int & field; // Reference to the feildname
  int def; // Default value
  string name; // identifier for the option

public:
  COptionInt(string option_field_name, int & option_field, int default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionInt() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    string out = optionCheckMultipleValues(option_value, "int", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    istringstream is(option_value[0]);
    int val;
    if (is >> val) {
      this->field = val;
      return "";
    }
    return badValue(option_value, "int", this->name);
  }
  void SetDefault() {
    this->field = this->def;
  }
};

class COptionULong : public COptionBase {
  unsigned long & field; // Reference to the feildname
  unsigned long def; // Default value
  string name; // identifier for the option

public:
  COptionULong(string option_field_name, unsigned long & option_field, unsigned long default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionULong() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    string out = optionCheckMultipleValues(option_value, "unsigned long", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    istringstream is(option_value[0]);
    unsigned long val;
    if (is >> val) {
      this->field = val;
      return "";
    }
    return badValue(option_value, "unsigned long", this->name);
  }
  void SetDefault() {
    this->field = this->def;
  }
};

class COptionUShort : public COptionBase {
  unsigned short & field; // Reference to the feildname
  unsigned short def; // Default value
  string name; // identifier for the option

public:
  COptionUShort(string option_field_name, unsigned short & option_field, unsigned short default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionUShort() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    string out = optionCheckMultipleValues(option_value, "unsigned short", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    istringstream is(option_value[0]);
    unsigned short val;
    if (is >> val) {
      this->field = val;
      return "";
    }
    return badValue(option_value, "unsigned short", this->name);
  }
  void SetDefault() {
    this->field = this->def;
  }
};

class COptionLong : public COptionBase {
  long & field; // Reference to the feildname
  long def; // Default value
  string name; // identifier for the option

public:
  COptionLong(string option_field_name, long & option_field, long default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionLong() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    string out = optionCheckMultipleValues(option_value, "long", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    istringstream is(option_value[0]);
    long val;
    if (is >> val) {
      this->field = val;
      return "";
    }
    return badValue(option_value, "long", this->name);
  }
  void SetDefault() {
    this->field = this->def;
  }
};

class COptionBool : public COptionBase {
  bool & field; // Reference to the feildname
  bool def; // Default value
  string name; // identifier for the option

public:
  COptionBool(string option_field_name, bool & option_field, bool default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionBool() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // check if there is more than one value
    string out = optionCheckMultipleValues(option_value, "bool", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    if (option_value[0].compare("YES") == 0) {
      this->field = true;
      return "";
    }
    if (option_value[0].compare("NO") == 0) {
      this->field = false;
      return "";
    }
    return badValue(option_value, "bool", this->name);
  }
  void SetDefault() {
    this->field = this->def;
  }
};

template <class Tenum>
class COptionEnumList : public COptionBase {

  map<string, Tenum> m;
  unsigned short * & field; // Reference to the feildname
  string name; // identifier for the option
  unsigned short & size;

public:
  COptionEnumList(string option_field_name, const map<string, Tenum> m, unsigned short * & option_field, unsigned short & list_size) : field(option_field) , size(list_size) {
    this->m = m;
    this->name = option_field_name;
  }

  ~COptionEnumList() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    if (option_value.size() == 1 && option_value[0].compare("NONE") == 0) {
      this->size = 0;
      return "";
    }
    // size is the length of the option list
    this->size = option_value.size();
    unsigned short * enums = new unsigned short[size];
    for (int i  = 0; i < this->size; i++) {
      // Check to see if the enum value is in the map
      if (this->m.find(option_value[i]) == m.end()) {
        string str;
        str.append(this->name);
        str.append(": invalid option value ");
        str.append(option_value[i]);
        str.append(". Check current SU2 options in config_template.cfg.");
        return str;
      }
      // If it is there, set the option value
      enums[i] = this->m[option_value[i]];
    }
    this->field = enums;
    return "";
  }

  void SetDefault() {
    // No default to set
    size = 0;
  }
};

class COptionDoubleArray : public COptionBase {
  su2double * & field; // Reference to the feildname
  string name; // identifier for the option
  const int size;
  su2double * def;
  su2double * vals;
  su2double * default_value;

public:
  COptionDoubleArray(string option_field_name, const int list_size, su2double * & option_field, su2double * default_value) : field(option_field), size(list_size) {
    this->name = option_field_name;
    this->default_value = default_value;
    def  = NULL;
    vals = NULL;
  }

  ~COptionDoubleArray() {
     if(def  != NULL) delete [] def;
     if(vals != NULL) delete [] vals;
  };
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // Check that the size is correct
    if (option_value.size() != (unsigned long)this->size) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": wrong number of arguments: ");
      stringstream ss;
      ss << this->size;
      newstring.append(ss.str());
      newstring.append(" expected, ");
      stringstream ss2;
      ss2 << option_value.size();
      newstring.append(ss2.str());
      newstring.append(" found");
      return newstring;
    }
    vals = new su2double[this->size];
    for (int i  = 0; i < this->size; i++) {
      istringstream is(option_value[i]);
      su2double val;
      if (!(is >> val)) {
        delete [] vals;
        return badValue(option_value, "su2double array", this->name);
      }
      vals[i] = val;
    }
    this->field = vals;
    return "";
  }

  void SetDefault() {
    def = new su2double [size];
    for (int i = 0; i < size; i++) {
      def[i] = default_value[i];
    }
    this->field = def;
  }
};

class COptionDoubleList : public COptionBase {
  su2double * & field; // Reference to the feildname
  string name; // identifier for the option
  unsigned short & size;

public:
  COptionDoubleList(string option_field_name, unsigned short & list_size, su2double * & option_field) : field(option_field), size(list_size) {
    this->name = option_field_name;
  }

  ~COptionDoubleList() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // The size is the length of option_value
    unsigned short option_size = option_value.size();
    if (option_size == 1 && option_value[0].compare("NONE") == 0) {
      // No options
      this->size = 0;
      return "";
    }

    this->size = option_size;

    // Parse all of the options
    su2double * vals = new su2double[option_size];
    for (unsigned long i  = 0; i < option_size; i++) {
      istringstream is(option_value[i]);
      su2double val;
      if (!(is >> val)) {
        delete [] vals;
        return badValue(option_value, "su2double list", this->name);
      }
      vals[i] = val;
    }
    this->field = vals;
    return "";
  }

  void SetDefault() {
    this->size = 0; // There is no default value for list
  }
};

class COptionShortList : public COptionBase {
  short * & field; // Reference to the feildname
  string name; // identifier for the option
  unsigned short & size;

public:
  COptionShortList(string option_field_name, unsigned short & list_size,  short * & option_field) : field(option_field), size(list_size) {
    this->name = option_field_name;
  }

  ~COptionShortList() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // The size is the length of option_value
    unsigned short option_size = option_value.size();
    if (option_size == 1 && option_value[0].compare("NONE") == 0) {
      // No options
      this->size = 0;
      return "";
    }
    this->size = option_size;

    // Parse all of the options
    short * vals = new  short[option_size];
    for (unsigned long i  = 0; i < option_size; i++) {
      istringstream is(option_value[i]);
      unsigned short val;
      if (!(is >> val)) {
        delete [] vals;
        return badValue(option_value, "short", this->name);
      }
      vals[i] = val;
    }
    this->field = vals;
    return "";
  }

  void SetDefault() {
    this->size = 0; // There is no default value for list
  }
};

class COptionUShortList : public COptionBase {
  unsigned short * & field; // Reference to the feildname
  string name; // identifier for the option
  unsigned short & size;

public:
  COptionUShortList(string option_field_name, unsigned short & list_size, unsigned short * & option_field) : field(option_field), size(list_size) {
    this->name = option_field_name;
  }

  ~COptionUShortList() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // The size is the length of option_value
    unsigned short option_size = option_value.size();
    if (option_size == 1 && option_value[0].compare("NONE") == 0) {
      // No options
      this->size = 0;
      return "";
    }
    this->size = option_size;

    // Parse all of the options
    unsigned short * vals = new unsigned short[option_size];
    for (unsigned long i  = 0; i < option_size; i++) {
      istringstream is(option_value[i]);
      unsigned short val;
      if (!(is >> val)) {
        delete [] vals;
        return badValue(option_value, "unsigned short", this->name);
      }
      vals[i] = val;
    }
    this->field = vals;
    return "";
  }

  void SetDefault() {
    this->size = 0; // There is no default value for list
  }
};

class COptionStringList : public COptionBase {
  string * & field; // Reference to the feildname
  string name; // identifier for the option
  unsigned short & size;

public:
  COptionStringList(string option_field_name, unsigned short & list_size, string * & option_field) : field(option_field), size(list_size) {
    this->name = option_field_name;
  }

  ~COptionStringList() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // The size is the length of option_value
    unsigned short option_size = option_value.size();
    if (option_size == 1 && option_value[0].compare("NONE") == 0) {
      this->size = 0;
      return "";
    }
    this->size = option_size;

    // Parse all of the options
    string * vals = new string[option_size];
    for (unsigned long i  = 0; i < option_size; i++) {
      vals[i].assign(option_value[i]);
    }
    this->field = vals;
    return "";
  }

  void SetDefault() {
    this->size = 0; // There is no default value for list
  }
};

class COptionConvect : public COptionBase {
  string name; // identifier for the option
  unsigned short & space;
  unsigned short & centered;
  unsigned short & upwind;

public:
  COptionConvect(string option_field_name, unsigned short & space_field, unsigned short & centered_field, unsigned short & upwind_field)
    : name(option_field_name), space(space_field), centered(centered_field), upwind(upwind_field) { }

  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);

    string out = optionCheckMultipleValues(option_value, "unsigned short", this->name);
    if (out.compare("") != 0) {
      return out;
    }

    if (Centered_Map.count(option_value[0])) {
      this->space = Space_Map.find("SPACE_CENTERED")->second;
      this->centered = Centered_Map.find(option_value[0])->second;
      this->upwind = NO_UPWIND;
      return "";
    }
    if (Upwind_Map.count(option_value[0])) {
      this->space = Space_Map.find("SPACE_UPWIND")->second;
      this->upwind = Upwind_Map.find(option_value[0])->second;
      this->centered = NO_CENTERED;
      return "";
    }
    // Make them defined in case something weird happens
    SetDefault();
    return badValue(option_value, "convect", this->name);

  }

  void SetDefault() {
    this->centered = NO_CENTERED;
    this->upwind = NO_UPWIND;
    this->space = NO_CONVECTIVE;
  }
};

class COptionFEMConvect : public COptionBase{
  string name; // identifier for the option
  unsigned short & space;
  unsigned short & fem;

public:
  COptionFEMConvect(string option_field_name, unsigned short & space_field, unsigned short & fem_field) : space(space_field), fem(fem_field) {
    this->name = option_field_name;
  }

  ~COptionFEMConvect() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);

    string out = optionCheckMultipleValues(option_value, "unsigned short", this->name);
    if (out.compare("") != 0) {
      return out;
    }

    if (FEM_Map.count(option_value[0])) {
      this->space = Space_Map.find("FINITE_ELEMENT")->second;
      this->fem = FEM_Map.find(option_value[0])->second;
      return "";
    }

    // Make them defined in case something weird happens
    this->fem = NO_FEM;
    return badValue(option_value, "convect", this->name);

  }

  void SetDefault() {
    this->fem = NO_FEM;
  }
};

class COptionMathProblem : public COptionBase {
  string name; // identifier for the option
  bool & cont_adjoint;
  bool cont_adjoint_def;
  bool & disc_adjoint;
  bool disc_adjoint_def;
  bool & restart;
  bool restart_def;

public:
  COptionMathProblem(string option_field_name, bool & cont_adjoint_field, bool cont_adjoint_default, bool & disc_adjoint_field, bool disc_adjoint_default, bool & restart_field, bool restart_default) : cont_adjoint(cont_adjoint_field), disc_adjoint(disc_adjoint_field), restart(restart_field) {
    this->name = option_field_name;
    this->cont_adjoint_def = cont_adjoint_default;
    this->disc_adjoint_def = disc_adjoint_default;
    this->restart_def = restart_default;
  }

  ~COptionMathProblem() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    string out = optionCheckMultipleValues(option_value, "unsigned short", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    if (option_value[0] == "ADJOINT") {
      return badValue(option_value, "math problem (try CONTINUOUS_ADJOINT)", this->name);
    }
    if (Math_Problem_Map.find(option_value[0]) == Math_Problem_Map.end()) {
      return badValue(option_value, "math problem", this->name);
    }
    if (option_value[0] == "DIRECT") {
      this->cont_adjoint = false;
      this->disc_adjoint = false;
      this->restart = false;
      return "";
    }
    if (option_value[0] == "CONTINUOUS_ADJOINT") {
      this->cont_adjoint= true;
      this->disc_adjoint = false;
      this->restart= true;
      return "";
    }
    if (option_value[0] == "DISCRETE_ADJOINT") {
      this->disc_adjoint = true;
      this->cont_adjoint= false;
      this->restart = true;
      return "";
    }
    return "option in math problem map not considered in constructor";
  }

  void SetDefault() {
    this->cont_adjoint = this->cont_adjoint_def;
    this->disc_adjoint = this->disc_adjoint_def;
    this->restart = this->restart_def;
  }

};

class COptionDVParam : public COptionBase {
  string name; // identifier for the option
  unsigned short & nDV;
  su2double ** & paramDV;
  string * & FFDTag;
  unsigned short* & design_variable;

public:
  COptionDVParam(string option_field_name, unsigned short & nDV_field, su2double** & paramDV_field, string* & FFDTag_field, unsigned short * & design_variable_field) : nDV(nDV_field), paramDV(paramDV_field), FFDTag(FFDTag_field), design_variable(design_variable_field) {
    this->name = option_field_name;
  }

  ~COptionDVParam() {};

  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    if ((option_value.size() == 1) && (option_value[0].compare("NONE") == 0)) {
      this->nDV = 0;
      return "";
    }

    // Cannot have ; at the beginning or the end
    if (option_value[0].compare(";") == 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": may not have beginning semicolon");
      return newstring;
    }
    if (option_value[option_value.size()-1].compare(";") == 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": may not have ending semicolon");
      return newstring;
    }


    // use the ";" token to determine the number of design variables
    // This works because semicolon is not one of the delimiters in tokenize string
    this->nDV = 0;
    //unsigned int num_semi = 0;
    for (unsigned int i = 0; i < static_cast<unsigned int>(option_value.size()); i++) {
      if (option_value[i].compare(";") == 0) {
        this->nDV++;
        //        num_semi++;
      }
    }

    // One more design variable than semicolon
    this->nDV++;

    if ( (this->nDV > 0) && (this->design_variable == NULL) ) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": Design_Variable array has not been allocated. Check that DV_KIND appears before DV_PARAM in configuration file.");
      return newstring;
    }

    this->paramDV = new su2double*[this->nDV];
    for (unsigned short iDV = 0; iDV < this->nDV; iDV++) {
      this->paramDV[iDV] = new su2double[MAX_PARAMETERS];
    }

    this->FFDTag = new string[this->nDV];

   vector<unsigned short> nParamDV(nDV, 0);
   unsigned short totalnParamDV = 0;
   stringstream ss;
   unsigned int i = 0;

    for (unsigned short iDV = 0; iDV < this->nDV; iDV++) {
      switch (this->design_variable[iDV]) {
        case NO_DEFORMATION:       nParamDV[iDV] = 0; break;
        case FFD_SETTING:          nParamDV[iDV] = 0; break;
        case FFD_CONTROL_POINT_2D: nParamDV[iDV] = 5; break;
        case FFD_CAMBER_2D:        nParamDV[iDV] = 2; break;
        case FFD_THICKNESS_2D:     nParamDV[iDV] = 2; break;
        case FFD_TWIST_2D:         nParamDV[iDV] = 3; break;
        case HICKS_HENNE:          nParamDV[iDV] = 2; break;
        case SURFACE_BUMP:         nParamDV[iDV] = 3; break;
        case CST:                  nParamDV[iDV] = 3; break;
        case ANGLE_OF_ATTACK:      nParamDV[iDV] = 1; break;
        case SCALE:                nParamDV[iDV] = 0; break;
        case TRANSLATION:          nParamDV[iDV] = 3; break;
        case ROTATION:             nParamDV[iDV] = 6; break;
        case NACA_4DIGITS:         nParamDV[iDV] = 3; break;
        case PARABOLIC:            nParamDV[iDV] = 2; break;
        case AIRFOIL:              nParamDV[iDV] = 2; break;
        case FFD_CONTROL_POINT:    nParamDV[iDV] = 7; break;
        case FFD_NACELLE:          nParamDV[iDV] = 6; break;
        case FFD_GULL:             nParamDV[iDV] = 2; break;
        case FFD_TWIST:            nParamDV[iDV] = 8; break;
        case FFD_ROTATION:         nParamDV[iDV] = 7; break;
        case FFD_CONTROL_SURFACE:  nParamDV[iDV] = 7; break;
        case FFD_CAMBER:           nParamDV[iDV] = 3; break;
        case FFD_THICKNESS:        nParamDV[iDV] = 3; break;
        case FFD_ANGLE_OF_ATTACK:  nParamDV[iDV] = 2; break;
        case SURFACE_FILE:         nParamDV[iDV] = 0; break;
        case DV_EFIELD:            nParamDV[iDV] = 2; break;
        case DV_YOUNG:             nParamDV[iDV] = 0; break;
        case DV_POISSON:           nParamDV[iDV] = 0; break;
        case DV_RHO:               nParamDV[iDV] = 0; break;
        case DV_RHO_DL:            nParamDV[iDV] = 0; break;
        case SCALE_GRID:           nParamDV[iDV] = 0; break;
        case TRANSLATE_GRID:       nParamDV[iDV] = 3; break;
        case ROTATE_GRID:          nParamDV[iDV] = 6; break;
        default : {
          string newstring;
          newstring.append(this->name);
          newstring.append(": undefined design variable type found in configuration file.");
          return newstring;
        }
      }
      totalnParamDV += nParamDV[iDV];
    }

    if (totalnParamDV > option_value.size()){
      SU2_MPI::Error("Wrong number of arguments for DV_PARAM!", CURRENT_FUNCTION);
    }

    for (unsigned short iDV = 0; iDV < this->nDV; iDV++) {
      for (unsigned short iParamDV = 0; iParamDV < nParamDV[iDV]; iParamDV++) {

        ss << option_value[i] << " ";

        if ((iParamDV == 0) &&
            ((this->design_variable[iDV] == NO_DEFORMATION) ||
             (this->design_variable[iDV] == FFD_SETTING) ||
             (this->design_variable[iDV] == FFD_ANGLE_OF_ATTACK)||
             (this->design_variable[iDV] == FFD_CONTROL_POINT_2D) ||
             (this->design_variable[iDV] == FFD_CAMBER_2D) ||
             (this->design_variable[iDV] == FFD_TWIST_2D) ||
             (this->design_variable[iDV] == FFD_THICKNESS_2D) ||
             (this->design_variable[iDV] == FFD_CONTROL_POINT) ||
             (this->design_variable[iDV] == FFD_NACELLE) ||
             (this->design_variable[iDV] == FFD_GULL) ||
             (this->design_variable[iDV] == FFD_TWIST) ||
             (this->design_variable[iDV] == FFD_ROTATION) ||
             (this->design_variable[iDV] == FFD_CONTROL_SURFACE) ||
             (this->design_variable[iDV] == FFD_CAMBER) ||
             (this->design_variable[iDV] == FFD_THICKNESS))) {
              ss >> this->FFDTag[iDV];
              this->paramDV[iDV][iParamDV] = 0;
            }
        else
          ss >> this->paramDV[iDV][iParamDV];

        i++;
      }
      if (iDV < (this->nDV-1)) {
        if (option_value[i].compare(";") != 0) {
          string newstring;
          newstring.append(this->name);
          newstring.append(": a design variable in the configuration file has the wrong number of parameters");
          return newstring;
        }
        i++;
      }
    }

    // Need to return something...
    return "";
  }

  void SetDefault() {
    this->nDV = 0;
    this->paramDV = NULL;
    this->FFDTag = NULL;
    // Don't mess with the Design_Variable because it's an input, not modified
  }
};

class COptionDVValue : public COptionBase {
  string name; // identifier for the option
  unsigned short* & nDV_Value;
  su2double ** & valueDV;
  unsigned short & nDV;
  su2double ** & paramDV;
  unsigned short* & design_variable;

public:
  COptionDVValue(string option_field_name, unsigned short* & nDVValue_field, su2double** & valueDV_field, unsigned short & nDV_field,  su2double** & paramDV_field, unsigned short * & design_variable_field) : nDV_Value(nDVValue_field), valueDV(valueDV_field), nDV(nDV_field), paramDV(paramDV_field), design_variable(design_variable_field) {
    this->name = option_field_name;
  }

  ~COptionDVValue() {};

  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    if ((option_value.size() == 1) && (option_value[0].compare("NONE") == 0)) {
      this->nDV_Value = NULL;
      return "";
    }

    if ( (this->nDV > 0) && (this->design_variable == NULL) ) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": Design_Variable array has not been allocated. Check that DV_KIND appears before DV_VALUE in configuration file.");
      return newstring;
    }
    if ( (this->nDV > 0) && (this->paramDV == NULL) ) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": Design_Parameter array has not been allocated. Check that DV_PARAM appears before DV_VALUE in configuration file.");
      return newstring;
    }

    this->valueDV = new su2double*[this->nDV];
    this->nDV_Value = new unsigned short[this->nDV];

    for (unsigned short iDV = 0; iDV < this->nDV; iDV++) {
      this->valueDV[iDV] = new su2double[3];
    }

    unsigned short nValueDV = 0;
    unsigned short totalnValueDV = 0;
    stringstream ss;
    unsigned int i = 0;
    for (unsigned short iDV = 0; iDV < this->nDV; iDV++) {
      switch (this->design_variable[iDV]) {
        case FFD_CONTROL_POINT:
          if((this->paramDV[iDV][4] == 0) &&
             (this->paramDV[iDV][5] == 0) &&
             (this->paramDV[iDV][6] == 0)) {
            nValueDV = 3;
          } else {
            nValueDV = 1;
          }
          break;
        case FFD_CONTROL_POINT_2D:
          if((this->paramDV[iDV][3] == 0) &&
             (this->paramDV[iDV][4] == 0)) {
            nValueDV = 2;
          } else {
            nValueDV = 1;
          }
          break;
        default :
          nValueDV = 1;
      }

      this->nDV_Value[iDV] = nValueDV;

      totalnValueDV += nValueDV;

      for (unsigned short iValueDV = 0; iValueDV < nValueDV; iValueDV++) {

        if (i >= option_value.size()) {
          string newstring;
          newstring.append(this->name);
          newstring.append(": DV_VALUE does not contain enough entries to match DV_KIND or DV_PARAM.");
          return newstring;
        }

        ss << option_value[i] << " ";

        ss >> this->valueDV[iDV][iValueDV];

        i++;
      }
    }

    if (i != totalnValueDV) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": a design variable in the configuration file has the wrong number of values");
      return newstring;
    }

    // Need to return something...
    return "";
  }

  void SetDefault() {
    this->nDV_Value = 0;
    this->valueDV = NULL;
    // Don't mess with the Design_Variable because it's an input, not modified
  }
};

class COptionFFDDef : public COptionBase {
  string name;
  unsigned short & nFFD;
  su2double ** & CoordFFD;
  string * & FFDTag;

public:
  COptionFFDDef(string option_field_name, unsigned short & nFFD_field, su2double** & coordFFD_field, string* & FFDTag_field) : nFFD(nFFD_field), CoordFFD(coordFFD_field), FFDTag(FFDTag_field) {
    this->name = option_field_name;
  }

  ~COptionFFDDef() {};

  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    if ((option_value.size() == 1) && (option_value[0].compare("NONE") == 0)) {
      this->nFFD = 0;
      return "";
    }

    // Cannot have ; at the beginning or the end
    if (option_value[0].compare(";") == 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": may not have beginning semicolon");
      return newstring;
    }
    if (option_value[option_value.size()-1].compare(";") == 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": may not have ending semicolon");
      return newstring;
    }


    // use the ";" token to determine the number of design variables
    // This works because semicolon is not one of the delimiters in tokenize string
    this->nFFD = 0;
    for (unsigned int i = 0; i < static_cast<unsigned int>(option_value.size()); i++) {
      if (option_value[i].compare(";") == 0) {
        this->nFFD++;
      }
    }

    // One more design variable than semicolon
    this->nFFD++;

    this->CoordFFD = new su2double*[this->nFFD];
    for (unsigned short iFFD = 0; iFFD < this->nFFD; iFFD++) {
      this->CoordFFD[iFFD] = new su2double[25];
    }

    this->FFDTag = new string[this->nFFD];

    unsigned short nCoordFFD = 0;
    stringstream ss;
    unsigned int i = 0;

    for (unsigned short iFFD = 0; iFFD < this->nFFD; iFFD++) {

      nCoordFFD = 25;

      for (unsigned short iCoordFFD = 0; iCoordFFD < nCoordFFD; iCoordFFD++) {

        ss << option_value[i] << " ";

        if (iCoordFFD == 0) ss >> this->FFDTag[iFFD];
        else ss >> this->CoordFFD[iFFD][iCoordFFD-1];

        i++;
      }

      if (iFFD < (this->nFFD-1)) {
        if (option_value[i].compare(";") != 0) {
          string newstring;
          newstring.append(this->name);
          newstring.append(": a FFD box in the configuration file has the wrong number of parameters");
          return newstring;
        }
        i++;
      }

    }

    // Need to return something...
    return "";
  }

  void SetDefault() {
    this->nFFD = 0;
    this->CoordFFD = NULL;
    this->FFDTag = NULL;
  }

};

class COptionFFDDegree : public COptionBase {
  string name;
  unsigned short & nFFD;
  unsigned short ** & DegreeFFD;

public:
  COptionFFDDegree(string option_field_name, unsigned short & nFFD_field, unsigned short** & degreeFFD_field) : nFFD(nFFD_field), DegreeFFD(degreeFFD_field) {
    this->name = option_field_name;
  }

  ~COptionFFDDegree() {};

  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    if ((option_value.size() == 1) && (option_value[0].compare("NONE") == 0)) {
      this->nFFD = 0;
      return "";
    }

    // Cannot have ; at the beginning or the end
    if (option_value[0].compare(";") == 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": may not have beginning semicolon");
      return newstring;
    }
    if (option_value[option_value.size()-1].compare(";") == 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": may not have ending semicolon");
      return newstring;
    }


    // use the ";" token to determine the number of design variables
    // This works because semicolon is not one of the delimiters in tokenize string
    this->nFFD = 0;
    for (unsigned int i = 0; i < static_cast<unsigned int>(option_value.size()); i++) {
      if (option_value[i].compare(";") == 0) {
        this->nFFD++;
      }
    }

    // One more design variable than semicolon
    this->nFFD++;

    this->DegreeFFD = new unsigned short*[this->nFFD];
    for (unsigned short iFFD = 0; iFFD < this->nFFD; iFFD++) {
      this->DegreeFFD[iFFD] = new unsigned short[3];
    }

    unsigned short nDegreeFFD = 0;
    stringstream ss;
    unsigned int i = 0;

    for (unsigned short iFFD = 0; iFFD < this->nFFD; iFFD++) {

      nDegreeFFD = 3;

      for (unsigned short iDegreeFFD = 0; iDegreeFFD < nDegreeFFD; iDegreeFFD++) {
        ss << option_value[i] << " ";
        ss >> this->DegreeFFD[iFFD][iDegreeFFD];
        i++;
      }

      if (iFFD < (this->nFFD-1)) {
        if (option_value[i].compare(";") != 0) {
          string newstring;
          newstring.append(this->name);
          newstring.append(": a FFD degree in the configuration file has the wrong number of parameters");
          return newstring;
        }
        i++;
      }

    }

    // Need to return something...
    return "";
  }

  void SetDefault() {
    this->nFFD = 0;
    this->DegreeFFD = NULL;
  }

};

// Class where the option is represented by (String, su2double, string, su2double, ...)
class COptionStringDoubleList : public COptionBase {
  string name; // identifier for the option
  unsigned short & size; // how many strings are there (same as number of su2doubles)

  string * & s_f; // Reference to the string fields
  su2double* & d_f; // reference to the su2double fields

public:
  COptionStringDoubleList(string option_field_name, unsigned short & list_size, string * & string_field, su2double* & double_field) : size(list_size), s_f(string_field), d_f(double_field) {
    this->name = option_field_name;
  }

  ~COptionStringDoubleList() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    // There must be an even number of entries (same number of strings and doubles
    unsigned short totalVals = option_value.size();
    if ((totalVals % 2) != 0) {
      if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)) {
        // It's okay to say its NONE
        this->size = 0;
        return "";
      }
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have an even number of entries");
      return newstring;
    }
    unsigned short nVals = totalVals / 2;
    this->size = nVals;
    this->s_f = new string[nVals];
    this->d_f = new su2double[nVals];

    for (unsigned long i = 0; i < nVals; i++) {
      this->s_f[i].assign(option_value[2*i]); // 2 because have su2double and string
      istringstream is(option_value[2*i + 1]);
      su2double val;
      if (!(is >> val)) {
        return badValue(option_value, "string su2double", this->name);
      }
      this->d_f[i] = val;
    }
    // Need to return something...
    return "";
  }

  void SetDefault() {
    this->size = 0; // There is no default value for list
  }
};

class COptionInlet : public COptionBase {
  string name; // identifier for the option
  unsigned short & size;
  string * & marker;
  su2double * & ttotal;
  su2double * & ptotal;
  su2double ** & flowdir;

public:
  COptionInlet(string option_field_name, unsigned short & nMarker_Inlet, string* & Marker_Inlet, su2double* & Ttotal, su2double* & Ptotal, su2double** & FlowDir) : size(nMarker_Inlet), marker(Marker_Inlet), ttotal(Ttotal), ptotal(Ptotal), flowdir(FlowDir) {
    this->name = option_field_name;
  }

  ~COptionInlet() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    unsigned short totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)) {
      this->size = 0;
      this->marker = NULL;
      this->ttotal = NULL;
      this->ptotal = NULL;
      this->flowdir = NULL;
      return "";
    }

    if (totalVals % 6 != 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have a number of entries divisible by 6");
      this->size = 0;
      this->marker = NULL;
      this->ttotal = NULL;
      this->ptotal = NULL;
      this->flowdir = NULL;
      return newstring;
    }

    unsigned short nVals = totalVals / 6;
    this->size = nVals;
    this->marker = new string[nVals];
    this->ttotal = new su2double[nVals];
    this->ptotal = new su2double[nVals];
    this->flowdir = new su2double*[nVals];
    for (unsigned long i = 0; i < nVals; i++) {
      this->flowdir[i] = new su2double[3];
    }

    for (unsigned long i = 0; i < nVals; i++) {
      this->marker[i].assign(option_value[6*i]);
      istringstream ss_1st(option_value[6*i + 1]);
      if (!(ss_1st >> this->ttotal[i])) {
        return badValue(option_value, "inlet", this->name);
      }
      istringstream ss_2nd(option_value[6*i + 2]);
      if (!(ss_2nd >> this->ptotal[i])) {
        return badValue(option_value, "inlet", this->name);
      }
      istringstream ss_3rd(option_value[6*i + 3]);
      if (!(ss_3rd >> this->flowdir[i][0])) {
        return badValue(option_value, "inlet", this->name);
      }
      istringstream ss_4th(option_value[6*i + 4]);
      if (!(ss_4th >> this->flowdir[i][1])) {
        return badValue(option_value, "inlet", this->name);
      }
      istringstream ss_5th(option_value[6*i + 5]);
      if (!(ss_5th >> this->flowdir[i][2])) {
        return badValue(option_value, "inlet", this->name);
      }
    }

    return "";
  }

  void SetDefault() {
    this->marker = NULL;
    this->ttotal = NULL;
    this->ptotal = NULL;
    this->flowdir = NULL;
    this->size = 0; // There is no default value for list
  }
};

template <class Tenum>
class COptionRiemann : public COptionBase {

protected:
  map<string, Tenum> m;
  string name; // identifier for the option
  unsigned short & size;
  string * & marker;
  unsigned short* & field; // Reference to the field name
  su2double * & var1;
  su2double * & var2;
  su2double ** & flowdir;

public:
  COptionRiemann(string option_field_name, unsigned short & nMarker_Riemann, string* & Marker_Riemann, unsigned short* & option_field, const map<string, Tenum> m, su2double* & var1, su2double* & var2, su2double** & FlowDir) : size(nMarker_Riemann),
                  marker(Marker_Riemann), field(option_field), var1(var1), var2(var2), flowdir(FlowDir) {
    this->name = option_field_name;
    this->m = m;
  }
  ~COptionRiemann() {};

  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    unsigned short totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)) {
      this->size = 0;
      this->marker = NULL;
      this->field = 0;
      this->var1 = NULL;
      this->var2 = NULL;
      this->flowdir = NULL;
      return "";
    }

    if (totalVals % 7 != 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have a number of entries divisible by 7");
      this->size = 0;
      this->marker = NULL;
      this->var1 = NULL;
      this->var2 = NULL;
      this->flowdir = NULL;
      this->field = NULL;
      return newstring;
    }

    unsigned short nVals = totalVals / 7;
    this->size = nVals;
    this->marker = new string[nVals];
    this->var1 = new su2double[nVals];
    this->var2 = new su2double[nVals];
    this->flowdir = new su2double*[nVals];
    this->field = new unsigned short[nVals];

    for (unsigned long i = 0; i < nVals; i++) {
      this->flowdir[i] = new su2double[3];
    }

    for (unsigned long i = 0; i < nVals; i++) {
      this->marker[i].assign(option_value[7*i]);
        // Check to see if the enum value is in the map
    if (this->m.find(option_value[7*i + 1]) == m.end()) {
      string str;
      str.append(this->name);
      str.append(": invalid option value ");
      str.append(option_value[0]);
      str.append(". Check current SU2 options in config_template.cfg.");
      return str;
    }
      Tenum val = this->m[option_value[7*i + 1]];
      this->field[i] = val;

      istringstream ss_1st(option_value[7*i + 2]);
      if (!(ss_1st >> this->var1[i])) {
        return badValue(option_value, "Riemann", this->name);
      }
      istringstream ss_2nd(option_value[7*i + 3]);
      if (!(ss_2nd >> this->var2[i])) {
        return badValue(option_value, "Riemann", this->name);
      }
      istringstream ss_3rd(option_value[7*i + 4]);
      if (!(ss_3rd >> this->flowdir[i][0])) {
        return badValue(option_value, "Riemann", this->name);
      }
      istringstream ss_4th(option_value[7*i + 5]);
      if (!(ss_4th >> this->flowdir[i][1])) {
        return badValue(option_value, "Riemann", this->name);
      }
      istringstream ss_5th(option_value[7*i + 6]);
      if (!(ss_5th >> this->flowdir[i][2])) {
        return badValue(option_value, "Riemann", this->name);
      }
    }

    return "";
  }

  void SetDefault() {
    this->marker = NULL;
    this->var1 = NULL;
    this->var2 = NULL;
    this->flowdir = NULL;
    this->size = 0; // There is no default value for list
  }
};

template <class Tenum>
class COptionGiles : public COptionBase{

  map<string, Tenum> m;
  unsigned short & size;
  string * & marker;
  unsigned short* & field; // Reference to the fieldname
  string name; // identifier for the option
  su2double * & var1;
  su2double * & var2;
  su2double ** & flowdir;
  su2double * & relfac1;
  su2double * & relfac2;

public:
  COptionGiles(string option_field_name, unsigned short & nMarker_Giles, string* & Marker_Giles, unsigned short* & option_field, const map<string, Tenum> m, su2double* & var1, su2double* & var2, su2double** & FlowDir, su2double* & relfac1, su2double* & relfac2) : size(nMarker_Giles),
               marker(Marker_Giles), field(option_field), var1(var1), var2(var2), flowdir(FlowDir), relfac1(relfac1), relfac2(relfac2) {
    this->name = option_field_name;
    this->m = m;
  }
  ~COptionGiles() {};

  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    unsigned long totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)) {
      this->size = 0;
      this->marker = NULL;
      this->field = 0;
      this->var1 = NULL;
      this->var2 = NULL;
      this->flowdir = NULL;
      this->relfac1 = NULL;
      this->relfac2 = NULL;
      return "";
    }

    if (totalVals % 9 != 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have a number of entries divisible by 9");
      this->size = 0;
      this->marker = NULL;
      this->var1 = NULL;
      this->var2 = NULL;
      this->flowdir = NULL;
      this->field = NULL;
      this->relfac1 = NULL;
      this->relfac2 = NULL;
      return newstring;
    }

    unsigned long nVals = totalVals / 9;
    this->size = nVals;
    this->marker = new string[nVals];
    this->var1 = new su2double[nVals];
    this->var2 = new su2double[nVals];
    this->flowdir = new su2double*[nVals];
    this->field = new unsigned short[nVals];
    this->relfac1 = new su2double[nVals];
    this->relfac2 = new su2double[nVals];

    for (unsigned int i = 0; i < nVals; i++) {
      this->flowdir[i] = new su2double[3];
    }

    for (unsigned int i = 0; i < nVals; i++) {
      this->marker[i].assign(option_value[9*i]);
        // Check to see if the enum value is in the map
    if (this->m.find(option_value[9*i + 1]) == m.end()) {
      string str;
      str.append(this->name);
      str.append(": invalid option value ");
      str.append(option_value[0]);
      str.append(". Check current SU2 options in config_template.cfg.");
      return str;
    }
      Tenum val = this->m[option_value[9*i + 1]];
      this->field[i] = val;

      istringstream ss_1st(option_value[9*i + 2]);
      if (!(ss_1st >> this->var1[i])) {
        return badValue(option_value, "Giles BC", this->name);
      }
      istringstream ss_2nd(option_value[9*i + 3]);
      if (!(ss_2nd >> this->var2[i])) {
        return badValue(option_value, "Giles BC", this->name);
      }
      istringstream ss_3rd(option_value[9*i + 4]);
      if (!(ss_3rd >> this->flowdir[i][0])) {
        return badValue(option_value, "Giles BC", this->name);
      }
      istringstream ss_4th(option_value[9*i + 5]);
      if (!(ss_4th >> this->flowdir[i][1])) {
        return badValue(option_value, "Giles BC", this->name);
      }
      istringstream ss_5th(option_value[9*i + 6]);
      if (!(ss_5th >> this->flowdir[i][2])) {
        return badValue(option_value, "Giles BC", this->name);
      }
      istringstream ss_6th(option_value[9*i + 7]);
      if (!(ss_6th >> this->relfac1[i])) {
        return badValue(option_value, "Giles BC", this->name);
      }
      istringstream ss_7th(option_value[9*i + 8]);
      if (!(ss_7th >> this->relfac2[i])) {
        return badValue(option_value, "Giles BC", this->name);
      }
    }

    return "";
  }

  void SetDefault() {
    this->marker = NULL;
    this->var1 = NULL;
    this->var2 = NULL;
    this->relfac1 = NULL;
    this->relfac2 = NULL;
    this->flowdir = NULL;
    this->size = 0; // There is no default value for list
  }
};

//Inlet condition where the input direction is assumed
class COptionExhaust : public COptionBase {
  string name; // identifier for the option
  unsigned short & size;
  string * & marker;
  su2double * & ttotal;
  su2double * & ptotal;

public:
  COptionExhaust(string option_field_name, unsigned short & nMarker_Exhaust, string* & Marker_Exhaust, su2double* & Ttotal, su2double* & Ptotal) : size(nMarker_Exhaust), marker(Marker_Exhaust), ttotal(Ttotal), ptotal(Ptotal) {
    this->name = option_field_name;
  }

  ~COptionExhaust() {};

  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    unsigned short totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)) {
      this->size = 0;
      this->marker = NULL;
      this->ttotal = NULL;
      this->ptotal = NULL;
      return "";
    }

    if (totalVals % 3 != 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have a number of entries divisible by 3");
      this->size = 0;
      this->marker = NULL;
      this->ttotal = NULL;
      this->ptotal = NULL;
      return newstring;
    }

    unsigned short nVals = totalVals / 3;
    this->size = nVals;
    this->marker = new string[nVals];
    this->ttotal = new su2double[nVals];
    this->ptotal = new su2double[nVals];

    for (unsigned long i = 0; i < nVals; i++) {
      this->marker[i].assign(option_value[3*i]);
      istringstream ss_1st(option_value[3*i + 1]);
      if (!(ss_1st >> this->ttotal[i]))
        return badValue(option_value, "exhaust fixed", this->name);
      istringstream ss_2nd(option_value[3*i + 2]);
      if (!(ss_2nd >> this->ptotal[i]))
        return badValue(option_value, "exhaust fixed", this->name);
    }

    return "";
  }

  void SetDefault() {
    this->marker = NULL;
    this->ttotal = NULL;
    this->ptotal = NULL;
    this->size = 0; // There is no default value for list
  }

};

class COptionPeriodic : public COptionBase {
  string name; // identifier for the option
  unsigned short & size;
  string * & marker_bound;
  string * & marker_donor;
  su2double ** & rot_center;
  su2double ** & rot_angles;
  su2double ** & translation;

public:
  COptionPeriodic(const string option_field_name, unsigned short & nMarker_PerBound,
                  string* & Marker_PerBound, string* & Marker_PerDonor,
                  su2double** & RotCenter, su2double** & RotAngles, su2double** & Translation) : size(nMarker_PerBound), marker_bound(Marker_PerBound), marker_donor(Marker_PerDonor), rot_center(RotCenter), rot_angles(RotAngles), translation(Translation) {
    this->name = option_field_name;
  }

  ~COptionPeriodic() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    const int mod_num = 11;

    unsigned short totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)) {
      this->size = 0;
      this->marker_bound = NULL;
      this->marker_donor = NULL;
      this->rot_center = NULL;
      this->rot_angles = NULL;
      this->translation = NULL;
      return "";
    }

    if (totalVals % mod_num != 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have a number of entries divisible by 11");
      this->size = 0;
      this->marker_bound = NULL;
      this->marker_donor = NULL;
      this->rot_center = NULL;
      this->rot_angles = NULL;
      this->translation = NULL;
      return newstring;
    }

    unsigned short nVals = 2 * (totalVals / mod_num); // To account for periodic and donor
    this->size = nVals;
    this->marker_bound = new string[nVals];
    this->marker_donor = new string[nVals];
    this->rot_center = new su2double*[nVals];
    this->rot_angles = new su2double*[nVals];
    this->translation = new su2double*[nVals];
    for (unsigned long i = 0; i < nVals; i++) {
      this->rot_center[i] = new su2double[3];
      this->rot_angles[i] = new su2double[3];
      this->translation[i] = new su2double[3];
    }

    su2double deg2rad = PI_NUMBER/180.0;

    for (unsigned long i = 0; i < (nVals/2); i++) {
      this->marker_bound[i].assign(option_value[mod_num*i]);
      this->marker_donor[i].assign(option_value[mod_num*i+1]);
      istringstream ss_1st(option_value[mod_num*i + 2]);
      if (!(ss_1st >> this->rot_center[i][0])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_2nd(option_value[mod_num*i + 3]);
      if (!(ss_2nd >> this->rot_center[i][1])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_3rd(option_value[mod_num*i + 4]);
      if (!(ss_3rd >> this->rot_center[i][2])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_4th(option_value[mod_num*i + 5]);
      if (!(ss_4th >> this->rot_angles[i][0])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_5th(option_value[mod_num*i + 6]);
      if (!(ss_5th >> this->rot_angles[i][1])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_6th(option_value[mod_num*i + 7]);
      if (!(ss_6th >> this->rot_angles[i][2])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_7th(option_value[mod_num*i + 8]);
      if (!(ss_7th >> this->translation[i][0])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_8th(option_value[mod_num*i + 9]);
      if (!(ss_8th >> this->translation[i][1])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_9th(option_value[mod_num*i + 10]);
      if (!(ss_9th >> this->translation[i][2])) {
        return badValue(option_value, "periodic", this->name);
      }
      this->rot_angles[i][0] *= deg2rad;
      this->rot_angles[i][1] *= deg2rad;
      this->rot_angles[i][2] *= deg2rad;
    }

    for (unsigned long i = (nVals/2); i < nVals; i++) {
      this->marker_bound[i].assign(option_value[mod_num*(i-nVals/2)+1]);
      this->marker_donor[i].assign(option_value[mod_num*(i-nVals/2)]);
      istringstream ss_1st(option_value[mod_num*(i-nVals/2) + 2]);
      if (!(ss_1st >> this->rot_center[i][0])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_2nd(option_value[mod_num*(i-nVals/2) + 3]);
      if (!(ss_2nd >> this->rot_center[i][1])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_3rd(option_value[mod_num*(i-nVals/2) + 4]);
      if (!(ss_3rd >> this->rot_center[i][2])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_4th(option_value[mod_num*(i-nVals/2) + 5]);
      if (!(ss_4th >> this->rot_angles[i][0])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_5th(option_value[mod_num*(i-nVals/2) + 6]);
      if (!(ss_5th >> this->rot_angles[i][1])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_6th(option_value[mod_num*(i-nVals/2) + 7]);
      if (!(ss_6th >> this->rot_angles[i][2])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_7th(option_value[mod_num*(i-nVals/2) + 8]);
      if (!(ss_7th >> this->translation[i][0])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_8th(option_value[mod_num*(i-nVals/2) + 9]);
      if (!(ss_8th >> this->translation[i][1])) {
        return badValue(option_value, "periodic", this->name);
      }
      istringstream ss_9th(option_value[mod_num*(i-nVals/2) + 10]);
      if (!(ss_9th >> this->translation[i][2])) {
        return badValue(option_value, "periodic", this->name);
      }
      /*--- Mirror the rotational angles and translation vector (rotational
       center does not need to move) ---*/
      this->rot_center[i][0] *= 1.0;
      this->rot_center[i][1] *= 1.0;
      this->rot_center[i][2] *= 1.0;
      this->rot_angles[i][0] *= -deg2rad;
      this->rot_angles[i][1] *= -deg2rad;
      this->rot_angles[i][2] *= -deg2rad;
      this->translation[i][0] *= -1.0;
      this->translation[i][1] *= -1.0;
      this->translation[i][2] *= -1.0;
    }

    return "";
  }

  void SetDefault() {
    this->size = 0;
    this->marker_bound = NULL;
    this->marker_donor = NULL;
    this->rot_center = NULL;
    this->rot_angles = NULL;
    this->translation = NULL;
  }
};

class COptionTurboPerformance : public COptionBase {
  string name; // identifier for the option
  unsigned short & size;
  string * & marker_turboIn;
  string * & marker_turboOut;

public:
  COptionTurboPerformance(const string option_field_name, unsigned short & nMarker_TurboPerf,
                          string* & Marker_TurboBoundIn, string* & Marker_TurboBoundOut) : size(nMarker_TurboPerf), marker_turboIn(Marker_TurboBoundIn), marker_turboOut(Marker_TurboBoundOut){
    this->name = option_field_name;
  }

  ~COptionTurboPerformance() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    const int mod_num = 2;

    unsigned long totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)) {
      this->size = 0;
      this->marker_turboIn= NULL;
      this->marker_turboOut = NULL;
      return "";
    }

    if (totalVals % mod_num != 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have a number of entries divisible by 2");
      this->size = 0;
      this->marker_turboIn= NULL;
      this->marker_turboOut = NULL;;
      return newstring;
    }

    unsigned long nVals = totalVals / mod_num;
    this->size = nVals;
    this->marker_turboIn = new string[nVals];
    this->marker_turboOut = new string[nVals];
    for (unsigned long i = 0; i < nVals; i++) {
      this->marker_turboIn[i].assign(option_value[mod_num*i]);
      this->marker_turboOut[i].assign(option_value[mod_num*i+1]);
     }


    return "";
  }

  void SetDefault() {
    this->size = 0;
    this->marker_turboIn= NULL;
    this->marker_turboOut = NULL;
  }
};

class COptionPython : public COptionBase {
  string name;
public:
  COptionPython(const string name) {
    this->name = name;
  }
  ~COptionPython() {};
  // No checking happens with python options
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    return "";
  }
  // No defaults with python options
  void SetDefault() {
    return;
  };
};

class COptionActDisk : public COptionBase {
  string name; // identifier for the option
  unsigned short & inlet_size;
  unsigned short & outlet_size;
  string * & marker_inlet;
  string * & marker_outlet;
  su2double ** & press_jump;
  su2double ** & temp_jump;
  su2double ** & omega;

public:
  COptionActDisk(const string name,
                 unsigned short & nMarker_ActDiskInlet, unsigned short & nMarker_ActDiskOutlet, string * & Marker_ActDiskInlet, string * & Marker_ActDiskOutlet,
                 su2double ** & ActDisk_PressJump, su2double ** & ActDisk_TempJump, su2double ** & ActDisk_Omega) :
  inlet_size(nMarker_ActDiskInlet), outlet_size(nMarker_ActDiskOutlet), marker_inlet(Marker_ActDiskInlet), marker_outlet(Marker_ActDiskOutlet),
  press_jump(ActDisk_PressJump), temp_jump(ActDisk_TempJump), omega(ActDisk_Omega) {
    this->name = name;
  }

  ~COptionActDisk() {};
  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    const int mod_num = 8;
    unsigned short totalVals = option_value.size();
    if ((totalVals == 1) && (option_value[0].compare("NONE") == 0)) {
      this->SetDefault();
      return "";
    }

    if (totalVals % mod_num != 0) {
      string newstring;
      newstring.append(this->name);
      newstring.append(": must have a number of entries divisible by 8");
      this->SetDefault();
      return newstring;
    }

    unsigned short nVals = totalVals / mod_num;
    this->inlet_size = nVals;
    this->outlet_size = nVals;
    this->marker_inlet = new string[this->inlet_size];
    this->marker_outlet = new string[this->outlet_size];

    this->press_jump = new su2double*[this->inlet_size];
    this->temp_jump = new su2double*[this->inlet_size];
    this->omega = new su2double*[this->inlet_size];
    for (int i = 0; i < this->inlet_size; i++) {
      this->press_jump[i] = new su2double[2];
      this->temp_jump[i] = new su2double[2];
      this->omega[i] = new su2double[2];
    }

    string tname = "actuator disk";

    for (int i = 0; i < this->inlet_size; i++) {
      this->marker_inlet[i].assign(option_value[mod_num*i]);
      this->marker_outlet[i].assign(option_value[mod_num*i+1]);
      istringstream ss_1st(option_value[mod_num*i + 2]);
      if (!(ss_1st >> this->press_jump[i][0])) {
        return badValue(option_value, tname, this->name);
      }
      istringstream ss_2nd(option_value[mod_num*i + 3]);
      if (!(ss_2nd >> this->temp_jump[i][0])) {
        return badValue(option_value, tname, this->name);
      }
      istringstream ss_3rd(option_value[mod_num*i + 4]);
      if (!(ss_3rd >> this->omega[i][0])) {
        return badValue(option_value, tname, this->name);
      }
      istringstream ss_4th(option_value[mod_num*i + 5]);
      if (!(ss_4th >> this->press_jump[i][1])) {
        return badValue(option_value, tname, this->name);
      }
      istringstream ss_5th(option_value[mod_num*i + 6]);
      if (!(ss_5th >> this->temp_jump[i][1])) {
        return badValue(option_value, tname, this->name);
      }
      istringstream ss_6th(option_value[mod_num*i + 7]);
      if (!(ss_6th >> this->omega[i][1])) {
        return badValue(option_value, tname, this->name);
      }
    }
    return "";
  }
  void SetDefault() {
    this->inlet_size = 0;
    this->outlet_size = 0;
    this->marker_inlet = NULL;
    this->marker_outlet = NULL;
    this->press_jump = NULL;
    this->temp_jump = NULL;
    this->omega = NULL;
  }
};

class COptionWallFunction : public COptionBase {
  string name; // identifier for the option
  unsigned short &nMarkers;
  string* &markers;
  unsigned short*  &walltype;
  unsigned short** &intInfo;
  su2double**      &doubleInfo;

public:
  COptionWallFunction(const string name, unsigned short &nMarker_WF,
                      string* &Marker_WF, unsigned short* &type_WF,
                      unsigned short** &intInfo_WF, su2double** &doubleInfo_WF) :
  nMarkers(nMarker_WF), markers(Marker_WF), walltype(type_WF),
  intInfo(intInfo_WF), doubleInfo(doubleInfo_WF) {
    this->name = name;
  }

  ~COptionWallFunction(){}

  string SetValue(vector<string> option_value) {
    COptionBase::SetValue(option_value);
    /*--- First check if NONE is specified. ---*/
    unsigned short totalSize = option_value.size();
    if ((totalSize == 1) && (option_value[0].compare("NONE") == 0)) {
      this->SetDefault();
      return "";
    }

    /*--- Determine the number of markers, for which a wall
          function treatment has been specified. ---*/
    unsigned short counter = 0, nVals = 0;
    while (counter < totalSize ) {

      /* Update the counter for the number of markers specified
         and store the current index for possible error messages. */
      ++nVals;
      const unsigned short indMarker = counter;

      /* Check if a wall function type has been specified for this marker.
         If not, create an error message and return. */
      ++counter;
      const unsigned short indWallType = counter;
      unsigned short typeWF = NO_WALL_FUNCTION;
      bool validWF = true;
      if (counter == totalSize) validWF = false;
      else {
        map<string, ENUM_WALL_FUNCTIONS>::const_iterator it;
        it = Wall_Functions_Map.find(option_value[counter]);
        if(it == Wall_Functions_Map.end()) validWF = false;
        else                               typeWF  = it->second;
      }

      if (!validWF ) {
        string newstring;
        newstring.append(this->name);
        newstring.append(": Invalid wall function type, ");
        newstring.append(option_value[counter]);
        newstring.append(", encountered for marker ");
        newstring.append(option_value[indMarker]);
        return newstring;
      }

      /* Update the counter, as the wall function type is valid. */
      ++counter;

      /*--- For some wall function types some additional info
            must be specified. Hence the counter must be updated
            accordingly. ---*/
      switch( typeWF ) {
        case EQUILIBRIUM_WALL_MODEL: counter += 3; break;
        case LOGARITHMIC_WALL_MODEL: counter += 3; break;
        default: break;
      }

      /* In case the counter is larger than totalSize, the data for
         this wall function type has not been specified correctly. */
      if (counter > totalSize) {
        string newstring;
        newstring.append(this->name);
        newstring.append(", marker ");
        newstring.append(option_value[indMarker]);
        newstring.append(", wall function type ");
        newstring.append(option_value[indWallType]);
        newstring.append(": Additional information is missing.");
        return newstring;
      }
    }

    /* Allocate the memory to store the data for the wall function markers. */
    this->nMarkers   = nVals;
    this->markers    = new string[nVals];
    this->walltype   = new unsigned short[nVals];
    this->intInfo    = new unsigned short*[nVals];
    this->doubleInfo = new su2double*[nVals];

    for (unsigned short i=0; i<nVals; i++) {
      this->intInfo[i]    = NULL;
      this->doubleInfo[i] = NULL;
    }

    /*--- Loop over the wall markers and store the info in the
          appropriate arrays. ---*/
    counter = 0;
    for (unsigned short i=0; i<nVals; i++) {

      /* Set the name of the wall function marker. */
      this->markers[i].assign(option_value[counter++]);

      /* Determine the wall function type. As their validaties have
         already been tested, there is no need to do so again. */
      map<string, ENUM_WALL_FUNCTIONS>::const_iterator it;
      it = Wall_Functions_Map.find(option_value[counter++]);

      this->walltype[i] = it->second;

      /*--- For some wall function types, some additional info
            is needed, which is extracted from option_value. ---*/
      switch( this->walltype[i] ) {

        case EQUILIBRIUM_WALL_MODEL: {

          /* LES equilibrium wall model. The exchange distance, stretching
             factor and number of points in the wall model must be specified. */
          this->intInfo[i]    = new unsigned short[1];
          this->doubleInfo[i] = new su2double[2];

          istringstream ss_1st(option_value[counter++]);
          if (!(ss_1st >> this->doubleInfo[i][0])) {
            return badValue(option_value, "su2double", this->name);
          }

          istringstream ss_2nd(option_value[counter++]);
          if (!(ss_2nd >> this->doubleInfo[i][1])) {
            return badValue(option_value, "su2double", this->name);
          }

          istringstream ss_3rd(option_value[counter++]);
          if (!(ss_3rd >> this->intInfo[i][0])) {
            return badValue(option_value, "unsigned short", this->name);
          }

          break;
        }

        case LOGARITHMIC_WALL_MODEL: {

          /* LES Logarithmic law-of-the-wall model. The exchange distance, stretching
           factor and number of points in the wall model must be specified. */
          this->intInfo[i]    = new unsigned short[1];
          this->doubleInfo[i] = new su2double[2];

          istringstream ss_1st(option_value[counter++]);
          if (!(ss_1st >> this->doubleInfo[i][0])) {
            return badValue(option_value, "su2double", this->name);
          }

          istringstream ss_2nd(option_value[counter++]);
          if (!(ss_2nd >> this->doubleInfo[i][1])) {
            return badValue(option_value, "su2double", this->name);
          }

          istringstream ss_3rd(option_value[counter++]);
          if (!(ss_3rd >> this->intInfo[i][0])) {
            return badValue(option_value, "unsigned short", this->name);
          }

          break;
        }

        default: // Just to avoid a compiler warning.
          break;
      }
    }

    // Need to return something...
    return "";
  }

  void SetDefault() {
    this->nMarkers   = 0;
    this->markers    = NULL;
    this->walltype   = NULL;
    this->intInfo    = NULL;
    this->doubleInfo = NULL;
  }
};
