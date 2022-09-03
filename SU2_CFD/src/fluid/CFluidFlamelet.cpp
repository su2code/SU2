
#include "../include/fluid/CFluidFlamelet.hpp"
#include "../../../Common/include/containers/CLookUpTable.hpp"

CFluidFlamelet::CFluidFlamelet(CConfig *config, su2double value_pressure_operating) : CFluidModel() {

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    cout << "n_scalars="<<n_scalars<<endl;
    // this is already done?
    n_scalars = 2;
    config->SetNScalars(n_scalars);

    cout << "*************************************************" << endl;
    cout << "***** initializing the lookup table         *****" << endl;
    cout << "*************************************************" << endl;

    table_scalar_names.resize(n_scalars);
    table_scalar_names.at(I_ENTH) = "EnthalpyTot";
    table_scalar_names.at(I_PROGVAR) = "Progress Variable";

    config->SetLUTScalarNames(table_scalar_names);

    /*--- we currently only need one source term from the LUT ---*/
    n_table_sources = 1;
    config->SetNLUTSources(n_table_sources);

    table_source_names.resize(n_table_sources);
    table_source_names.at(I_SRC_TOT_PROGVAR) = "ProdRateTot-PV";

    config->SetLUTSourceNames(table_source_names);

    look_up_table = new CLookUpTable(config->GetFileNameLUT(),table_scalar_names.at(I_PROGVAR),table_scalar_names.at(I_ENTH));

    n_lookups = config->GetNLookups();
    table_lookup_names.resize(n_lookups);
    for (int i_lookup=0; i_lookup < n_lookups; ++i_lookup) {
      table_lookup_names.at(i_lookup) = config->GetLUTLookupName(i_lookup);
    }

    source_scalar.resize(n_scalars);
    lookup_scalar.resize(n_lookups);

    cout << "nr of scalar sources = "<< n_scalars << endl;
    cout << "nr of lookup variables = "<< n_lookups << endl;

    Pressure = value_pressure_operating;
}

CFluidFlamelet::~CFluidFlamelet() {

  if (look_up_table!=NULL) delete   look_up_table;
}

unsigned long CFluidFlamelet::SetScalarLookups(su2double *val_scalars){

  su2double enth   = val_scalars[I_ENTH];
  su2double prog   = val_scalars[I_PROGVAR];

  string name_enth = table_scalar_names.at(I_ENTH);
  string name_prog = table_scalar_names.at(I_PROGVAR);

  /* perform table look ups */
  unsigned long exit_code = look_up_table->LookUp_ProgEnth(table_lookup_names, lookup_scalar, prog, enth, name_prog, name_enth);

  return exit_code;
}

unsigned long CFluidFlamelet::SetScalarSources(su2double *val_scalars){

  su2double*         table_sources = new su2double[n_table_sources];
  vector<string>     look_up_tags;
  vector<su2double*> look_up_data;

  table_sources[0] = 0.0;



  su2double enth   = val_scalars[I_ENTH];
  su2double prog   = val_scalars[I_PROGVAR];

  string name_enth = table_scalar_names.at(I_ENTH);
  string name_prog = table_scalar_names.at(I_PROGVAR);

  for (int i_source=0; i_source < n_table_sources; ++i_source) {
    look_up_tags.push_back(table_source_names.at(i_source));
    look_up_data.push_back(&table_sources[i_source]);
  }
 
  /* perform table look ups */
  unsigned long exit_code = look_up_table->LookUp_ProgEnth(look_up_tags, look_up_data, prog, enth, name_prog, name_enth);

  source_scalar.at(I_ENTH) = 0;
  source_scalar.at(I_PROGVAR) = table_sources[I_SRC_TOT_PROGVAR];

  /*--- we clip at a small positive value --- */
  if (source_scalar.at(I_PROGVAR)<EPS){
  /* --- clip negative values of progress variable source term (this should not happen for a good lookup table) ---*/
    source_scalar.at(I_PROGVAR) = 0.0;
  }

  delete [] table_sources;

  /*--- not used at the moment ---*/ 
  return exit_code;
}

unsigned long CFluidFlamelet::SetTDState_T(su2double val_temperature, const su2double* val_scalars){

  su2double val_enth = val_scalars[I_ENTH];
  su2double val_prog = val_scalars[I_PROGVAR];

  string name_enth = table_scalar_names.at(I_ENTH);
  string name_prog = table_scalar_names.at(I_PROGVAR);

  vector<string>     look_up_tags;
  vector<su2double*> look_up_data;

  unsigned long exit_code;

  /* add all quantities and their address to the look up vectors */
  // nijso TODO: check if these exist in the lookup table

  
  look_up_tags.push_back("Temperature");
  look_up_data.push_back(&Temperature);
  look_up_tags.push_back("Density");
  look_up_data.push_back(&Density);
  look_up_tags.push_back("Cp");
  look_up_data.push_back(&Cp);
  look_up_tags.push_back("ViscosityDyn");
  look_up_data.push_back(&Mu);
  look_up_tags.push_back("Conductivity");
  look_up_data.push_back(&Kt);
  look_up_tags.push_back("Diffusivity");
  look_up_data.push_back(&mass_diffusivity);

  //look_up_tags.push_back("HeatRelease");
  //look_up_data.push_back(&source_energy);
  //Temperature = 500.0;
  //Density = 1.2;
  //Cp = 1000.0;
  //Mu = 1.5e-5;
  //Kt = 0.025;
  //mass_diffusivity = 1e-6;

  /* perform table look ups */
  exit_code = look_up_table->LookUp_ProgEnth(look_up_tags,look_up_data, val_prog, val_enth,name_prog,name_enth);

  // nijso: is Cv used somewhere?
  // according to cristopher, yes!
  // we could check for the existence of molar_weight_mix in the lookup table, and else we just use gamma
  // default value is 1.4
  Cv = Cp/1.4;
  //Cv = Cp - UNIVERSAL_GAS_CONSTANT / (molar_weight_mix / 1000.);
  return exit_code;
}

unsigned long CFluidFlamelet::GetEnthFromTemp(su2double *val_enth, su2double val_prog, su2double val_temp){

  su2double          delta_temp_final = 0.01 ; /* convergence criterion for temperature in [K] */
  su2double          enth_iter        = 0. ;   /* in CH4/Air flames, 0 is usually a good initial value for the iteration */
  su2double          delta_enth;
  su2double          delta_temp_iter  = 1e10;
  unsigned long      exit_code        = 0;
  vector<string>     look_up_tags;
  vector<su2double*> look_up_data;
  int                counter_limit    = 50 ;
  string name_prog = table_scalar_names.at(I_PROGVAR);
  string name_enth = table_scalar_names.at(I_ENTH);

  /* set up look up vectors */
  su2double  temp_iter;
  look_up_tags.push_back("Temperature");
  look_up_data.push_back(&temp_iter);

  su2double  cp_iter;
  look_up_tags.push_back("Cp");
  look_up_data.push_back(&cp_iter);

  
  int counter = 0;
  while ( (abs(delta_temp_iter) > delta_temp_final) && (counter++ < counter_limit) ){

    /* look up temperature and heat capacity */
    look_up_table->LookUp_ProgEnth(look_up_tags, look_up_data, val_prog, enth_iter, name_prog, name_enth);

    /* calculate delta_temperature */
    delta_temp_iter = val_temp - temp_iter;

    /* calculate delta_enthalpy following dh = cp * dT */
    delta_enth      = cp_iter * delta_temp_iter;

    /* update enthalpy */
    enth_iter      += delta_enth;

  }

  /* set enthalpy value */
  *val_enth = enth_iter;

  if (counter >= counter_limit) {
    exit_code = 1;

  }

  return exit_code;
}
