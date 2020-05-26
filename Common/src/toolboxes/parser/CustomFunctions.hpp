#pragma once

struct SurfaceSum : public exprtk::igeneric_function<su2double>
{
  typedef typename exprtk::igeneric_function<su2double>::parameter_list_t
  parameter_list_t;
  typedef typename generic_type::scalar_view scalar_t;
  typedef typename generic_type::vector_view vector_t;
  typedef typename generic_type::string_view string_t;
  SymbolTable* table;
  std::string name;

  SurfaceSum(SymbolTable *table, const std::string& name) :
    table(table), name(name) {}

  inline su2double operator()(parameter_list_t parameters) {
    su2double val = 0.0;
    for (std::size_t i = 0; i < parameters.size(); ++i) {
      generic_type& gt = parameters[i];

      if (generic_type::e_string == gt.type) {
        string_t st(gt);
        auto varName = table->get_variable(name + exprtk::to_str(st));
        if (!varName){
          SU2_MPI::Error("Cannot find variable " + name + exprtk::to_str(st), CURRENT_FUNCTION);
        }
        val += varName->value();
      } else {
        SU2_MPI::Error("In " + name + ": argument is not a string", CURRENT_FUNCTION);
      }

    }

    return val;
  }
};

