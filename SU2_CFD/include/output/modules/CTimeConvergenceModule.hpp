#ifndef CTIMECONVERGENCEMODULE_HPP
#define CTIMECONVERGENCEMODULE_HPP

#include "COutputModule.hpp"
#include "../tools/CWindowingTools.hpp"

class CTimeConvergenceModule : public CSolverOutputModule {

  std::vector<std::string> wndConvFields;                /*!< \brief Name of the field to be monitored for convergence. */
  std::vector<std::vector<su2double> > WndCauchy_Serie;  /*!< \brief Complete Cauchy serial. */
  unsigned long nWndCauchy_Elems = 0;              /*!< \brief Total number of cauchy elems to monitor */
  su2double wndCauchyEps = 0;                      /*!< \brief Defines the threshold when to stop the solver. */

  std::vector<su2double> WndOld_Func;  /*!< \brief Old value of the objective function (the function which is monitored). */
  std::vector<su2double> WndNew_Func;  /*!< \brief Current value of the objective function (the function which is monitored). */
  su2double WndCauchy_Func = 0;       /*!< \brief Current value of the convergence indicator at one iteration. */
  su2double WndCauchy_Value = 0;      /*!< \brief Summed value of the convergence indicator. */
  bool TimeConvergence = false;   /*!< \brief To indicate, if the windowed time average of the time loop has converged*/

  map<std::string, CWindowedAverage> windowedTimeAverages;


public:

  CTimeConvergenceModule(CConfig* config, int nDim);

  void DefineHistoryFields(CHistoryOutFieldManager& historyFields) override;

  void DefineHistoryFieldModifier(CHistoryOutFieldManager& historyFields) override;

  void LoadHistoryData(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                       const IterationInfo& iterationInfo) override;

};
#endif // CTIMECONVERGENCEMODULE_HPP
