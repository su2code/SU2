#ifndef CCONVERGENCEMODULE_HPP
#define CCONVERGENCEMODULE_HPP

#include "COutputModule.hpp"

class CConvergenceModule : public CModifierModule {

  su2double cauchyValue = 0,         /*!< \brief Summed value of the convergence indicator. */
  cauchyFunc = 0;                    /*!< \brief Current value of the convergence indicator at one iteration. */
  unsigned short Cauchy_Counter = 0; /*!< \brief Number of elements of the Cauchy serial. */
  std::vector<std::vector<su2double> > cauchySerie;        /*!< \brief Complete Cauchy serial. */
  unsigned long nCauchy_Elems = 0;   /*!< \brief Total number of cauchy elems to monitor */
  su2double cauchyEps = 0;           /*!< \brief Defines the threshold when to stop the solver. */
  su2double minLogResidual = 0;      /*!< \brief Minimum value of the residual to reach */
  std::vector<su2double> oldFunc,     /*!< \brief Old value of the coefficient. */
  newFunc;                       /*!< \brief Current value of the coefficient. */
  bool convergence = false;              /*!< \brief To indicate if the solver has converged or not. */
  su2double initResidual = 0;        /*!< \brief Initial value of the residual to evaluate the convergence level. */
  std::vector<std::string> convFields;     /*!< \brief Name of the field to be monitored for convergence. */
public:

  CConvergenceModule(CConfig* config, int nDim);

  void DefineHistoryFields(CHistoryOutFieldManager& historyFields) override;

  void DefineHistoryFieldModifier(CHistoryOutFieldManager& historyFields) override;

  void LoadHistoryDataModifier(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                               const IterationInfo& iterationInfo) override;

};

#endif // CCONVERGENCEMODULE_HPP
