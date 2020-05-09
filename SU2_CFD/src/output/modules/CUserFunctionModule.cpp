#include "../../../include/output/modules/CUserFunctionModule.hpp"
#include "../../../../Common/include/CConfig.hpp"


CUserFunctionModule::CUserFunctionModule(CConfig *config){

  historyScope = interpreter::globalScope.getChild();
  volumeScope  = interpreter::globalScope.getChild();

  historyUserFunctions = interpreter::GetUserFunctions(interpreter::globalScope, {interpreter::FunctionType::HISTFIELD});
  volumeUserFunctions  = interpreter::GetUserFunctions(interpreter::globalScope, {interpreter::FunctionType::VOLUMEFIELD,
                                                                                  interpreter::FunctionType::VOLUMEINTEGRAL,
                                                                                  interpreter::FunctionType::SURFACEINTEGRAL});

}

void CUserFunctionModule::DefineVolumeFields(COutFieldCollection &fieldCollection){

  funcVolFieldRefs.clear();

  for (const auto function : volumeUserFunctions){
    if (function->getType() == interpreter::FunctionType::SURFACEINTEGRAL){
      const std::string funcName = function->name();
      auto newField = fieldCollection.AddItem(funcName, COutputField(funcName, -1, "CUSTOM", "User-defined field", FieldType::SURFACE_INTEGRATE));
      funcVolFieldRefs.push_back({newField, function});
    }
  }

  volFieldTokenRefs.clear();

  for (const auto& field : fieldCollection.GetReferencesAll()){
    packToken* ref = &volumeScope[field->first];
    volFieldTokenRefs.push_back({field, ref});
  }
}

void CUserFunctionModule::LoadSurfaceData(COutFieldCollection &fieldCollection){
  EvalUserFunctions(volFieldTokenRefs, funcVolFieldRefs, volumeScope);
}
void CUserFunctionModule::DefineHistoryFields(COutFieldCollection &fieldCollection){

  funcHistFieldRefs.clear();

  for (const auto function : historyUserFunctions){
    const std::string funcName = function->name();
    auto newField = fieldCollection.AddItem(funcName, COutputField(funcName, ScreenOutputFormat::SCIENTIFIC, "CUSTOM", FieldType::DEFAULT, "User-defined field"));
    funcHistFieldRefs.push_back({newField, function});
  }
}

void CUserFunctionModule::DefineHistoryFieldModifier(COutFieldCollection &fieldCollection){

  histFieldTokenRefs.clear();

  for (const auto& field : fieldCollection.GetReferencesAll()){
    packToken* ref = &historyScope[field->first];
    histFieldTokenRefs.push_back({field,ref});
  }
}

void CUserFunctionModule::LoadHistoryDataModifier(COutFieldCollection &fieldCollection){
  EvalUserFunctions(histFieldTokenRefs, funcHistFieldRefs, historyScope);
}


void CUserFunctionModule::EvalUserFunctions(const std::vector<FieldRefTokenPair> &fieldTokenRef,
                                            const std::vector<FieldRefFuncPair> &funcFieldRef, TokenMap &scope){
  for (const auto& fieldToken : fieldTokenRef){
    decltype (fieldToken.second) token;
    decltype (fieldToken.first) field;
    std::tie(field, token) = fieldToken;
    if (token != nullptr){
      (*token) = field->second.value;
    }
  }
  for (const auto& funcField : funcFieldRef){
    decltype (funcField.second) func;
    decltype (funcField.first) field;
    std::tie(field, func) = funcField;
    try {
      field->second.value = func->exec(scope).asDouble();
    }  catch (msg_exception &err) {
      SU2_MPI::Error(std::string("In expression ") + field->first
                     + std::string(": ") + std::string(err.what()), CURRENT_FUNCTION);
    }
  }
}

std::string CUserFunctionModule::GetName(const std::string& baseName){
  std::string newName{baseName};
  replace_if(newName.begin(),newName.end(),::ispunct,'_');
  replace_if(newName.begin(),newName.end(),::isblank,'_');
  return newName;
}