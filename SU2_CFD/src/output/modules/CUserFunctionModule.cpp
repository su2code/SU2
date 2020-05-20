#include "../../../include/output/modules/CUserFunctionModule.hpp"
#include "../../../../Common/include/CConfig.hpp"

CUserFunctionModule::CUserFunctionModule(CConfig *config, int nDim) :
  CSolverOutputModule(nDim),
  historyScope(interpreter::globalScope.getChild()),
  volumeScope(interpreter::globalScope.getChild()),
  historyUserFunctions(interpreter::GetUserFunctions(interpreter::globalScope, {interpreter::FunctionType::HISTFIELD})),
  volumeUserFunctions(interpreter::GetUserFunctions(interpreter::globalScope, {interpreter::FunctionType::VOLUMEFIELD,
                                                                               interpreter::FunctionType::VOLUMEINTEGRAL,
                                                                               interpreter::FunctionType::SURFACEINTEGRAL}))
{

  auto addInlineUserFunction = [&](std::string field,
                                   std::vector<string>& userFuncsVector){
    if (*field.begin() == '{' && *(field.end()-1) == '}'){
      userFuncsVector.push_back(field);
    }
  };

  for (int iField = 0; iField < config->GetnHistoryOutput(); iField++){
    addInlineUserFunction(config->GetHistoryOutput_Field(iField), inlineHistoryUserFunctions);
  }
  for (int iField = 0; iField < config->GetnScreenOutput(); iField++){
    addInlineUserFunction(config->GetScreenOutput_Field(iField), inlineHistoryUserFunctions);
  }
  for (int iField = 0; iField < config->GetnVolumeOutput(); iField++){
    addInlineUserFunction(config->GetVolumeOutput_Field(iField), inlineVolumeUserFunctions);
  }
  for (int iMarker_CfgFile = 0; iMarker_CfgFile <  config->GetnMarker_CfgFile(); iMarker_CfgFile++){
    const string markerNameCfg = config->GetMarker_CfgFile_TagBound(iMarker_CfgFile);
    historyScope[markerNameCfg] = markerNameCfg;
  }

}

void CUserFunctionModule::DefineVolumeFields(CVolumeOutFieldManager &volumeFields){

  funcVolFieldRefs.clear();

  for (const auto function : volumeUserFunctions){
    if (function->getType() == interpreter::FunctionType::SURFACEINTEGRAL){
      const std::string funcName = function->name();
      auto newField = volumeFields.AddField(funcName, funcName, "CUSTOM", "User-defined field", FieldType::SURFACE_INTEGRATE);
      funcVolFieldRefs.emplace_back(newField, function);
    }
  }

  for (const auto& inlineFunction : inlineVolumeUserFunctions){
    auto newField = volumeFields.AddField(inlineFunction, inlineFunction, "CUSTOM", "User-defined field", FieldType::DEFAULT);
    funcVolFieldRefs.emplace_back(newField, CreateInlineUserFunction(inlineFunction, interpreter::FunctionType::VOLUMEFIELD));
  }

  volFieldTokenRefs.clear();

  for (const auto& field : volumeFields.GetCollection().GetReferencesAll()){

    /*--- Create a packToken that contains a double value and store in our scope using the field name ---*/

    volumeScope[field->first] = packToken(0.0);

    /*--- Since we explicitely created a packToken with double values, we know the underlying type.
     * To avoid expensive new/delete operations and to have direct access to the underlying double variable
     * we store a reference to the true type here ---*/

    auto *ref = static_cast<Token<su2double>*>(volumeScope[field->first].token());
    volFieldTokenRefs.emplace_back(field, ref);
  }
}

void CUserFunctionModule::LoadSurfaceData(CVolumeOutFieldManager &volumeFields){
  EvalUserFunctions(volFieldTokenRefs, funcVolFieldRefs, volumeScope, interpreter::FunctionType::SURFACEINTEGRAL);
}

void CUserFunctionModule::LoadVolumeData(CVolumeOutFieldManager &volumeFields){
  EvalUserFunctions(volFieldTokenRefs, funcVolFieldRefs, volumeScope, interpreter::FunctionType::VOLUMEFIELD);
}

void CUserFunctionModule::DefineHistoryFields(CHistoryOutFieldManager &historyFields){

  funcHistFieldRefs.clear();

  for (const auto function : historyUserFunctions){
    const std::string funcName = function->name();
    auto newField = historyFields.AddField(funcName, funcName, ScreenOutputFormat::SCIENTIFIC, "CUSTOM", "User-defined field", FieldType::DEFAULT);
    funcHistFieldRefs.emplace_back(newField, function);
  }

  for (const auto& inlineFunction : inlineHistoryUserFunctions){
    auto newField = historyFields.AddField(inlineFunction, inlineFunction, ScreenOutputFormat::SCIENTIFIC, "CUSTOM", "User-defined field", FieldType::DEFAULT);
    funcHistFieldRefs.emplace_back(newField, CreateInlineUserFunction(inlineFunction, interpreter::FunctionType::HISTFIELD));
  }
}

void CUserFunctionModule::DefineHistoryFieldModifier(CHistoryOutFieldManager &historyFields){

  histFieldTokenRefs.clear();

  for (const auto& field : historyFields.GetCollection().GetReferencesAll()){

    /*--- Create a packToken that contains a double value and store in our scope using the field name ---*/

    historyScope[field->first] = packToken(0.0);

    /*--- Since we explicitely created a packToken with double values, we know the underlying type.
     * To avoid expensive new/delete operations and to have direct access to the underlying double variable
     * we store a reference to the true type here ---*/

    auto* ref = static_cast<Token<su2double>*>(historyScope[field->first].token());
    histFieldTokenRefs.emplace_back(field, ref);
  }

  for (const auto field : historyFields.GetCollection().GetFieldsByType({FieldType::COEFFICIENT})){
    historyScope[field->first + "@"] = CppFunction(&CUserFunctionModule::getSurfaceValue, {}, field->first);
  }

}

void CUserFunctionModule::LoadHistoryDataModifier(CHistoryOutFieldManager &historyFields){
  EvalUserFunctions(histFieldTokenRefs, funcHistFieldRefs, historyScope, interpreter::FunctionType::HISTFIELD);
}


void CUserFunctionModule::EvalUserFunctions(const std::vector<FieldRefTokenPair> &fieldTokenRef,
                                            const std::vector<FieldRefFuncPair> &funcFieldRef, TokenMap &scope, interpreter::FunctionType type){
  for (const auto& fieldToken : fieldTokenRef){
    decltype (fieldToken.second) token;
    decltype (fieldToken.first) field;
    std::tie(field, token) = fieldToken;
    if (token != nullptr){
      token->val = field->second.value;
    }
  }
  for (const auto& funcField : funcFieldRef){
    decltype (funcField.second) func;
    decltype (funcField.first) field;
    std::tie(field, func) = funcField;
    if (func->getType() == type){
      try {
        field->second.value = func->exec(scope).asDouble();
      }  catch (msg_exception &err) {
        SU2_MPI::Error(std::string("In expression ") + field->first
                       + std::string(": ") + std::string(err.what()), CURRENT_FUNCTION);
      }
    }
  }
}

interpreter::UserFunction* CUserFunctionModule::CreateInlineUserFunction(string name, interpreter::FunctionType type){
  interpreter::BlockStatement body;
  string rawName = name;
  rawName.pop_back();
  rawName.erase(0,1);

  string funcName = rawName;
  replace_if(funcName.begin(),funcName.end(),::ispunct,'_');
  replace_if(funcName.begin(),funcName.end(),::isblank,'_');
  std::string func = "{"
                     " var = " + rawName + ";"
                     " return var; }";

  body.compile(func.c_str());

  return new interpreter::UserFunction({},body, funcName, type);

}

packToken CUserFunctionModule::getSurfaceValue(TokenMap scope, const std::string &name){

  vector<std::string> markerName;
  TokenList extraArgs = scope["args"].asList();

  if (extraArgs.list().empty()){
    throw(msg_exception("Expecting a marker name."));
  }

  for (const auto& item : extraArgs.list()){
    markerName.push_back(item.asString());
  }

  su2double val = 0.0;
  for (const auto& marker : markerName){
    packToken* token = scope.find(name + "@" + marker);
    if (token != nullptr){
      val += token->asDouble();
    } else {
      throw(msg_exception(string("Unkown marker ") + "\"" + marker + "\""));
    }
  }

  return val;

}
