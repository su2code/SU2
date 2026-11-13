/* Shared field selection for legacy VTK, VTU and PVTU output. LGPL-2.1-or-later. */
#pragma once
#include <algorithm>
#include <string>
#include <vector>

struct CParaviewField {
  std::size_t first;
  unsigned short components;
  std::string name;
};

inline std::vector<CParaviewField> GetParaviewFields(const std::vector<std::string>& fieldNames,
                                                  unsigned short nDim) {
  auto names = fieldNames;
  for (auto& name : names) name.erase(std::remove(name.begin(), name.end(), '"'), name.end());
  std::vector<CParaviewField> fields;
  for (std::size_t i = nDim; i < names.size();) {
    const auto& name = names[i];
    const bool startsVector = name.size() > 2 && name.compare(name.size() - 2, 2, "_x") == 0;
    const auto base = startsVector ? name.substr(0, name.size() - 2) : name;
    // A component selected on its own remains a scalar; never read missing or unrelated fields.
    const bool complete = (nDim == 2 || nDim == 3) && startsVector && i + nDim <= names.size() && names[i + 1] == base + "_y" &&
                          (nDim == 2 || names[i + 2] == base + "_z");
    fields.push_back({i, static_cast<unsigned short>(complete ? 3 : 1), complete ? base : name});
    i += complete ? nDim : 1;
  }
  return fields;
}
