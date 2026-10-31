/* Native regression checks for output field selection, failures, limits and accounting. LGPL-2.1-or-later. */
#include "../../../SU2_CFD/include/output/filewriter/CParaviewFieldLayout.hpp"
#include "../../../SU2_CFD/include/output/filewriter/CParaviewFileWriter.hpp"
#include "../../../SU2_CFD/include/output/filewriter/CParaviewBinaryFileWriter.hpp"
#include "../../../SU2_CFD/include/output/filewriter/CParaviewXMLFileWriter.hpp"
#include "../../../SU2_CFD/include/output/filewriter/CParaviewPVTUFileWriter.hpp"
#include "../../../SU2_CFD/include/output/filewriter/CSU2BinaryFileWriter.hpp"
#include <array>
#include <chrono>
#include <cmath>
#include <climits>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <thread>
#ifndef _WIN32
#include <csignal>
#include <sys/resource.h>
#endif

static void require(bool value, const char* message) {
  if (!value) SU2_MPI::Error(message, CURRENT_FUNCTION);
}

class Fixture final : public CParallelDataSorter {
 public:
  Fixture(unsigned short dim, const std::vector<std::string>& selected, bool empty = false)
      : CParallelDataSorter(nullptr, Names(dim, selected)) {
    nDim = dim;
    nPoints = empty ? 0 : (dim == 3 ? 8 : 4);
    nPointsGlobal = nPoints * size;
    linearPartitioner.Initialize(nPointsGlobal, 0);
    dataBuffer = new passivedouble[nPoints * GlobalField_Counter];
    for (unsigned long point = 0; point < nPoints; ++point)
      for (unsigned short field = 0; field < GlobalField_Counter; ++field)
        dataBuffer[point * GlobalField_Counter + field] = rank * 1000 + point * 10 + field;
    const auto type = dim == 3 ? HEXAHEDRON : QUADRILATERAL;
    auto** conn = dim == 3 ? &Conn_Hexa_Par : &Conn_Quad_Par;
    *conn = new int[nPoints];
    for (unsigned long p = 0; p < nPoints; ++p) (*conn)[p] = rank * nPoints + p + 1;
    nElem = empty ? 0 : 1; nElemGlobal = nElem * size;
    nConn = nPoints; nConnGlobal = nPointsGlobal;
    nElemPerType[TypeMap.at(type)] = nElem;
    nElemPerTypeGlobal[TypeMap.at(type)] = nElemGlobal;
    for (int r = 0; r <= size; ++r) { nElem_Cum[r] = r * nElem; nElemConn_Cum[r] = r * nConn; }
    connectivitySorted = true;
  }
  static std::vector<std::string> Names(unsigned short dim, const std::vector<std::string>& selected) {
    std::vector<std::string> names = {"x", "y"};
    if (dim == 3) names.push_back("z");
    names.insert(names.end(), selected.begin(), selected.end());
    return names;
  }
  void OversizedRestart() { nPointsGlobal = static_cast<unsigned long>(INT_MAX) + 1; }
};

class AccountingWriter final : public CFileWriter {
 public:
  AccountingWriter() : CFileWriter(".bin") {}
  using CFileWriter::OpenMPIFile;
  using CFileWriter::CloseMPIFile;
  using CFileWriter::WriteMPIString;
  using CFileWriter::WriteMPIBinaryData;
  using CFileWriter::WriteMPIBinaryDataAll;
};

int main(int argc, char** argv) {
  SU2_MPI::Init(&argc, &argv);
  const int rank = SU2_MPI::GetRank(), size = SU2_MPI::GetSize();
  const std::string mode = argc > 1 ? argv[1] : "normal";
  if (mode == "limit-node" || mode == "limit-product") {
    if (mode == "limit-node") CParallelDataSorter::CheckIntCount(static_cast<unsigned long>(INT_MAX) + 1);
    else CParallelDataSorter::CheckIntCount(static_cast<unsigned long>(INT_MAX) / 8 + 1, 8);
    SU2_MPI::Finalize();
    return 0;
  }
  if (mode == "limit-restart") {
    Fixture data(3, {"Density"}); data.OversizedRestart();
    CSU2BinaryFileWriter writer(&data); writer.WriteData("oversized");
    SU2_MPI::Finalize();
    return 0;
  }
  if (mode == "piece-error") {
    if (!rank) {
      std::filesystem::create_directory("failed-pieces");
      std::filesystem::create_symlink("/dev/full", "failed-pieces/failed-pieces_" + std::to_string(size - 1) + ".vtu");
    }
    SU2_MPI::Barrier(SU2_MPI::GetComm());
    Fixture data(3, {"Density"});
    CParaviewPVTUFileWriter writer(&data, true); writer.WriteData("failed-pieces");
    SU2_MPI::Finalize();
    return 0;
  }
#ifndef _WIN32
  if (mode == "shared-error") {
    AccountingWriter writer;
    writer.OpenMPIFile("failed-shared");
    struct rlimit limit;
    require(getrlimit(RLIMIT_FSIZE, &limit) == 0, "getrlimit failed");
    limit.rlim_cur = 1;
    std::signal(SIGXFSZ, SIG_IGN);
    require(setrlimit(RLIMIT_FSIZE, &limit) == 0, "setrlimit failed");
    const std::array<char, 256> data{};
    writer.WriteMPIBinaryDataAll(data.data(), data.size(), data.size() * size, data.size() * rank);
    writer.CloseMPIFile();
    SU2_MPI::Finalize();
    return 0;
  }
#endif
  CParallelDataSorter::CheckIntCount(INT_MAX);
  CParallelDataSorter::CheckIntCount(INT_MAX / 8, 8);
  AccountingWriter accounting;
  accounting.OpenMPIFile("accounting");
  const int header = 123;
  accounting.WriteMPIString("H", 0);
  accounting.WriteMPIBinaryData(&header, sizeof(header), 0);
  const char bytes[] = {'a', 'b', 'c'};
  accounting.WriteMPIBinaryDataAll(bytes, 3, 3 * size, 3 * rank);
  if (rank == size - 1) std::this_thread::sleep_for(std::chrono::milliseconds(80));
  accounting.CloseMPIFile();
  require(SU2_TYPE::GetValue(accounting.GetFilesize()) == 5 + 3 * size, "Replicated header bytes were counted twice");
  require(SU2_TYPE::GetValue(accounting.GetUsedTime()) >= .075, "Writer timer omitted work before file close");
  const auto seconds = SU2_TYPE::GetValue(accounting.GetUsedTime());
  require(std::abs(SU2_TYPE::GetValue(accounting.GetBandwidth()) - (5 + 3 * size) / (1e6 * seconds)) < 1e-12,
          "Incorrect writer bandwidth");
  const std::vector<std::vector<std::string>> selections = {
      {"Momentum_x"}, {"Momentum_y"}, {"Momentum_z"}, {"Momentum_x", "Momentum_z"},
      {"Momentum_x", "Other_y", "Momentum_z"}, {"Flux_x_extra"},
      {"Momentum_x", "Momentum_y", "Momentum_z", "Energy"}};
  for (unsigned short dim : {2, 3}) for (size_t i = 0; i < selections.size(); ++i) {
    auto selected = selections[i];
    if (dim == 2 && i == selections.size() - 1) selected.erase(selected.begin() + 2);
    Fixture data(dim, selected);
    const auto prefix = "d" + std::to_string(dim) + "s" + std::to_string(i);
    auto stats = [&](const CFileWriter& writer, const std::string& name, const std::string& extension) {
      if (!rank) std::cout << "WRITER_BYTES " << name << extension << " " << SU2_TYPE::GetValue(writer.GetFilesize()) << "\n";
    };
    CParaviewXMLFileWriter xml(&data, true); xml.WriteData(prefix + "-xml");
    stats(xml, prefix + "-xml", ".vtu");
    CParaviewPVTUFileWriter pvtu(&data, true); pvtu.WriteData(prefix + "-pieces");
    stats(pvtu, prefix + "-pieces", ".pvtu");
    CParaviewFileWriter ascii(&data); ascii.WriteData(prefix + "-ascii");
    stats(ascii, prefix + "-ascii", ".vtk");
    CParaviewBinaryFileWriter binary(&data); binary.WriteData(prefix + "-binary");
    stats(binary, prefix + "-binary", ".vtk");
  }
  Fixture empty(3, {"Momentum_y"}, true);
  CParaviewPVTUFileWriter emptyWriter(&empty); emptyWriter.WriteData("empty");
  if (!rank) std::cout << "PASS: writer completion, accounting and field fixtures\n";
  SU2_MPI::Finalize();
}
