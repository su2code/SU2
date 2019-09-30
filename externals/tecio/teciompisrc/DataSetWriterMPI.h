#include "DataSetWriter.h"
#include "ThirdPartyHeadersBegin.h"
#include <vector>
#include <mpi.h>
#include "ThirdPartyHeadersEnd.h"
namespace tecplot { namespace ___3933 { class ___37; class FileWriterInterface; class ___1844; }} namespace tecplot { namespace teciompi { class DataSetWriterMPI : public ___3933::DataSetWriter { public: DataSetWriterMPI( ___3933::___37*        ___36, ___3501                    vars, ___3501                    ___4671, ___3933::___1844 const&         maxIJKSubzoneSize, ___2090::ItemOffset_t maxFESubzoneSize, MPI_Comm                  communicator, int                       mainProcess, int                       localProcess, bool                      flushToDisk = false); virtual ~DataSetWriterMPI(); virtual ___372 writeDataSet( tecplot::___3933::FileWriterInterface& szpltFile, tecplot::___3933::___1392&        szpltZoneHeaderFileLocs); private: MPI_Comm m_communicator; int m_mainProcess; int m_localProcess; }; }}
