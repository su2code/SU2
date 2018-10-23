#include "SzlFileLoader.h"
#include "fileStuff.h"
#include "ThirdPartyHeadersBegin.h"
#  include <string>
#  include <stdarg.h>
#include "ThirdPartyHeadersEnd.h"
using std::string; namespace tecplot { namespace ___3933 { std::string getFileNameSansFolder(std::string const& ___1394) { REQUIRE(!___1394.empty()); size_t beginBaseFileNamePos; size_t backslashPos = ___1394.rfind('\\'); size_t forwardSlashPos = ___1394.rfind('/'); if ( backslashPos != string::npos ) if ( forwardSlashPos != string::npos ) beginBaseFileNamePos = std::max(backslashPos,forwardSlashPos)+1; else beginBaseFileNamePos = backslashPos+1; else if ( forwardSlashPos != string::npos ) beginBaseFileNamePos = forwardSlashPos+1; else beginBaseFileNamePos = 0; std::string baseFileName = ___1394.substr(beginBaseFileNamePos, string::npos); ENSURE(!baseFileName.empty()); return baseFileName; } }}
