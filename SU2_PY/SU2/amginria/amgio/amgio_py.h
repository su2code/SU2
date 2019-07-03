
int py_ConvertSU2toInria( char *MshNam, char *SolNam, char *OutNam ) ;
int py_ConvertInriatoSU2( char *MshNam, char *SolNam, char *OutNam ) ;
int py_SplitSolution(char *SolNam, int dim, char *prefix, char *adap_sensor);

//void py_ReadMesh (char *MshNam, char *SolNam, PyObject *pyVer, PyObject *pyTri, PyObject *pyTet, PyObject *pyEdg, PyObject *pySol, PyObject *pySolHeader, PyObject *pyMarkers);
void py_ReadMesh (char *MshNam, char *SolNam, PyObject *pyVer, PyObject *pyTri, PyObject *pyTet, PyObject *pyEdg, PyObject *pyHex, 
PyObject *pyQua, PyObject *pyPyr, PyObject *pyPri, 
 PyObject *pySol, PyObject *pySolHeader,  PyObject *pyMarkers);

void py_WriteMesh(char *MshNam, char *SolNam, PyObject *pyVer, PyObject *pyTri, PyObject *pyTet, PyObject *pyEdg,  PyObject *pyHex, 
 PyObject *pyQua, PyObject *pyPyr, PyObject *pyPri, PyObject *pySol, PyObject *pyMarkers, int Dim);

//void py_WriteMesh(char *MshNam, char *SolNam, PyObject *pyVer, PyObject *pyTri, PyObject *pyTet, PyObject *pyEdg, PyObject *pySol, PyObject *pyMarkers, int Dim);
void py_WriteSolution(char *SolNam, PyObject *pyVer, PyObject *pySol, PyObject *pySolHeader, int NbrVer, int Dim);
