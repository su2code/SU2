void InitHashTables ();
int ParseGLFunc _ANSI_ARGS_((Tcl_Interp *interp,
			     int argc,
			     char *argv [],
			     int *nArg));
int SearchEnumVal _ANSI_ARGS_((Tcl_Interp* interp,
			       char * name,
			       GLenum* val));
int SearchEnumName _ANSI_ARGS_((Tcl_Interp* interp,
				GLenum val,
				char ** name));



