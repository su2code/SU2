import sys,os
from distutils.core import setup, Extension
from distutils import sysconfig


def setup_amgio(argv=[]):
    
    sav_argv = sys.argv;    
    sys.argv = ['', 'build_ext', '--inplace'];
    
    working_dir = os.getcwd();
    file_dir    = os.path.dirname(os.path.abspath(__file__));
    
    os.chdir(file_dir);
    ext_opts = {'extra_compile_args': ['-Wl,-soname,units.so', '-Isrc']}
    setup(name = '_amgio',
          package_data={'': ['_amgio.so']},
          ext_modules=[ \
          Extension("_amgio",
            sources=[ "./_amgio/amgio_py_wrap.c", \
                      "./_amgio/amgio_py.c", \
                      "./_amgio/mesh.c", \
                      "./_amgio/GMFio.c", \
                      "./_amgio/SU2io.c", \
                      "./_amgio/option.c", \
                      "./_amgio/libmesh6.c", \
                      "./_amgio/amgio_py.i", \
                      "./_amgio/convert.c"],
            extra_compile_args=["-std=c99",
                                "-Wno-unused-variable",
                                "-Wno-unused-result",
                                "-Wl,-soname,_amgio.so"]),
          ],);

    # os.rename("_amgio"+sysconfig.get_config_var('EXT_SUFFIX'),"_amgio.so")
    os.chdir(working_dir);
    sys.argv = sav_argv;

if __name__ == '__main__':
    setup_amgio(sys.argv)