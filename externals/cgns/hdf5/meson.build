hdf5_default_warnings = []

cc = meson.get_compiler('c')

if target_machine.system() != 'windows' and cc.get_id() != 'intel'
  desired_warnings = ['-Wno-unused-result',
                      '-Wno-unused-parameter',
                      '-Wno-unused-variable',
                      '-Wno-unused-but-set-variable',
                      '-Wno-sign-compare',
                      '-Wno-stringop-overflow',
                      '-Wno-error=discarded-qualifiers',
                      '-Wno-error=maybe-uninitialized',
                      '-Wno-error=unused-function',
                      '-Wno-error=unused-label',
                      '-Wno-pedantic',
                      '-Wno-error=stringop-truncation',
                      '-Wno-error=cast-function-type',
                      '-Wno-stringop-truncation']
  foreach flag : desired_warnings
    if cc.has_argument(flag)
      hdf5_default_warnings += flag
    endif
  endforeach
endif

# Check if we can run tests
cc_can_run = not meson.is_cross_build() or meson.has_exe_wrapper()


hdf5_conf_data = configuration_data()
hdf5_conf_data.set('version', '1.12.1')

ext_deps = []

opt_zlib = dependency('zlib', required: false, static: get_option('static-cgns-deps'))
opt_szip = dependency('szip', required: false, static: get_option('static-cgns-deps'))

if opt_zlib.found() and cc.has_header('zlib.h')
  if cc.has_function('inflate', dependencies: opt_zlib, prefix: '#include <zlib.h>')
    hdf5_conf_data.set('H5_HAVE_FILTER_DEFLATE', 1)
    hdf5_conf_data.set('H5_HAVE_ZLIB_H', 1)
    hdf5_conf_data.set('H5_HAVE_LIBZ', 1)
    ext_deps += opt_zlib
  endif
endif

if opt_szip.found()
  hdf5_conf_data.set('H5_HAVE_FILTER_SZIP', 1)
  hdf5_conf_data.set('H5_HAVE_SZLIB_H', 1)
  hdf5_conf_data.set('H5_HAVE_LIBSZ', 1)
  hdf5_conf_data.set('H5_HAVE_SZIP_ENCODER', 1)
  ext_deps += opt_szip
endif


if target_machine.system() == 'windows'
  hdf5_conf_data.set('H5_HAVE_WINDOWS', 1)
  hdf5_conf_data.set('H5_HAVE_WINDOW_PATH', 1)
  hdf5_conf_data.set('H5_HAVE_MINGW', 1)
  hdf5_conf_data.set('H5_HAVE_WIN32_API', 1)
  #x
  hdf5_conf_data.set('H5_HAVE_LIBM', 1)
  hdf5_conf_data.set('H5_HAVE_STRDUP', 1)
  hdf5_conf_data.set('H5_HAVE_SYSTEM', 1)
  hdf5_conf_data.set('H5_HAVE_LONGJMP', 1)

  #
  hdf5_conf_data.set('H5_HAVE_GETCONSOLESCREENBUFFERINFO', 1)
  hdf5_conf_data.set('H5_GETTIMEOFDAY_GIVES_TZ', 1)
  hdf5_conf_data.set('H5_HAVE_TIMEZONE', 1)
  hdf5_conf_data.set('H5_HAVE_GETTIMEOFDAY', 1)
  hdf5_conf_data.set('H5_HAVE_LIBWS2_32', 1)
  hdf5_conf_data.set('H5_HAVE_LIBWSOCK32', 1)
endif

if target_machine.system() != 'windows'
  if target_machine.endian() == 'big'
    hdf5_conf_data.set('H5_WORDS_BIGENDIAN', 1)
  else
    hdf5_conf_data.set('H5_WORDS_BIGENDIAN', 0)
  endif
endif


if target_machine.system() == 'windows'
  hdf5_conf_data.set('H5_DEFAULT_PLUGINDIR', '%ALLUSERSPROFILE%\\\\hdf5\\\\lib\\\\plugin')
else
  hdf5_conf_data.set('H5_DEFAULT_PLUGINDIR', '/usr/local/hdf5/lib/plugin')
endif

if target_machine.system()=='darwin'
  hdf5_conf_data.set('H5_HAVE_DARWIN', 1)
endif

#if cc.has_function('difftime', prefix : '#include <time.h>')
#  # function exists, do whatever is required.
#endif



if cc.has_header('sys/file.h')
  hdf5_conf_data.set('H5_HAVE_SYS_FILE_H', 1)
endif

if cc.has_header('sys/ioctl.h')
  hdf5_conf_data.set('H5_HAVE_SYS_IOCTL_H', 1)
endif

if cc.has_header('sys/resource.h')
  hdf5_conf_data.set('H5_HAVE_SYS_RESOURCE_H', 1)
endif

if cc.has_header('sys/socket.h')
  hdf5_conf_data.set('H5_HAVE_SYS_SOCKET_H', 1)
endif

if cc.has_header('sys/stat.h')
  hdf5_conf_data.set('H5_HAVE_SYS_STAT_H', 1)
endif

if cc.has_header('sys/time.h')
  hdf5_conf_data.set('H5_HAVE_SYS_TIME_H', 1)
endif

if cc.has_header('sys/types.h')
  hdf5_conf_data.set('H5_HAVE_SYS_TYPES_H', 1)
endif

if cc.has_header('features.h')
  hdf5_conf_data.set('H5_HAVE_FEATURES_H', 1)
endif

if cc.has_header('dirent.h')
  hdf5_conf_data.set('H5_HAVE_DIRENT_H', 1)
endif

if cc.has_header('setjmp.h')
  hdf5_conf_data.set('H5_HAVE_SETJMP_H', 1)
endif

if cc.has_header('stddef.h')
  hdf5_conf_data.set('H5_HAVE_STDDEF_H', 1)
endif

if cc.has_header('stdint.h')
  hdf5_conf_data.set('H5_HAVE_STDINT_H', 1)
endif

if cc.has_header('unistd.h')
  hdf5_conf_data.set('H5_HAVE_UNISTD_H', 1)
endif

if cc.has_header('io.h')
  hdf5_conf_data.set('H5_HAVE_IO_H', 1)
endif

if target_machine.system() != 'windows'
  hdf5_conf_data.set('H5_HAVE_WINSOCK2_H', cc.has_header('winsock2.h'))
endif


if cc.has_header('globus/common.h')
  hdf5_conf_data.set('H5_HAVE_GLOBUS_COMMON_H', 1)
endif
if cc.has_header('pdb.h')
  hdf5_conf_data.set('H5_HAVE_PDB_H', 1)
endif
if cc.has_header('pthread.h')
  hdf5_conf_data.set('H5_HAVE_PTHREAD_H', 1)
endif

if cc.has_header('srbclient.h')
  hdf5_conf_data.set('H5_HAVE_SRBCLIENT_H', 1)
endif
if cc.has_header('string.h')
  hdf5_conf_data.set('H5_HAVE_STRING_H', 1)
endif
if cc.has_header('strings.h')
  hdf5_conf_data.set('H5_HAVE_STRINGS_H', 1)
endif
if cc.has_header('stdlib.h')
  hdf5_conf_data.set('H5_HAVE_STDLIB_H', 1)
endif
if cc.has_header('memory.h')
  hdf5_conf_data.set('H5_HAVE_MEMORY_H', 1)
endif
if cc.has_header('dlfcn.h')
  hdf5_conf_data.set('H5_HAVE_DLFCN_H', 1)
endif
if cc.has_header('inttypes.h')
  hdf5_conf_data.set('H5_HAVE_INTTYPES_H', 1)
endif
if cc.has_header('netinet/in.h')
  hdf5_conf_data.set('H5_HAVE_NETINET_IN_H', 1)
endif
if cc.has_header('netdb.h')
  hdf5_conf_data.set('H5_HAVE_NETDB_H', 1)
endif
if cc.has_header('arpa/inet.h')
  hdf5_conf_data.set('H5_HAVE_ARPA_INET_H', 1)
endif

if cc.has_header('stdbool.h')
  hdf5_conf_data.set('H5_HAVE_STDBOOL_H', 1)
endif

if cc.has_header('quadmath.h')
  hdf5_conf_data.set('H5_HAVE_QUADMATH_H', 1)
else
  hdf5_conf_data.set('H5_HAVE_QUADMATH_H', 0)
endif


libm = []
if target_machine.system() != 'windows'
  libm = [ cc.find_library('m') ]
  if not(cc.has_function('ceil', dependencies: libm, prefix: '#include <math.h>'))
    hdf5_conf_data.set('H5_HAVE_LIBM', 0)
  else
    hdf5_conf_data.set('H5_HAVE_LIBM', 1)
    ext_deps += libm
  endif
endif

libdl = []
if target_machine.system() != 'windows'
  libdl = [ cc.find_library('dl') ]
  if not(cc.has_function('dlopen', dependencies: libdl, prefix: '#include <dlfcn.h>'))
    hdf5_conf_data.set('H5_HAVE_LIBDL', 0)
  else
    hdf5_conf_data.set('H5_HAVE_LIBDL', 1)
    ext_deps += libdl
  endif
endif



libws2_32 = []
if target_machine.system() == 'windows'
  libws2_32 = [ cc.find_library('ws2_32') ]
  if not(cc.has_function('WSAStartup', dependencies: libws2_32, prefix: '#include <winsock2.h>'))
    hdf5_conf_data.set('H5_HAVE_LIBWS2_32', 0)
  else
    hdf5_conf_data.set('H5_HAVE_LIBWS2_32', 1)
  endif
endif


libwsock32 = []
if target_machine.system() == 'windows'
  libwsock32 = [ cc.find_library('wsock32') ]
  if not(cc.has_function('gethostbyname', dependencies: libwsock32, prefix: '#include <winsock.h>'))
    hdf5_conf_data.set('H5_HAVE_LIBWSOCK32', 0)
  else
    hdf5_conf_data.set('H5_HAVE_LIBWSOCK32', 1)
  endif
endif


libucb = [ cc.find_library('ucb', required:false) ]
if not(cc.has_function('gethostbyname', dependencies: libucb, prefix: '#include <unistd.h>'))
  hdf5_conf_data.set('H5_HAVE_LIBUCB', 0)
else
  hdf5_conf_data.set('H5_HAVE_LIBUCB', 1)
endif

defs = []
#if  not windows or mingw
if target_machine.system()!='darwin'
  defs += [
        '-D_FILE_OFFSET_BITS=64',
        '-D_LARGEFILE64_SOURCE',
        '-D_LARGEFILE_SOURCE',]

  defs += ['-D_POSIX_C_SOURCE=200809L', '-D_GNU_SOURCE',]
endif
#endif

if not(cc.has_function('lseek64', prefix: '''#include <sys/types.h>
#include <unistd.h>'''))
  hdf5_conf_data.set('H5_HAVE_LSEEK64', 0)
else
  hdf5_conf_data.set('H5_HAVE_LSEEK64', 1)
endif


if not(cc.has_function('fseeko', prefix: '''#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>'''))
  hdf5_conf_data.set('H5_HAVE_FSEEKO', 0)
else
  hdf5_conf_data.set('H5_HAVE_FSEEKO', 1)
endif

if not(cc.has_function('stat64', prefix: '''#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>'''))
  hdf5_conf_data.set('H5_HAVE_STAT64', 0)
else
  hdf5_conf_data.set('H5_HAVE_STAT64', 1)
endif






####
####

hdf5_conf_data.set('H5_SIZEOF_CHAR', cc.sizeof('char'))

hdf5_conf_data.set('H5_SIZEOF_SHORT', cc.sizeof('short'))

hdf5_conf_data.set('H5_SIZEOF_INT', cc.sizeof('int'))

hdf5_conf_data.set('H5_SIZEOF_UNSIGNED', cc.sizeof('unsigned'))


if target_machine.system() != 'darwin'
  hdf5_conf_data.set('H5_SIZEOF_LONG', cc.sizeof('long'))
endif


hdf5_conf_data.set('H5_SIZEOF_LONG_LONG', cc.sizeof('long long'))
hdf5_conf_data.set('H5_SIZEOF___INT64', cc.sizeof('__int64'))
if hdf5_conf_data.get('H5_SIZEOF___INT64') < 1
  hdf5_conf_data.set('H5_SIZEOF___INT64', 0)
endif

hdf5_conf_data.set('H5_SIZEOF_FLOAT',            cc.sizeof('float'))
hdf5_conf_data.set('H5_SIZEOF_DOUBLE',           cc.sizeof('double'))
hdf5_conf_data.set('H5_SIZEOF_LONG_DOUBLE',      cc.sizeof('long double'))

hdf5_conf_data.set('H5_SIZEOF_INT8_T',           cc.sizeof ('int8_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_UINT8_T',          cc.sizeof ('uint8_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_INT_LEAST8_T',     cc.sizeof ('int_least8_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_UINT_LEAST8_T',    cc.sizeof ('uint_least8_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_INT_FAST8_T',      cc.sizeof ('int_fast8_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_UINT_FAST8_T',     cc.sizeof ('uint_fast8_t', prefix:'#include <stdint.h>'))

hdf5_conf_data.set('H5_SIZEOF_INT16_T',          cc.sizeof ('int16_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_UINT16_T',         cc.sizeof ('uint16_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_INT_LEAST16_T',    cc.sizeof ('int_least16_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_UINT_LEAST16_T',   cc.sizeof ('uint_least16_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_INT_FAST16_T',     cc.sizeof ('int_fast16_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_UINT_FAST16_T',    cc.sizeof ('uint_fast16_t', prefix:'#include <stdint.h>'))

hdf5_conf_data.set('H5_SIZEOF_INT32_T',          cc.sizeof ('int32_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_UINT32_T',         cc.sizeof ('uint32_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_INT_LEAST32_T',    cc.sizeof ('int_least32_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_UINT_LEAST32_T',   cc.sizeof ('uint_least32_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_INT_FAST32_T',     cc.sizeof ('int_fast32_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_UINT_FAST32_T',    cc.sizeof ('uint_fast32_t', prefix:'#include <stdint.h>'))

hdf5_conf_data.set('H5_SIZEOF_INT64_T',          cc.sizeof ('int64_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_UINT64_T',         cc.sizeof ('uint64_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_INT_LEAST64_T',    cc.sizeof ('int_least64_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_UINT_LEAST64_T',   cc.sizeof ('uint_least64_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_INT_FAST64_T',     cc.sizeof ('int_fast64_t', prefix:'#include <stdint.h>'))
hdf5_conf_data.set('H5_SIZEOF_UINT_FAST64_T',    cc.sizeof ('uint_fast64_t', prefix:'#include <stdint.h>'))
if target_machine.system() != 'darwin'
  hdf5_conf_data.set('H5_SIZEOF_SIZE_T', cc.sizeof('size_t'))
  hdf5_conf_data.set('H5_SIZEOF_SSIZE_T', cc.sizeof('ssize_t'))
  if hdf5_conf_data.get('H5_SIZEOF_SSIZE_T') < 1
    hdf5_conf_data.set('H5_SIZEOF_SSIZE_T', 0)
  endif
  hdf5_conf_data.set('H5_SIZEOF_PTRDIFF_T', cc.sizeof('ptrdiff_t'))
endif

hdf5_conf_data.set('H5_SIZEOF_OFF_T', cc.sizeof('off_t', prefix: '#include <sys/types.h>'))
hdf5_conf_data.set('H5_SIZEOF_OFF64_T', cc.sizeof('off64_t', prefix: '#include <sys/types.h>'))
if hdf5_conf_data.get('H5_SIZEOF_OFF64_T') < 1
  hdf5_conf_data.set('H5_SIZEOF_OFF64_T', 0)
endif

hdf5_conf_data.set('H5_SIZEOF_TIME_T', cc.sizeof('time_t', prefix: '#include <sys/time.h>'))

if hdf5_conf_data.get('H5_HAVE_STDBOOL_H') == 1
  hdf5_conf_data.set('H5_SIZEOF_BOOL', cc.sizeof('bool', prefix: '#include <stdbool.h>'))
else
  hdf5_conf_data.set('H5_SIZEOF_BOOL', cc.sizeof('_Bool'))
endif

if hdf5_conf_data.get('H5_HAVE_QUADMATH_H') == 1
  hdf5_conf_data.set('H5_SIZEOF__QUAD', cc.sizeof('_Quad', prefix: '#include <quadmath.h>'))
  if hdf5_conf_data.get('H5_SIZEOF__QUAD') < 1
    hdf5_conf_data.set('H5_SIZEOF__QUAD', 0)
  endif
  hdf5_conf_data.set('H5_SIZEOF___FLOAT128', cc.sizeof('__float128', prefix: '#include <quadmath.h>'))
  if hdf5_conf_data.get('H5_SIZEOF___FLOAT128') < 1
    hdf5_conf_data.set('H5_SIZEOF___FLOAT128', 0)
  endif
else
  hdf5_conf_data.set('H5_SIZEOF__QUAD', 0)
  hdf5_conf_data.set('H5_SIZEOF___FLOAT128', 0)
endif

code_dev_t ='''
#include <sys/types.h>
int main ()
{
  dev_t d1, d2;
  if(d1==d2)
    return 0;
  return 1;
}
'''

result = cc.compiles(code_dev_t, name : 'check dev_t')
if result and cc_can_run
  result_dev = cc.run(code_dev_t, name : 'check dev_t scalar')
  hdf5_conf_data.set('H5_DEV_T_IS_SCALAR', result_dev.returncode())
else
  hdf5_conf_data.set('H5_DEV_T_IS_SCALAR', 1)
endif

#_POSIX_C_SOURCE couls possibly be updated to 200112L
hdf5_conf_data.set('HAVE_CLOCK_MONOTONIC', cc.compiles('''
       #define _POSIX_C_SOURCE 200809L
       #include <time.h>
       #include <unistd.h>
       #ifndef CLOCK_MONOTONIC
       #error CLOCK_MONOTONIC not defined
       #endif
   ''',
       name: 'CLOCK_MONOTONIC'))


if cc.has_member('struct tm', 'tm_gmtoff', prefix : '#include<time.h>')
  hdf5_conf_data.set('H5_HAVE_TM_GMTOFF', 1)
endif

if cc.has_member('struct tm', '__tm_gmtoff', prefix : '#include<time.h>')
  hdf5_conf_data.set('H5_HAVE___TM_GMTOFF', 1)
endif

if cc.has_member('struct tm', 'tm_sec', prefix : '''#include<sys/types.h>
  #include<sys/time.h>
  #include<time.h>''')
  hdf5_conf_data.set('H5_TIME_WITH_SYS_TIME', 1)
endif

if hdf5_conf_data.get('H5_HAVE_SYS_TIME_H') == 1
  if cc.has_member('struct tm', 'tz_minuteswest', prefix : '''#include<sys/types.h>
  #include<sys/time.h>
  #include<time.h>''')
    hdf5_conf_data.set('H5_HAVE_STRUCT_TIMEZONE', 1)
  endif
else
  if cc.has_member('struct tm', 'tz_minuteswest', prefix : '''#include<sys/types.h>
  #include<time.h>''')
    hdf5_conf_data.set('H5_HAVE_STRUCT_TIMEZONE', 1)
  endif
endif

if cc.has_function('gettimeofday')
  hdf5_conf_data.set('H5_HAVE_GETTIMEOFDAY', 1)
endif

code_gettimeofday_give_tz='''
#include <time.h>
#include <sys/time.h>
int main(void)
{
 struct timeval tv;
 struct timezone tz;
 tz.tz_minuteswest = 7777;  /* Initialize to an unreasonable number */
 tz.tz_dsttime = 7;
 gettimeofday(&tv, &tz);
    /* Check whether the function returned any value at all */
 if(tz.tz_minuteswest == 7777 && tz.tz_dsttime == 7)
     return 1;
 else return 0;
}
'''

hdf5_conf_data.set('H5_GETTIMEOFDAY_GIVE_TZ', cc.compiles(code_gettimeofday_give_tz))
hdf5_conf_data.set('H5_HAVE_TM_ZONE', 1)
hdf5_conf_data.set('H5_HAVE_STRUCT_TM_TM_ZONE', 1)

if target_machine.system() != 'windows'
  hdf5_conf_data.set('H5_HAVE_TIMEZONE', 1)
endif

# Have non POSIX st_blocks ?
if cc.has_member('struct stat', 'st_blocks', prefix : '''#include<sys/types.h>
  #include<sys/stat.h>''')
  hdf5_conf_data.set('H5_HAVE_STAT_ST_BLOCKS', 1)
endif

if cc.has_function('ioctl', prefix: '#include <sys/ioctl.h>')
  hdf5_conf_data.set('H5_HAVE_IOCTL', 1)
endif

if cc.has_header_symbol('sys/ioctl.h', 'TIOCGWINSZ')
  hdf5_conf_data.set('H5_HAVE_TIOCGWINSZ', 1)
endif
if cc.has_header_symbol('sys/ioctl.h', 'TIOCGETD')
  hdf5_conf_data.set('H5_HAVE_TIOCGETD', 1)
endif

if cc.has_function('getpwuid')
  hdf5_conf_data.set('H5_HAVE_GETPWUID', 1)
endif

#-----------------------------------------------------------------------------
# Check for some functions that are used
#
hdf5_conf_data.set('H5_HAVE_ALARM',           cc.has_function ('alarm') )
hdf5_conf_data.set('H5_HAVE_FCNTL',           cc.has_function ('fcntl')          )
hdf5_conf_data.set('H5_HAVE_FLOCK',           cc.has_function ('flock')         )
hdf5_conf_data.set('H5_HAVE_FORK',            cc.has_function ('fork')            )
hdf5_conf_data.set('H5_HAVE_FREXPF',          cc.has_function ('frexpf')        )
hdf5_conf_data.set('H5_HAVE_FREXPL',          cc.has_function ('frexpl')        )

hdf5_conf_data.set('H5_HAVE_GETHOSTNAME',     cc.has_function ('gethostname') )
hdf5_conf_data.set('H5_HAVE_GETRUSAGE',       cc.has_function ('getrusage') )
hdf5_conf_data.set('H5_HAVE_LLROUND',         cc.has_function ('llround') )
hdf5_conf_data.set('H5_HAVE_LLROUNDF',        cc.has_function ('llroundf') )
hdf5_conf_data.set('H5_HAVE_LROUND',          cc.has_function ('lround') )
hdf5_conf_data.set('H5_HAVE_LROUNDF',         cc.has_function ('lroundf'))
hdf5_conf_data.set('H5_HAVE_LSTAT',           cc.has_function ('lstat'))

hdf5_conf_data.set('H5_HAVE_PREAD',           cc.has_function ('pread') )
hdf5_conf_data.set('H5_HAVE_PWRITE',          cc.has_function ('pwrite'))
hdf5_conf_data.set('H5_HAVE_RAND_R',          cc.has_function ('rand_r'))
hdf5_conf_data.set('H5_HAVE_RANDOM',          cc.has_function ('random'))
hdf5_conf_data.set('H5_HAVE_ROUND',           cc.has_function ('round') )
hdf5_conf_data.set('H5_HAVE_ROUNDF',          cc.has_function ('roundf'))
hdf5_conf_data.set('H5_HAVE_SETSYSINFO',      cc.has_function ('setsysinfo'))

hdf5_conf_data.set('H5_HAVE_SIGNAL',          cc.has_function ('signal')         )
hdf5_conf_data.set('H5_HAVE_LONGJMP',         cc.has_function ('longjmp')        )
hdf5_conf_data.set('H5_HAVE_SETJMP',          cc.has_function ('setjmp')         )
hdf5_conf_data.set('H5_HAVE_SIGLONGJMP',      cc.has_function ('siglongjmp')     )
hdf5_conf_data.set('H5_HAVE_SIGSETJMP',       cc.has_function ('sigsetjmp')      )
hdf5_conf_data.set('H5_HAVE_SIGPROCMASK',     cc.has_function ('sigprocmask')    )

hdf5_conf_data.set('H5_HAVE_SNPRINTF',        cc.has_function ('snprintf')       )
hdf5_conf_data.set('H5_HAVE_SRANDOM',         cc.has_function ('srandom')        )
hdf5_conf_data.set('H5_HAVE_STRDUP',          cc.has_function ('strdup')         )
hdf5_conf_data.set('H5_HAVE_STRTOLL',         cc.has_function ('strtoll')        )
hdf5_conf_data.set('H5_HAVE_STRTOULL',        cc.has_function ('strtoull')      )
hdf5_conf_data.set('H5_HAVE_SYMLINK',         cc.has_function ('symlink')        )
hdf5_conf_data.set('H5_HAVE_SYSTEM',          cc.has_function ('system')         )

hdf5_conf_data.set('H5_HAVE_TMPFILE',         cc.has_function ('tmpfile')        )
hdf5_conf_data.set('H5_HAVE_ASPRINTF',        cc.has_function ('asprintf')       )
hdf5_conf_data.set('H5_HAVE_VASPRINTF',       cc.has_function ('vasprintf')      )
hdf5_conf_data.set('H5_HAVE_WAITPID',         cc.has_function ('waitpid')        )

hdf5_conf_data.set('H5_HAVE_VSNPRINTF',       cc.has_function ('vsnprintf')      )

if cc.has_header_symbol('setjmp.h', 'sigsetjmp')
  hdf5_conf_data.set('H5_HAVE_MACRO_SIGSETJMP', 1)
  hdf5_conf_data.set('H5_HAVE_SIGSETJMP', 1)
endif

#-----------------------------------------------------------------------------
# Check a bunch of other functions TODO
#-----------------------------------------------------------------------------

hdf5_conf_data.set('H5_STDC_HEADERS', 1)
hdf5_conf_data.set('H5_HAVE_ATTRIBUTE', 1)
hdf5_conf_data.set('H5_HAVE_C99_FUNC', 1)
hdf5_conf_data.set('H5_HAVE_FUNCTION', 1)
hdf5_conf_data.set('H5_HAVE_C99_DESIGNATED_INITIALIZER', 1)
hdf5_conf_data.set('H5_HAVE_SYSTEM_SCOPE_THREADS', 1)
hdf5_conf_data.set('H5_HAVE_SOCKLEN_T', 1)
#

# Check IOEO
if target_machine.system() == 'windows'
  code_ioeo='''
#include <windows.h>
typedef void (WINAPI *PGNSI)(LPSYSTEM_INFO);
int main ()
{
    PGNSI pGNSI;
    pGNSI = (PGNSI) GetProcAddress(
      GetModuleHandle(TEXT("kernel32.dll")),
      "InitOnceExecuteOnce");
    if(NULL == pGNSI)
        return 1;
    else
        return 0;
}
'''
  result = cc.compiles(code_ioeo, name : 'check ioeo')
  if result and cc_can_run
    result_ioeo = cc.run(code_ioeo, name : 'check ioeo works')
    hdf5_conf_data.set('H5_HAVE_IOEO', result_ioeo.returncode())
  else
    hdf5_conf_data.set('H5_HAVE_IOEO', 1)
  endif

endif

# Check inline
if not cc.compiles('inline void func(); void func() { } int main() { func(); return 0; }', name: 'inline keyword check')

  if cc.compiles('__inline void func(); void func() { } int main() { func(); return 0; }', name: '__inline keyword check')
    hdf5_conf_data.set('H5_HAVE___INLINE', 1)
  else
    if cc.compiles('__inline__ void func(); void func() { } int main() { func(); return 0; }', name: '__inline__ keyword check')
      hdf5_conf_data.set('H5_HAVE__INLINE__', 1)
    endif
  endif
else
  hdf5_conf_data.set('H5_HAVE_INLINE', 1)
endif

# TODO check
hdf5_conf_data.set('H5_PRINTF_LL_WIDTH', '"ll"')


hdf5_conf_data.set('H5_CONVERT_DENORMAL_FLOAT', 1)


hdf5_conf_data.set('H5_DEFAULT_VFD', 'H5FD_SEC2')

hdf5_conf_data.set('H5_WANT_DATA_ACCURACY', 1)
hdf5_conf_data.set('H5_WANT_DCONV_EXCEPTION', 1)
hdf5_conf_data.set('H5_USE_FILE_LOCKING', 1)
hdf5_conf_data.set('H5_IGNORE_DISABLED_FILE_LOCKS', 1)
hdf5_conf_data.set('H5_HAVE_LARGE_HSIZET', 1)

if mpi
  hdf5_conf_data.set('H5_HAVE_PARALLEL', 1)
  ext_deps += mpi_dep
else
  threads_dep = dependency('threads')
  ext_deps += threads_dep
  #hdf5_conf_data.set('H5_HAVE_THREADSAFE', 1) # could be 1 if linux
endif

hdf5_conf_data.set('H5_USE_112_API_DEFAULT', 1)

#PACKAGE
hdf5_conf_data.set('HDF5_PACKAGE', 'hdf5')
hdf5_conf_data.set('HDF5_PACKAGE_NAME', 'HDF5')
hdf5_conf_data.set('HDF5_PACKAGE_VERSION_STRING', '1.12.1')
hdf5_conf_data.set('HDF5_PACKAGE_TARNAME', 'hdf5')
hdf5_conf_data.set('HDF5_PACKAGE_URL', '')
hdf5_conf_data.set('HDF5_PACKAGE_STRING', 'HDF5 1.12.1')
hdf5_conf_data.set('HDF5_PACKAGE_BUGREPORT', 'help@hdfgroup.com')


# No Fortran
hdf5_conf_data.set('CMAKE_Fortran_COMPILER_ID', 'GNU')
hdf5_conf_data.set('H5_FC_FUNC', 'H5_FC_FUNC')
hdf5_conf_data.set('H5_FC_FUNC_', 'H5_FC_FUNC_')
hdf5_conf_data.set('H5_FORTRAN_C_LONG_DOUBLE_IS_UNIQUE', '1')
hdf5_conf_data.set('H5_FORTRAN_HAVE_C_LONG_DOUBLE', '1')
hdf5_conf_data.set('H5_PAC_FC_MAX_REAL_PRECISION', 36)

# Check TODO
#should be determined through LDBL_DIG ...
hdf5_conf_data.set('H5_PAC_C_MAX_REAL_PRECISION', 36)

hdf5_conf_data.set('H5_LDOUBLE_TO_LONG_SPECIAL', 1)
hdf5_conf_data.set('H5_LONG_TO_LDOUBLE_SPECIAL', 1)
hdf5_conf_data.set('H5_LDOUBLE_TO_LLONG_ACURRATE', 1)
hdf5_conf_data.set('H5_LLONG_TO_LDOUBLE_CORRECT', 1)
hdf5_conf_data.set('H5_DISABLE_SOME_LDOUBLE_CONV', 1)
hdf5_conf_data.set('H5_NO_ALIGNEMENT_RESTRICTIONS', 1)




configure_file(input: 'H5pubconf.h.in',
  output: 'H5pubconf.h',
  format: 'cmake@',
  configuration : hdf5_conf_data
  )


hdf5_include = include_directories('./')



H5_SOURCES = [
     'H5.c',
     'H5checksum.c',
     'H5dbg.c',
     'H5mpi.c',
     'H5system.c',
     'H5timer.c',
     'H5trace.c' ]

H5A_SOURCES = [
    'H5A.c',
    'H5Abtree2.c',
    'H5Adense.c',
    'H5Adeprec.c',
    'H5Aint.c',
    'H5Atest.c']

H5AC_SOURCES = [
    'H5AC.c',
    'H5ACdbg.c',
    'H5ACmpio.c',
    'H5ACproxy_entry.c']

H5B_SOURCES = [
 'H5B.c',
 'H5Bcache.c',
 'H5Bdbg.c',
]


H5B2_SOURCES = [
    'H5B2.c',
    'H5B2cache.c',
    'H5B2dbg.c',
    'H5B2hdr.c',
    'H5B2int.c',
    'H5B2internal.c',
    'H5B2leaf.c',
    'H5B2stat.c',
    'H5B2test.c']

H5C_SOURCES = [
    'H5C.c',
    'H5Cdbg.c',
    'H5Cepoch.c',
    'H5Cimage.c',
    'H5Clog.c',
    'H5Clog_json.c',
    'H5Clog_trace.c',
    'H5Cmpio.c',
    'H5Cprefetched.c',
    'H5Cquery.c',
    'H5Ctag.c',
    'H5Ctest.c' ]

H5CS_SOURCES = [
    'H5CS.c' ]

H5CX_SOURCES = [ 'H5CX.c' ]


H5D_SOURCES = [
    'H5D.c',
    'H5Dbtree.c',
    'H5Dbtree2.c',
    'H5Dchunk.c',
    'H5Dcompact.c',
    'H5Dcontig.c',
    'H5Ddbg.c',
    'H5Ddeprec.c',
    'H5Dearray.c',
    'H5Defl.c',
    'H5Dfarray.c',
    'H5Dfill.c',
    'H5Dint.c',
    'H5Dio.c',
    'H5Dlayout.c',
    'H5Dmpio.c',
    'H5Dnone.c',
    'H5Doh.c',
    'H5Dscatgath.c',
    'H5Dselect.c',
    'H5Dsingle.c',
    'H5Dtest.c',
    'H5Dvirtual.c']

H5E_SOURCES = [
   'H5E.c',
   'H5Edeprec.c',
   'H5Eint.c' ]

H5EA_SOURCES = [
    'H5EA.c',
    'H5EAcache.c',
    'H5EAdbg.c',
    'H5EAdblkpage.c',
    'H5EAdblock.c',
    'H5EAhdr.c',
    'H5EAiblock.c',
    'H5EAint.c',
    'H5EAsblock.c',
    'H5EAstat.c',
    'H5EAtest.c' ]

H5F_SOURCES = [
    'H5F.c',
    'H5Faccum.c',
    'H5Fcwfs.c',
    'H5Fdbg.c',
    'H5Fdeprec.c',
    'H5Fefc.c',
    'H5Ffake.c',
    'H5Fint.c',
    'H5Fio.c',
    'H5Fmount.c',
    'H5Fmpi.c',
    'H5Fquery.c',
    'H5Fsfile.c',
    'H5Fspace.c',
    'H5Fsuper.c',
    'H5Fsuper_cache.c',
    'H5Ftest.c' ]

H5FA_SOURCES = [
    'H5FA.c',
    'H5FAcache.c',
    'H5FAdbg.c',
    'H5FAdblkpage.c',
    'H5FAdblock.c',
    'H5FAhdr.c',
    'H5FAint.c',
    'H5FAstat.c',
    'H5FAtest.c'
    ]

H5FD_SOURCES = [
    'H5FD.c',
    'H5FDcore.c',
    'H5FDdirect.c',
    'H5FDfamily.c',
    'H5FDhdfs.c',
    'H5FDint.c',
    'H5FDlog.c',
    'H5FDmirror.c',
    'H5FDmpi.c',
    'H5FDmpio.c',
    'H5FDmulti.c',
    'H5FDros3.c',
    'H5FDs3comms.c',
    'H5FDsec2.c',
    'H5FDspace.c',
    'H5FDsplitter.c',
    'H5FDstdio.c',
    'H5FDtest.c',
    'H5FDwindows.c']

H5FL_SOURCES = [ 'H5FL.c' ]

H5FO_SOURCES = [ 'H5FO.c' ]

H5FS_SOURCES = [
    'H5FS.c',
    'H5FScache.c',
    'H5FSdbg.c',
    'H5FSint.c',
    'H5FSsection.c',
    'H5FSstat.c',
    'H5FStest.c' ]

H5G_SOURCES =[
    'H5G.c',
    'H5Gbtree2.c',
    'H5Gcache.c',
    'H5Gcompact.c',
    'H5Gdense.c',
    'H5Gdeprec.c',
    'H5Gent.c',
    'H5Gint.c',
    'H5Glink.c',
    'H5Gloc.c',
    'H5Gname.c',
    'H5Gnode.c',
    'H5Gobj.c',
    'H5Goh.c',
    'H5Groot.c',
    'H5Gstab.c',
    'H5Gtest.c',
    'H5Gtraverse.c']

H5HF_SOURCES = [
    'H5HF.c',
    'H5HFbtree2.c',
    'H5HFcache.c',
    'H5HFdbg.c',
    'H5HFdblock.c',
    'H5HFdtable.c',
    'H5HFhdr.c',
    'H5HFhuge.c',
    'H5HFiblock.c',
    'H5HFiter.c',
    'H5HFman.c',
    'H5HFsection.c',
    'H5HFspace.c',
    'H5HFstat.c',
    'H5HFtest.c',
    'H5HFtiny.c' ]

H5HG_SOURCES = [
    'H5HG.c',
    'H5HGcache.c',
    'H5HGdbg.c',
    'H5HGquery.c' ]

H5HL_SOURCES = [
    'H5HL.c',
    'H5HLcache.c',
    'H5HLdbg.c',
    'H5HLdblk.c',
    'H5HLint.c',
    'H5HLprfx.c' ]

H5HP_SOURCES = [ 'H5HP.c' ]

H5I_SOURCES = [
    'H5I.c',
    'H5Idbg.c',
    'H5Iint.c',
    'H5Itest.c' ]

H5L_SOURCES = [
    'H5L.c',
    'H5Ldeprec.c',
    'H5Lexternal.c' ]

H5M_SOURCES = [ 'H5M.c' ]

H5MF_SOURCES = [
    'H5MF.c',
    'H5MFaggr.c',
    'H5MFdbg.c',
    'H5MFsection.c' ]

H5MM_SOURCES = [ 'H5MM.c' ]

H5MP_SOURCES = [
  'H5MP.c',
  'H5MPtest.c' ]

H5O_SOURCES = [
    'H5O.c',
    'H5Oainfo.c',
    'H5Oalloc.c',
    'H5Oattr.c',
    'H5Oattribute.c',
    'H5Obogus.c',
    'H5Obtreek.c',
    'H5Ocache.c',
    'H5Ocache_image.c',
    'H5Ochunk.c',
    'H5Ocont.c',
    'H5Ocopy.c',
    'H5Ocopy_ref.c',
    'H5Odbg.c',
    'H5Odeprec.c',
    'H5Odrvinfo.c',
    'H5Odtype.c',
    'H5Oefl.c',
    'H5Ofill.c',
    'H5Oflush.c',
    'H5Ofsinfo.c',
    'H5Oginfo.c',
    'H5Oint.c',
    'H5Olayout.c',
    'H5Olinfo.c',
    'H5Olink.c',
    'H5Omessage.c',
    'H5Omtime.c',
    'H5Oname.c',
    'H5Onull.c',
    'H5Opline.c',
    'H5Orefcount.c',
    'H5Osdspace.c',
    'H5Oshared.c',
    'H5Oshmesg.c',
    'H5Ostab.c',
    'H5Otest.c',
    'H5Ounknown.c' ]

H5P_SOURCES = [
    'H5P.c',
    'H5Pacpl.c',
    'H5Pdapl.c',
    'H5Pdcpl.c',
    'H5Pdeprec.c',
    'H5Pdxpl.c',
    'H5Pencdec.c',
    'H5Pfapl.c',
    'H5Pfcpl.c',
    'H5Pfmpl.c',
    'H5Pgcpl.c',
    'H5Pint.c',
    'H5Plapl.c',
    'H5Plcpl.c',
    'H5Pmapl.c',
    'H5Pmcpl.c',
    'H5Pocpl.c',
    'H5Pocpypl.c',
    'H5Pstrcpl.c',
    'H5Ptest.c' ]

H5PB_SOURCES = [ 'H5PB.c' ]

H5PL_SOURCES = [
   'H5PL.c',
   'H5PLint.c',
   'H5PLpath.c',
   'H5PLplugin_cache.c' ]

H5R_SOURCES = [
    'H5R.c',
    'H5Rdeprec.c',
    'H5Rint.c' ]

H5RS_SOURCES = [
    'H5RS.c'
    ]

H5S_SOURCES = [
    'H5S.c',
    'H5Sall.c',
    'H5Sdbg.c',
    'H5Sdeprec.c',
    'H5Shyper.c',
    'H5Smpio.c',
    'H5Snone.c',
    'H5Spoint.c',
    'H5Sselect.c',
    'H5Stest.c' ]


H5SL_SOURCES = ['H5SL.c',]

H5SM_SOURCES = [
   'H5SM.c',
   'H5SMbtree2.c',
   'H5SMcache.c',
   'H5SMmessage.c',
   'H5SMtest.c' ]


H5T_SOURCES = [
   'H5T.c',
   'H5Tarray.c',
   'H5Tbit.c',
   'H5Tcommit.c',
   'H5Tcompound.c',
   'H5Tconv.c',
   'H5Tcset.c',
   'H5Tdbg.c',
   'H5Tdeprec.c',
   'H5Tenum.c',
   'H5Tfields.c',
   'H5Tfixed.c',
   'H5Tfloat.c',
   'H5Tnative.c',
   'H5Toffset.c',
   'H5Toh.c',
   'H5Topaque.c',
   'H5Torder.c',
   'H5Tpad.c',
   'H5Tprecis.c',
   'H5Tref.c',
   'H5Tstrpad.c',
   'H5Tvisit.c',
   'H5Tvlen.c',
]

H5TS_SOURCES = [ 'H5TS.c' ]

H5UC_SOURCES = [ 'H5UC.c' ]

H5VL_SOURCES = [
   'H5VL.c',
   'H5VLcallback.c',
   'H5VLint.c',
   'H5VLnative.c',
   'H5VLnative_attr.c',
   'H5VLnative_blob.c',
   'H5VLnative_dataset.c',
   'H5VLnative_datatype.c',
   'H5VLnative_file.c',
   'H5VLnative_group.c',
   'H5VLnative_link.c',
   'H5VLnative_introspect.c',
   'H5VLnative_object.c',
   'H5VLnative_token.c',
   'H5VLpassthru.c',
    ]


H5VM_SOURCES = [ 'H5VM.c' ]

H5WB_SOURCES = [ 'H5WB.c' ]

H5Z_SOURCES = [
   'H5Z.c',
   'H5Zdeflate.c',
   'H5Zfletcher32.c',
   'H5Znbit.c',
   'H5Zscaleoffset.c',
   'H5Zshuffle.c',
   'H5Zszip.c',
   'H5Ztrans.c',
    ]

common_SRCS = []
common_SRCS += H5_SOURCES
common_SRCS += H5A_SOURCES
common_SRCS += H5AC_SOURCES
common_SRCS += H5B_SOURCES
common_SRCS += H5B2_SOURCES
common_SRCS += H5C_SOURCES
common_SRCS += H5CS_SOURCES
common_SRCS += H5CX_SOURCES
common_SRCS += H5D_SOURCES
common_SRCS += H5E_SOURCES
common_SRCS += H5EA_SOURCES
common_SRCS += H5F_SOURCES
common_SRCS += H5FA_SOURCES
common_SRCS += H5FD_SOURCES
common_SRCS += H5FL_SOURCES
common_SRCS += H5FO_SOURCES
common_SRCS += H5FS_SOURCES
common_SRCS += H5G_SOURCES
common_SRCS += H5HF_SOURCES
common_SRCS += H5HG_SOURCES
common_SRCS += H5HL_SOURCES
common_SRCS += H5HP_SOURCES
common_SRCS += H5I_SOURCES
common_SRCS += H5L_SOURCES
common_SRCS += H5M_SOURCES
common_SRCS += H5MF_SOURCES
common_SRCS += H5MM_SOURCES
common_SRCS += H5MP_SOURCES
common_SRCS += H5O_SOURCES
common_SRCS += H5P_SOURCES
common_SRCS += H5PB_SOURCES
common_SRCS += H5PL_SOURCES
common_SRCS += H5R_SOURCES
common_SRCS += H5UC_SOURCES
common_SRCS += H5RS_SOURCES
common_SRCS += H5S_SOURCES
common_SRCS += H5SL_SOURCES
common_SRCS += H5SM_SOURCES
common_SRCS += H5T_SOURCES
common_SRCS += H5TS_SOURCES
common_SRCS += H5VL_SOURCES
common_SRCS += H5VM_SOURCES
common_SRCS += H5WB_SOURCES
common_SRCS += H5Z_SOURCES

H5_GENERATED_HEADERS = [
   'H5Edefin.h',
   'H5Einit.h',
   'H5Epubgen.h',
   'H5Eterm.h',
   'H5version.h',
   'H5overflow.h',
    ]


# Generate 'H5detect to generate H5Tinit.c
if cc_can_run
  H5detect = executable('H5detect', 'H5detect.c',
                      c_args: hdf5_default_warnings + defs,
		      dependencies: ext_deps)

  H5Tinit_c = custom_target('H5T_init.c',
	       output: 'H5T_init.c',
	       command: [H5detect, '@OUTPUT@'])

  settings_comp = executable('H5make_libsettings', 'H5make_libsettings.c',
                      c_args: hdf5_default_warnings + defs,
		      dependencies: ext_deps)

  H5lib_settings_c = custom_target('H5lib_settings.c',
	       output: 'H5lib_settings.c',
	       command: [settings_comp, '@OUTPUT@'])

  gen_SRCS = [ H5lib_settings_c, H5Tinit_c ]
else
  H5Tinit_c = configure_file(input: 'H5T_init.cross',
                 output: 'H5T_init.c',
                 copy: true)

  H5lib_settings_c = configure_file(input: 'H5lib_settings.cross',
                 output: 'H5lib_settings.c',
                 copy: true)

  gen_SRCS = [ H5lib_settings_c, H5Tinit_c ]
endif


lib_hdf5_embedded = static_library('su2hdf5',
                          common_SRCS + gen_SRCS,
                          install : false,
                          c_args: hdf5_default_warnings + defs,
                          dependencies : ext_deps)

hdf5_dep = declare_dependency(link_with: lib_hdf5_embedded,
                              dependencies: ext_deps,
                              include_directories: hdf5_include)
