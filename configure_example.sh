# /usr/bin.bash

./configure -prefix=$SU2_RUN/.. --with-MPI=mpicxx --with-Metis-lib=$SU2_HOME/metis-4.0.3/ --with-Metis-include=$SU2_HOME/metis-4.0.3/Lib/

# ./configure --with-MPI=mpicxx --with-Metis-lib=/usr/local/metis-4.0.3 --with-Metis-include=/usr/local/metis-4.0.3/Lib --with-Tecio-lib=/usr/local/tecio-2013 --with-Tecio-include=/usr/local/tecio-2013/tecsrc --with-CGNS-lib=/usr/local/cgnslib_3.1.4/lib --with-CGNS-include=/usr/local/cgnslib_3.1.4/include --prefix=${SU2_HOME}

#./configure -prefix=$SU2_RUN/.. --with-Metis-lib=$SU2_HOME/metis-4.0.3/ --with-Metis-include=$SU2_HOME/metis-4.0.3/Lib/ --with-CGNS-lib=/usr/local/lib --with-CGNS-include=/usr/local/include

# ./configure -prefix=$SU2_RUN/.. --with-MPI=mpicxx --with-Metis-lib=$SU2_HOME/metis-4.0.3/ --with-Metis-include=$SU2_HOME/metis-4.0.3/Lib/ --disable-CFD --disable-GPC --disable-MDC --disable-SMC --disable-PBC --disable-MAC