FROM ubuntu:14.04

MAINTAINER Krishna Prasad <praskrishna@gmail.com>

RUN apt-get update
RUN apt-get -y install build-essential g++
RUN apt-get -y install unzip wget make
RUN wget https://github.com/su2code/SU2/archive/master.zip

RUN unzip master.zip
WORKDIR SU2-master
RUN ./configure 
RUN make
RUN make install

ENV SU2_RUN "/usr/local/bin"
ENV SU2_HOME "~/SU2-master"
ENV PATH $PATH:$SU2_RUN
ENV PYTHONPATH $PYTHONPATH:$SU2_RUN

CMD cp /usr/local/bin/SU2_CFD /home/krishna/su2-binaries/.
CMD cp /usr/local/bin/SU2_DEF /home/krishna/su2-binaries/.
CMD cp /usr/local/bin/SU2_DOT /home/krishna/su2-binaries/.
CMD cp /usr/local/bin/SU2_GEO /home/krishna/su2-binaries/.
CMD cp /usr/local/bin/SU2_MSH /home/krishna/su2-binaries/.
CMD cp /usr/local/bin/SU2_SOL /home/krishna/su2-binaries/.

#sudo docker build -t krishna/scfd-su2 .
#sudo docker run -v ~/su2-binaries/.:~/su2-binaries/. -t krishna/scfd-su2
