FROM centos:7

MAINTAINER Krishna Prasad <praskrishna@gmail.com>

RUN yum -y update; yum clean all
RUN yum -y install gcc-c++
RUN yum -y install unzip wget make
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

CMD cp /usr/local/bin/* /tmp/su2-binaries/.

#sudo docker build -t krishna/scfd-su2 .
#sudo docker run -v /tmp/su2-binaries/.:/tmp/su2-binaries/. -t krishna/scfd-su2
