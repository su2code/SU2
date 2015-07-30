# Create SU2 binaries in various platforms like Ubuntu 14.04, CentOS

* Install Docker on your machine, in ubuntu
```
sudo apt-get update
sudo apt-get install -y docker
```
* Run docker build command and wait for 45min
```
sudo docker build -t krishna/scfd-su2 .
```
* Create a folder called su2-binaries under /tmp folder
* Run docker run command
```
sudo docker run -v /tmp/su2-binaries/.:/tmp/su2-binaries/. -t krishna/scfd-su2
```

Go to /tmp/su2-binaries and you will have all the SU2 binaries, copy to a location, add them to .bashrc PATH variable and start using them.
