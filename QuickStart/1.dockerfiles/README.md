# Create SU2 binaries in Ubuntu 14.04

* As a first step, create a new folder ~/su2-binaries 
* Install Docker on your machine, in ubuntu
```
sudo apt-get update
sudo apt-get install -y docker
```
* Run docker build command and wait for 45min
```
sudo docker build -t krishna/scfd-su2 .
```
* Open Docker file and go to line 21 thru 26 and change /home/krishna to your <HOME_DIR>
* Run docker run command
```
sudo docker run -v <HOME_DIR>/su2-binaries/.:<HOME_DIR>/su2-binaries/. -t krishna/scfd-su2
```
