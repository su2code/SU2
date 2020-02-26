FROM su2code/build-su2:20191029

# Copies your code file from your action repository to the filesystem path `/` of the container
COPY runTests.sh /runTests.sh

# Code file to execute when the docker container starts up (`entrypoint.sh`)
ENTRYPOINT ["/runTests.sh"]
