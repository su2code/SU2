from glob import glob
from shutil import copyfile

testfile_names = glob("*_test.cpp")
test_names = [testfile_name[:-9] for testfile_name in testfile_names]

# Write the tests
copyfile("Makefile_header.txt", "Makefile.am")
with open("Makefile.am", 'a') as outfile:
    outfile.write("\n")
    outfile.write("# Programs to be built for testing\n")
    outfile.write("check_PROGRAMS =\n")
    for name in test_names:
        outfile.write("check_PROGRAMS += {0}\n".format(name))
    outfile.write("\n")

    outfile.write("# Test sources\n")
    for name in test_names:
        outfile.write("{0}_SOURCES = {0}_test.cpp\n".format(name))
    outfile.write("\n")

    outfile.write("# Tests to be run\n")
    outfile.write("TESTS =\n")
    for name in test_names:
        outfile.write("TESTS += {0}\n".format(name))
