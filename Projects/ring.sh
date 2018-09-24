#!/bin/sh

# compile the indicated program
java -jar ../bin/react.bin.jar $@ || exit 1

# produces output.txt
mv output.txt output.cpp

# ensure the cpp files are all built
(cd ../cpp && make) || exit 1

# compile the final program
g++ -g -I ../cpp/ -I../cpp/idas-1.1.0/include output.cpp ../cpp/*.o

