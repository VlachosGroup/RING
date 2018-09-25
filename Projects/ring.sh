#!/bin/sh

# compile the indicated program
java -jar ../silver/react.bin.jar $@ || exit 1

# produces output.txt
mv output.txt output.cpp

# ensure the cpp files are all built
(cd ../cpp && make) || exit 1

# compile the final program
g++ -g -I ../cpp/ -I../cpp/idas/include output.cpp ../cpp/*.o

