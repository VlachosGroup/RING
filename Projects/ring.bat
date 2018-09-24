@ECHO OFF

java -jar ../silver/react.bin.jar %*

IF %ERRORLEVEL% EQU 0 (

 copy /Y output.txt output.cpp

 call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat"

 pushd ..\cpp
 nmake /F Makefile.nmake
 cd ../examples
 popd
 
 cl /EHsc /MD output.cpp ..\cpp\*.obj /I ..\cpp /I ..\cpp\idas\include /Feprogram.exe

 IF %ERRORLEVEL% EQU 0 (
  echo.
  echo Build Succeeded.
 ) ELSE (
  echo.
  echo Build failed.
 )

) ELSE (
 echo.
 echo Build failed.
)
