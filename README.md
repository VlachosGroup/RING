# Rule Input Network Generator (RING)
Rule Input Network Generator (RING) - An automated reaction network generation tool

Instructions to install RING

1) Clone or download the repository on your local drive. RING is distributed in a zip file containing the following folders:

    - silver – contains a jar file of the compiler “react.bin.jar” that translates an input program.
    - cpp – contains C++ implementations of all functions in RING.
    - examples – contains a set of examples written for different kinds of systems.
    - doc – contains some more specific examples with explanations and a documentation of the syntax of
          the language.
    - Projects – contains a linux script and a windows batch file for running RING. This space is also
               where the outputs of RING.
  
  The parent folder contains the license text. RING is available under GNU Lesser GPL 2.1
  
2) Install that Java JRE (http://www.java.com/ ) is installed

3) Download and install gcc (g++ for C++) or Microsoft Visual Studio 2019
   (https://visualstudio.microsoft.com/, choose to download
   just the Visual Studio Community 2019 edition under Visual Studio IDE) depending on the operating system
   
   Additional installation instructions for Windows users:
      - The download will start automatically
      - Double-click on the downloaded installer
      - Click on Yes on the popup followed by clicking Continue
      - Checkmark the second option under Desktop & Mobile (Desktop development with C++)
      - On the right-hand side under Optional checkmark the following

          Optional
          - Just-In-Time debugger
          - VC++ 2017 version 15.7 v14.14 latest v141 tools
          - C++ profiling tools
          - Windows 10 SDK (10.0.17134.0)
          - Visual C++ tools for CMake
          - Visual C++ ATL for x86 and x64
          - Test Adapter for Boost.Test
          - Test Adapter for Google Test
          - Windows 8.1 SDK and UCRT SDK
          - Windows XP support for C++
          - Visual C++ MFC for x86 and x64
          - C++/CLI support
          - Windows 10 SDK (10.0.16299.0) for Desktop C++ [x86 and x64]

      - Wait for the completion
   
4) A good text editor, for example, notepad++ ( http://notepad-plus-plus.org/)
