## Setup
1. Change the hardcoded full project root path for 'root' in `.\microfluidics-optimization_files\user_files\particle-tracking-repeated-steps-with-bootstrap.scm` (be sure to escape backslashes)

## Optimize
To run the optimization:
1. Open ANSYS Workbench
2. Load the workbench file in the root directory `.\microfluidics-optimization.wb`
3. Open Fluent by double clicking on 'Setup' in the workspace
3. Run the optimization with File > Scripting > Run Script File... > `.\design_search.wbjn`
4. ???
5. Profit!

## One off simulation
To run the one off simulation in Fluent:
1. File > Read > Scheme
2. read in `.\microfluidics-optimization_files\user_files\particle-tracking-repeated-steps-with-bootstrap.scm`

## Updating UDFs
### Prerequisites
- Install Visual Studio 2017 (UDF compilation does not seem to work with VS2022)

### Steps
If you make changes to the udf files at `.\microfluidics-optimization_files\dp0\FFF\Fluent\inertial-lift-udfs.c`
1. ALWAYS delete the old 'libudf' folder before trying to recompile
    - `del /f .\microfluidics-optimization_files\dp0\FFF\Fluent\libudf`
    - `rmdir /f .\microfluidics-optimization_files\dp0\FFF\Fluent\libudf`
2. Make sure the latest Microsoft Visual Studio is installed. During install, make sure to select the add on for 'Desktop C++ development' to install the required compilers
3. Open Fluent via the 'Cross Tools Command Window' (Windows Menu > search 'Cross Tools' > open > run the following commands:)
    - `cd C:\Program Files\ANSYS Inc\ANSYS Student\v232\fluent\ntbin\win64\`
    - `fluent`
4. Load the case file `.\microfluidics-optimization_files\dp0\FFF\Fluent\FFF.3-Setup-Output.cas.h5`
5. Go to the 'User Defined' Tab
6. Click 'Functions' > 'Compiled...'
7. Delete `inertial-lift-udfs.c` from the 'Source Files' list and then click 'Add...' and select `.\microfluidics-optimization_files\dp0\FFF\Fluent\inertial-lift-udfs.c`
8. Build (leaving 'Use Built-In Compiler' unchecked)
    - If there are any errors related to the compilation, try checking the 'Use Built-In Compiler' box and click Build again
9. Click the 'Load' button. You should see a list of the udf functions as the last output if the operation was successful.
10. Attach functions to corresponding hooks in Fluent by clicking 'User Defined' Tab > 'Function Hooks...'

## Debugging in Visual Studio

### Instructions
- Install Visual Studio 2017
    - Install all C and C++ tools for Desktop listed in reference [2]
- Update compiler flags
    - Duplicate the makefile template file (`C:\Program Files\ANSYS Inc\ANSYS Student\v232\fluent\fluent23.2.0\src\udf\makefile_nt.udf`). Rename the old version to `makefile_nt.udf.old`
    - Open `makefile_nt.udf`
    - Change the following lines [1]:
        ```
        CFLAGS = /c /Za /Zi /DUDF_EXPORTING /DUDF_NT /DWIN64 /EHa /wd4224
        link $(LIBS) /dll /debug /assemblydebug /out:$(TARGET)
        ```
    - From eeroi post [1]: "In words, you add the compiler flag /Zi and linker flags /debug /assemblydebug, nothing else. These options generate the program database (.pdb) files that you'll later find in the binary output folders."
    - In ANSYS 2023R2, I needed to add /Za /Zi to the CFLAGS line in makefile_nt.udf
- Open Fluent via the Cross Tools Command Window as above
- In Visual Studio, Open the project folder and open the udf .c file you loaded above
- Attach Visual studio debugger to the Fluent process
    - 'Debug' > 'Attach to Process' > 'Available processes' > 'Filter processes' Search Box
    - Type 'fl' to filter for Fluent associated processes
    - Select 'fl_mpi*.exe'
    - Click Attach
- Set a breakpoint in the function you want to debug (The breakpoint should be a filled red circle.
    - If it is a circle with a red outline, see 'Troubleshooting' > 'No symbol file loaded for libudf.dll' below
- Trigger the relevant function call in Fluent
    - e.g. for a DEFINE_ON_DEMAND UDF macro function, run the function in 'User-Defined' Tab > 'Execute on Demand...'
    - e.g. for a DEFINE_INIT UDF macro function, run 'Outline View' Panel > right click on 'Initialization' > 'Initialize'

### Troubleshooting
Always delete the \libudf folder before recompiling. If the folder cannot be deleted because it is in use, close fluent and stop all 'fl*.exe' and 'cx*.exe' processes with the task manager (This can happen if fluent crashes on a previous run)
### Fluent Compile Error: Fatal error LNK1112
[StackOverflow answer](https://stackoverflow.com/questions/3563756/fatal-error-lnk1112-module-machine-type-x64-conflicts-with-target-machine-typ)
### Visual Studio: The breakpoint will not currently be hit. No symbols have been loaded for this document.
Try: [CFD Online Forum Post - Simplest way debug fluent UDF](https://web.archive.org/web/20211017011134/https://www.cfd-online.com/Forums/fluent-udf/206603-simplest-way-debug-fluent-udf.html)
### Fluent Compile Error: Fatal Error LNK1201
Make sure that Visual Studio is not attached to the fluent process when trying to compile. If so, detach VS from the fluent process and try to compile again.
### Fluent Compile Error: Fatal Error LNK1168
Current solution seems to be to restart Fluent and again?
### Fluent Compile Error: Fatal Error LNK2019, LNK1120 (And many others): 'unresolved external symbol'
This indicates that you have a variable that has not been declared in your code and was not included from an external file or library
This may also be because the files for an external library you have tried to include were not accessible. Include the files in the Compile window (both the .c and .h files). Note that there is a column on the right hand side for header files. [3]

### References
[1] [CFD Online - Simplest way debug fluent UDF, post by user `eeroi`](https://web.archive.org/web/20211017011134/https://www.cfd-online.com/Forums/fluent-udf/206603-simplest-way-debug-fluent-udf.html)

[2] [ANSYS Forum - How can I properly link Visual Studio with ANSYS Fluent, to load in UDFs, post by user `ANSYS_MMadore`](https://web.archive.org/web/20221128211933/https://forum.ansys.com/forums/topic/how-can-i-properly-link-visual-studio-with-ansys-fluent-to-load-in-udfs/)

[3] [ANSYS Forum - UDF compilation problem, post by user `DrAmine`](https://web.archive.org/web/20230127072511/https://forum.ansys.com/forums/topic/udf-compilation-problem/)
