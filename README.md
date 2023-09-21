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
- Install Visual Studio 2017
    - [Install all C and C++ tools for Desktop listed in this ANSYS forum post](https://web.archive.org/web/20221128211933/https://forum.ansys.com/forums/topic/how-can-i-properly-link-visual-studio-with-ansys-fluent-to-load-in-udfs/) and as shown below
![imagevs2017](https://github.com/williamstoy/inertial-microfluidics-optimization/assets/51916749/1bc1b992-ba58-4a1b-81fd-ae954140eff1)

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

## Troubleshooting
Always delete the \libudf folder before recompiling. If the folder cannot be deleted because it is in use, close fluent and stop all 'fl*.exe' and 'cx*.exe' processes with the task manager (This can happen if fluent crashes on a previous run)
### Fatal error LNK1112:
[StackOverflow answer](https://stackoverflow.com/questions/3563756/fatal-error-lnk1112-module-machine-type-x64-conflicts-with-target-machine-typ)
### The breakpoint will not currently be hit. No symbols have been loaded for this document.
Try: [CFD Online Forum Post - Simplest way debug fluent UDF](https://web.archive.org/web/20211017011134/https://www.cfd-online.com/Forums/fluent-udf/206603-simplest-way-debug-fluent-udf.html)

From post by eeroi (eero) on CFD Online Forum:
```
- Remember to back original files up before you make any changes.
- Change the makefile template makefile_nt.udf (in one of the Fluent installation folders):
    CFLAGS = /c /Za /Zi /DUDF_EXPORTING /DUDF_NT /DWIN64 /EHa /wd4224
    link $(LIBS) /dll /debug /assemblydebug /out:$(TARGET)

In words, you add the compiler flag /Zi and linker flags /debug /assemblydebug, nothing else. These options generate the program database (.pdb) files that you'll later find in the binary output folders.

- Compile the UDF code in Fluent GUI and load the libudf library. This seems to be the best (most stable and universal) approach for compilation these days, earlier you had to run nmake in command prompt.
- Open the source file in the "libudf/src" folder in Visual Studio (a free Express version will do fine).
- In Visual Studio, attach the code to fl_mpi1950.exe (Fluent 2019R3). This works in serial mode, in parallel you may have to attach to all Fluent processes (haven't used in years, don't remember anymore). You should see debug symbols loaded in VS.
- Try adding a breakpoint on a code row. The red circle on the left hand side should be filled, if the attach was successful.
- When you run the simulation, whenever the code with the breakpoint is called within Fluent, you should now break in Visual Studio and be able to see the variable content and track program flow. Please see https://ibb.co/K2pPh6K for a screenshot.
```
