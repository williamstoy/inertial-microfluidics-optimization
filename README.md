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
- Install Visual Studio 2019 (UDF compilation does not seem to work with VS2022)

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
- Install Visual Studio 2019
    - Install all C and C++ tools for Desktop
- Open Fluent via the Cross Tools Command Window as above
- In Visual Studio 2019, Open the project folder and open the udf .c file you loaded above
- Attach Visual studio debugger to the Fluent process
    - 'Debug' > 'Attach to Process' > 'Available processes' > 'Filter processes' Search Box
    - Type 'fl' to filter for Fluent associated processes
    - Select 'fl_mpi*.exe'
    - Click Attach
- Set a breakpoint in the function you want to debug (The breakpoint should be a filled red circle.
    - If it is a circle with a red outline, see 'Troubleshooting' > 'No symbol file loaded for libudf.dll' below

## Troubleshooting
Always delete the \libudf folder before recompiling. If the folder cannot be deleted because it is in use, close fluent and stop all 'fl*.exe' and 'cx*.exe' processes with the task manager (This can happen if fluent crashes on a previous run)
### Fatal error LNK1112:
[StackOverflow answer](https://stackoverflow.com/questions/3563756/fatal-error-lnk1112-module-machine-type-x64-conflicts-with-target-machine-typ)
### No symbol file loaded for libudf.dll
Not sure how to solve this yet
