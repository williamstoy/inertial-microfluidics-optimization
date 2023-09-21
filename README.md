Setup
1. Change the hardcoded full project root path for 'root' in `.\microfluidics-optimization_files\user_files\particle-tracking-repeated-steps-with-bootstrap.scm` (be sure to escape backslashes)

To run the optimization:
1. Open ANSYS Workbench
2. Load the workbench file in the root directory `.\microfluidics-optimization.wb`
3. Open Fluent by double clicking on 'Setup' in the workspace
3. Run the optimization with File > Scripting > Run Script File... > `.\design_search.wbjn`
4. ???
5. Profit!

To run the one off simulation in Fluent:
1. File > Read > Scheme
2. read in `.\microfluidics-optimization_files\user_files\particle-tracking-repeated-steps-with-bootstrap.scm`

If you make changes to the udf files at `.\microfluidics-optimization_files\dp0\FFF\Fluent\inertial-lift-udfs.c`
0. ALWAYS delete the old 'ibudf' folder before trying to recompile
    - `del /f .\microfluidics-optimization_files\dp0\FFF\Fluent\libudf`
    - `rmdir /f .\microfluidics-optimization_files\dp0\FFF\Fluent\libudf`
1. Make sure the latest Microsoft Visual Studio is installed. During install, make sure to select the add on for 'Desktop C++ development' to install the required compilers
2. Open Fluent via the 'Cross Tools Command Window' (Windows Menu > search 'Cross Tools' > open > run the following commands:)
    - `cd C:\Program Files\ANSYS Inc\ANSYS Student\v232\fluent\ntbin\win64\`
    - `fluent`
3. Load the case file `.\microfluidics-optimization_files\dp0\FFF\Fluent\FFF.3-Setup-Output.cas.h5`
4. Go to the 'User Defined' Tab
5. Click 'Functions' > 'Compiled...'
6. Delete `inertial-lift-udfs.c` from the 'Source Files' list and then click 'Add...' and select `.\microfluidics-optimization_files\dp0\FFF\Fluent\inertial-lift-udfs.c`
7. Build (leaving 'Use Built-In Compiler' unchecked)
   - If there are any errors related to the compilation, try checking the 'Use Built-In Compiler' box and click Build again
8. Click the 'Load' button. You should see a list of the udf functions as the last output if the operation was successful.

Troubleshooting:
Always delete the \libudf folder before recompiling. If the folder cannot be deleted because it is in use, close fluent and stop all 'fl*.exe' processes
