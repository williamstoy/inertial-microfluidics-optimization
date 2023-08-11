To run the one off simulation:
1. file > read > scheme
2. read in `.\microfluidics-optimization_files\user_files\particle-tracking-repeated-steps-with-bootstrap.scm`

To run the optimization:
1. Open ANSYS Workbench
2. Load the workbench file in the root directory `.\microfluidics-optimization.wb`
3. Change the parameters as necessary
4. Run the optimization by clicking 'Update All Design Points'
5. ???
6. Profit!

If you make changes to the udf files at `.\microfluidics-optimization_files\dp0\FFF\Fluent\inertial-lift-udfs.c`
1. Make sure the latest Microsoft Visual Studio is installed. During install, make sure to select the add on for 'Desktop C++ development' to install the required compilers
2. Open Fluent via the 'Cross Tools Command Window' (Windows Menu > search 'Cross Tools' > open > run the following commands:)
    - `cd C:\Program Files\ANSYS Inc\ANSYS Student\v232\fluent\ntbin\win64\`
    - `fluent`
3. Go to the 'User Defined' Tab
4. Click 'Functions' > 'Compiled...'
5. Delete and re-add the `inertial-lift-udfs.c` from the 'Source Files' list (`.\microfluidics-optimization_files\dp0\FFF\Fluent\inertial-lift-udfs.c`)
6. Build (leaving 'Use Built-In Compiler' unchecked)
7. If there are any errors related to the compilation, save the case data and re-open Fluent using the following instructions