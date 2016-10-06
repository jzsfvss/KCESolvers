Instructions:
 
1. Download the KCESolvers package from GitHub and unzip it into your "Documents\MATLAB" folder. Download the export_fig package from MathWorks, and unzip it into the "KCESolvers\Open Source" folder.
2. Open MatLab and navigate to the Direct Solver folder through the left sidebar menu called "Current Folder". The path at the top should now point to the selected folder.
3. Type reload at the command line and hit enter. (Good to do this after each run of the solver too.)
4. Run the direct / inverse solver by switching to the appropriate folder, running reload, then typing kced and hitting enter.
5. To enter more sets of parameters, edit the parameters.xls spreadsheet. Note that if you leave t_div and x_div as zero, then the program will prompt you to enter them. The cells must contain some element (zero or else).
6. To mass generate L+C for all columns, run kced and choose "2. All in file" when prompted for "Parameters". All columns must contain parameters for the same KCE method. The highest "t_div" value is used for all columns, the "x_div" value is irrelevant.

--------------------------------------------------------------------------

Notes:

1. You can use the bottom part of the parameters.xls file to enter your parameters in more natural (experimental) units. Open cells are white, generated are grey.
2. Typically use a large time division (~2000), but an x division of 1 (since only the detector is relevant; the solver will work with a large x-mesh though).