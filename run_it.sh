#//////////////////////////////////////////////////////////////////////////////////////////////////

# File: run_it.sh

# This script contains the terminal prompt commands to generate plots from .csv-files 
# It executes

#       * 'get_column_length.C'
#       * 'create_variable_file.C( int x )'
#       * 'make_variable_corrections.C'
#       * 'merge_root_files.C'
#       * 'make_plots.C( int x, std::string plot_start )' .

# Hereby 'x' is an integer labeling the x-th sensor that is readout (from first to last columns in .csv-files) and
# 'plot_start' is a string determining where the plotting starts. By default, when setting 

#       * 'plot_start = "user_defined"' -> Plotting starts from start of filling time
#       * 'plot_start = "anything else"' -> Plotting starts from beginning of dataset in .csv-files .

# The start of the plotting can be changed by modifying line 101 and uncommenting line 295 (and commenting line 294),
# so that 'start_index' is set to 'plots_index'
# Once ROOT is installed, this script can be run by typing in the terminal prompt: . run_it.sh
# The output is a .root-file called 'plots.root' in the 'root_file/' directory.

# Author: Lars Bathe-Peters <lars.bathe-peters@cern.ch>
# CERN Summer Student Programme 2019

#//////////////////////////////////////////////////////////////////////////////////////////////////


#!/bin/bash
rm 'txt_files/root_path_file.txt'
rm root_files/*.root
   root -l -q -b 'get_column_length.C'
   root -l -q -b 'create_variable_file.C(63)'
   root -l -q -b 'make_variable_corrections.C'
rm 'root_files/file_TT_ambiance_Top.root'
mv root_files/column_count.root .
mv root_files/*.root root_files/file_TT_ambiance_Top.root
printf "File 'root_files/file_TT_ambiance_Top.root' with TTree 'corrections' has been created for the temperature correction.\n"
mv column_count.root root_files/
rm 'txt_files/root_path_file.txt'
for m in {48..48}; do
   root -l -q -b 'create_variable_file.C('$m')'
done
   root -l -q -b 'make_variable_corrections.C'
   root -l -q -b 'merge_root_files.C'
for m in {48..48}; do
   root -l -q -b 'make_plots.C('$m', "user_defined")'
done
