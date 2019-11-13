////////////////////////////////////////////////////////////////////////////////////////////////////

// File: make_plots.C

// This script generates plots for sensors and stores them in a specified format (.pdf, .jpg etc.) 
// and to a .root-file called 'plots.root'.
// It takes 'make_plots_input.txt' as an input where

//      * 'txt_files/sensors_name_id.txt'
//      * 'root_files/merged_file.root'

// is specified as inputs to create a new .root-file called 'plots.root' in the 'root_files/' directory.
// After starting root, this script can be run by typing: .x make_plots.C( int x )
// Hereby x is an integer labeling the x-th sensor that is readout (from first to last columns in .csv-files).

// Author: Lars Bathe-Peters <lars.bathe-peters@cern.ch>
// CERN Summer Student Programme 2019

////////////////////////////////////////////////////////////////////////////////////////////////////


#include <fstream>
#include <sstream>

int ENTRIES = 4000;
int N_SENSORS_TOTAL = 200; // specifies the maximum number of sensors and has to greater than the actual number of sensors in .csv-files, i. e. for 64 sensors N_SENSORS_TOTAL has to be at least 64
int numm;
int i = 0;
int ii = 0;
std::string sensor_name_filename;

// Get sensor names and IDs
void get_sensor_names_and_IDs( std::string *sensor_name, std::string *sensor_ID) {

        std::ifstream inFile;
        inFile.open( "txt_files/sensors_name_id.txt" );

        std::vector<std::string> row;
        std::string line, word;

        while ( inFile ) {
                row.clear();
                getline( inFile, line );
                stringstream s( line );

                if ( inFile.eof() ) break;

                while (getline(s, word, ',') ) {
                        row.push_back(word);
                }

                sensor_name[i] = row[0];
                sensor_ID[i] = row[1];
//              cout << sensor_name[i] << "," << sensor_ID[i] << endl; // Check
                i = i + 1;
        }

        return 0;
}

// Get sensor names
void get_names( std::string* a ) {

        std::ifstream inFile;
        inFile.open( "txt_files/sensors_name_id.txt" );

        std::vector<std::string> row;
        std::string line, word;

        while ( inFile ) {
                row.clear();
                getline( inFile, line );
                stringstream s( line );

                if ( inFile.eof() ) break;

                while (getline(s, word, ',') ) {
                        row.push_back(word);
                }

                a[ii] = row[0];
                ii = ii + 1;
        }

        return 0;
}

// Class for plot modes (strain, displacement, temperature, pressure)
enum class plot_mode{
    strain,
    displacement,
    temperature,
    pressure,
    fplotmode
};

// Struct with variable information
struct VarInfo {
   VarInfo( const std::string& graph_name, int num_entries = ENTRIES ) : fGraphNameBase( graph_name ), fNumEntriess( num_entries ) {}

   std::string fGraphNameBase;
   int fNumEntriess;
};

// Class for generating plots
class plot_maker {
        public:

        plot_maker( plot_mode mode = plot_mode::strain )
                : fPlotMode( mode )
        {
                initialize_variable_info();
        }

        std::pair< std::string, std::vector<TGraph*> > make_plots( const std::string& data_filename, int m, const std::string sensor_name_filename, int plot_color, std::string start ) {

                // Get 'data' and 'corrections' from 'merged_file.root'
                TFile sensor_tree_file( data_filename.c_str(), "UPDATE" );
                fSensorTree = nullptr;
                sensor_tree_file.GetObject( "data", fSensorTree );
                assert( fSensorTree );

                TFile sensor_tree_correction_file( data_filename.c_str(), "UPDATE" );
                fSensorTreeCorr = nullptr;
                sensor_tree_correction_file.GetObject( "corrections", fSensorTreeCorr );
                assert( fSensorTreeCorr );
                fSensorTreeCorr->AddFriend( fSensorTree );

                // Get sensor names and IDs
                std::string sensor_name[N_SENSORS_TOTAL];
                std::string sensor_ID[N_SENSORS_TOTAL];

                get_sensor_names_and_IDs( sensor_name, sensor_ID);

                fSensorName = sensor_name[m];
                fSensorID = sensor_ID[m];

                // Check if TTree with option fSensorID contains values
                numm = fSensorTreeCorr->Draw( "time:strain:strain_0_corr:strain_0avg_corr:strain_0avgp_corr:strain_0avgpt_corr", ( "sensor_ID == " + fSensorID ).c_str() , "goff" );

                if ( numm == 0 ) {
                        cout << "Nothing to readout for this sensor. If this was not expected, check .csv-files. Continue with next sensor." << endl;
                        exit(0);
                }

                // Create vector of type TGraph
                std::vector<TGraph*> graphs;

                cout << "Making plot for sensor " << fSensorName << " (" << fSensorID << ")" << endl;

                graphs = this->build_graphs_from_root( sensor_name, sensor_ID, m, start ); // build graphs using function
                fPlotColor = plot_color;

                return std::pair< std::string, std::vector<TGraph*> >( fSensorName, graphs);
        }

        protected:

        // Declare protected members    
        std::map<std::string, VarInfo> fVariableInfo;
        void initialize_variable_info();
        plot_mode fPlotMode;
        std::string histo_name;
        TTree* fSensorTree;
        TTree* fSensorTreeCorr;
        int fNumEntries;
        int fPlotColor;
        std::string fSensorName;
        std::string fSensorID;

        // Function for graph settings
        void graph_settings( TGraph& temp_graph, std::string graph_name ) {

                temp_graph.SetMarkerStyle(21);
                temp_graph.SetMarkerSize(.2);
                temp_graph.GetXaxis()->SetTimeDisplay(1);
                temp_graph.GetXaxis()->SetNdivisions(503);
                temp_graph.GetXaxis()->SetTimeFormat("%Y-%m-%d %H:%M");
                temp_graph.GetXaxis()->SetTimeOffset(0,"UCT");
                temp_graph.SetTitle( graph_name.c_str() );
        }

        // Function to build TGraphs that are stored in a vector
        std::vector<TGraph*> build_graphs_from_root( std::string sensor_name[], std::string sensor_ID[], int m, std::string start ) {

        std::vector<TGraph*> graphs;

                for (const auto& pair : fVariableInfo) {

                        const std::string& var_name = pair.first;
                        const auto& info = pair.second;

                        std::string graph_name = info.fGraphNameBase;

                        if ( fSensorName.find("DS") != std::string::npos ) {
                                fPlotMode = plot_mode::displacement;
                        }

                        if ( fSensorName.find("TT") != std::string::npos ) {
                                fPlotMode = plot_mode::temperature;
                        }

                        if ( fSensorName.find("ID") != std::string::npos ) {
                               fPlotMode = plot_mode::pressure;
                        }

                        // Plot dataset from user defined start
                        if ( start == "user_defined" ) {
                                numm = fSensorTreeCorr->Draw( "time:strain:strain_0_corr:strain_0avg_corr:strain_0avgp_corr:strain_0avgpt_corr", ( "sensor_ID == " + fSensorID ).c_str() , "goff" );
                        }

                        cout << "Plotting starts from " << start << " value." << endl;

                        if ( start != "user_defined" ) {

                                // Plot from beginning of dataset
                                if ( fPlotMode == plot_mode::strain ) {
                                        numm = fSensorTreeCorr->Draw( "time_val:strain_val:strain_0avgp_corr:strain_0avg_corr:strain_0_corr:strain_0avgpt_corr", ( "sensor_ID == " + sensor_ID[m] ).c_str() , "goff" );
                                }

                                if ( fPlotMode == plot_mode::displacement ) {
                                        numm = fSensorTreeCorr->Draw( "time_val:displacement_val:strain_0_corr:strain_0avg_corr:strain_0avgp_corr:strain_0avgpt_corr", ( "sensor_ID == " + sensor_ID[m] ).c_str(), "goff" );
                                }

                                if ( fPlotMode == plot_mode::temperature ) {
                                        numm = fSensorTreeCorr->Draw( "time_val:temperature_val:strain_0_corr:strain_0avg_corr:strain_0avgp_corr:strain_0avgpt_corr", ( "sensor_ID == " + sensor_ID[m] ).c_str(), "goff" );
                                }

                                if ( fPlotMode == plot_mode::pressure ) {
                                        numm = fSensorTreeCorr->Draw( "time_val:pressure_val2:strain_0_corr:strain_0avg_corr:strain_0avgp_corr:strain_0avgpt_corr", ( "sensor_ID == " + sensor_ID[m] ).c_str(), "goff" );
                                }
                        }

                        cout << "Number of entries: " << numm << endl;

                        // Get values that are to be plotted
                        double *vx = fSensorTreeCorr->GetVal(0);
                        double *vxs = fSensorTreeCorr->GetVal(1);

                        // Create temporary TGraph for data
                        TGraph *temp_graph = new TGraph(numm, vx, vxs );

                        this->graph_settings( *temp_graph, graph_name ); // Apply function to set graph settings
                        graphs.push_back( temp_graph );

                        // There are no further corrections for values of temperature and pressure sensors -> skip these
                        if ( fPlotMode == plot_mode::temperature || fPlotMode == plot_mode::pressure ) continue;

                        double *vxs2 = fSensorTreeCorr->GetVal(4);

                        // Create temporary TGraph for zeroed, moving averaged and pressure corrected data
                        TGraph *temp_graph_pressure = new TGraph(numm, vx, vxs2);
                        this->graph_settings( *temp_graph_pressure, graph_name );
                        graphs.push_back( temp_graph_pressure );

                        double *vxs5 = fSensorTreeCorr->GetVal(5);

                        // Create temporary TGraph for zeroed, moving averaged and pressure and temperature corrected data
                        TGraph *temp_graph_temperature = new TGraph(numm, vx, vxs5);
                        this->graph_settings( *temp_graph_temperature, graph_name );
                        graphs.push_back( temp_graph_temperature );

                }
                return graphs;
        }
}; // End of plot_maker class


// Function for plot settings
void plot_maker::initialize_variable_info() {
        if ( fPlotMode == plot_mode::strain ) {
                fVariableInfo = {

                        { "strain_plot", VarInfo("Strain_vs_Time; Time [y-m-d-h]; Strain [#mum/m]")  },

                        // If you want to plot several variables, uncomment and add the variable information here 
//                      { "strain_presssure_plot", VarInfo("Strain_p_Time; Time [y-m-d-h]; Strain [#mum/m]") }
                };
        }

        else if ( fPlotMode == plot_mode::displacement ) {
                fVariableInfo = {

                        { "displacement_plot", VarInfo("Displacement_vs_Time; Time [y-m-d-h]; Displacement [#mum]") }
                };
        }

        else if ( fPlotMode == plot_mode::temperature ) {
                fVariableInfo = {

                        { "temperature_plot", VarInfo("Temperature_vs_Time; Time [y-m-d-h]; Temperature [#circC]") }
                };
        }

        else {
                fVariableInfo = {

                        { "pressure_plot", VarInfo("Pressure_vs_Time; Time [y-m-d-h]; Pressure [mbar]") }
                };
        }
}

// Main function
void make_plots( int m, std::string start ) {

        // Check if root_path_file.txt exists
        std::ifstream input_root_files("txt_files/root_path_file.txt");

        if (input_root_files.is_open()) {
                cout << "root_path_file.txt exists. Continue with making plots." << endl << endl;
        }
        else {
                cout << "root_path_file.txt does not exist. End of program." << endl;
                exit(0);
        }

        // Input file
        const std::string input_filename("make_plots_input.txt");

        // Determine beginning of plot
        if ( start != "user_defined" ) {
                cout << "Start of plotting is not set. Set start to the beginning of the dataset." << endl;
        }

        double max = 0;
        double min = 0;
        double temp_max = 0;
        double temp_min = 0;

        // Get sensor names
        std::string sensor_name[N_SENSORS_TOTAL];
        get_names(sensor_name);

        // Determine sensor and set corresponding plot mode
        plot_mode fplotmode, strain, displacement, temperature, pressure;

        if ( sensor_name[m].find("DS") != std::string::npos ) {
                fplotmode = plot_mode::displacement;
                cout << "Plot mode of sensor " << sensor_name[m] << " is 'displacement'." << endl;
        }


        else if ( sensor_name[m].find("TT") != std::string::npos ) {
                fplotmode = plot_mode::temperature;
                cout << "Plot mode of sensor " << sensor_name[m] << " is 'temperature'." << endl;
        }

        else if ( sensor_name[m].find("ID") != std::string::npos ) {
                fplotmode = plot_mode::pressure;
                cout << "Plot mode of sensor " << sensor_name[m] << " is 'pressure'." << endl;
        }

        else {
                fplotmode = plot_mode::strain;
                cout << "Plot mode of sensor " << sensor_name[m] << " is 'strain'." << endl;
        }

        plot_maker make_some_plots(fplotmode);

        std::map< std::string, std::vector<TGraph*> > plots_map;

        // Stream for input file
        std::ifstream input_file( input_filename );
        std::string data_filename;

        int plot_color = 1; // Initial plot color

        cout << "Column index of sensor " << sensor_name[m] << " is " << m << "." << endl;

        // Create file
        TFile *f1 = new TFile("root_files/plots.root","UPDATE");

        while ( input_file >> data_filename >> sensor_name_filename ) {
                plots_map.emplace( make_some_plots.make_plots( data_filename, m, sensor_name_filename, plot_color, start ) );
        }

        // Get number of plots
        size_t num_plots = plots_map.cbegin()->second.size();

        cout << "Number of plots: " <<  num_plots << endl;
        cout << "Sensor name: " << plots_map.cbegin()->first.c_str() << endl;

//      plots_map.cbegin()->second.at(0)->Draw("AP"); // Check

        // Create TCanvas
        TCanvas* canvas = new TCanvas;
        canvas->SetName( plots_map.cbegin()->first.c_str() );

        // Create TLegend
        TLegend* legend = new TLegend(0.12, 0.7, 0.35, 0.89);

        // Draw (several TGraphs for one sensor on one cavas) from root file
        for ( size_t p = 0; p < num_plots; ++p ) { // Loop over plots

                if ( p == 0 ) legend->AddEntry( plots_map.cbegin()->second.at(p), plots_map.cbegin()->first.c_str(), "p" );
                if ( p == 1 ) legend->AddEntry( plots_map.cbegin()->second.at(p), (  plots_map.cbegin()->first + "_ZERO_CORR"  ).c_str()  , "p" );
                if ( p == 2 ) legend->AddEntry( plots_map.cbegin()->second.at(p), (  plots_map.cbegin()->first + "_ZERO_AVG_CORR"  ).c_str()  , "p" );
                if ( p == 3 ) legend->AddEntry( plots_map.cbegin()->second.at(p), (  plots_map.cbegin()->first + "_ZERO_AVG_P_CORR"  ).c_str()  , "p" );
                if ( p == 4 ) legend->AddEntry( plots_map.cbegin()->second.at(p), (  plots_map.cbegin()->first + "_ZERO_AVG_PT_CORR"  ).c_str()  , "p" );

                if ( p == 0 ) {
                        plots_map.cbegin()->second.at(0)->Draw("AP");
                }
                else {
                        plots_map.cbegin()->second.at(p)->SetMarkerColor(plot_color);
                        plots_map.cbegin()->second.at(p)->Draw("P");
                }

                ++plot_color;
                if ( plot_color == 3 || plot_color == 5 || plot_color == 6 || plot_color == 7 ) ++plot_color;

        // Find maximum and minimum of plot with TGraphs
        temp_max = TMath::MaxElement( numm, plots_map.cbegin()->second.at(p)->GetY() );
        temp_min = TMath::MinElement( numm, plots_map.cbegin()->second.at(p)->GetY() );

        if ( temp_max > max && temp_max < 5000 ) max = temp_max;
        if ( temp_min < min && temp_min > -70 ) min = temp_min;
        }

        cout << "Maximum is " << max << " and minimum is " << min << "." << endl;

        // Set range of y-axis (so that TLegend also fits on the TGraph), concrete settings can be set after opening a TGraph in 'plots.root'
        plots_map.cbegin()->second.at(0)->GetYaxis()->SetRangeUser( min - 25, max + 70);

        // TLegend settings
        legend->SetTextSize(.035);
        legend->SetBorderSize(0);
        legend->Draw();

        // Save plot as '*.pdf' in 'Plots/' directory
        canvas->SaveAs( ( std::string("Plots/plot_") + ( plots_map.cbegin()->first ).c_str() + std::string(".pdf") ).c_str() );
        canvas->SaveAs( ( std::string("Plots/plot_") + ( plots_map.cbegin()->first ).c_str() + std::string(".jpg") ).c_str() ); // uncommend for '*.jpg'
        f1->cd();
        canvas->Write();

        delete canvas;
        delete legend;
}

