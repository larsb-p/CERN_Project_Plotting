////////////////////////////////////////////////////////////////////////////////////////////////////

// File: make_variable_corrections.C

// This script applies corrections to the data stored in 'root_files/file_*.root'.
// It takes 

//      * 'root_files/file_*.root'

// as an input to create another TBranch to 'file_*.root' called 'corrections'.
// After starting root, this script can be run by typing: .x make_variable_corrections.C()

// Author: Lars Bathe-Peters <lars.bathe-peters@cern.ch>
// CERN Summer Student Programme 2019

////////////////////////////////////////////////////////////////////////////////////////////////////


int plot_mode= 0; // Integer tagging the plot mode (strain, displacement, temperature, pressure)
int start_index; // Index of plot starting value

int P_REF = 23; // Reference value for pressure correction
int T_REF = 26; // Reference value for temperature correction
double delta_T; // For temperature correction
double offset; // For zeroing
int numm; // Number of entries with given conditions (see below)

int r = 0; // index for loop over values for moving average correction
int rr = 0; // index for loop over values for moving average correction (minus half a day until plus half a day from value with index r)

int temp_time_index; // Index of value that is timewise closest to the time at the start of the filling  
int time_difference = 999; // Calculate time difference to find closest time for every value

int time_difference_filling = 999; // Calculate time difference to find closest time to start of the filling...
int filling_index = 999; // ...and store respective index
int plots_index = 999; // Plotting start index


void make_variable_corrections( ) {

    // Set start time of filling (can be modified freely as long as this date is in between the start and end date of the .csv-files)
    TDatime start_of_filling_time_value;
    start_of_filling_time_value.Set( 2019, 7, 5, 18, 10, 0 );

    // Set time for beginning of plots that is different from the start time of filling (can be modified freely as long as this date
    // is in between the start and end date of the .csv-files)
    TDatime start_of_plots_time_value;
    start_of_plots_time_value.Set( 2019, 7, 7, 23, 03, 57 );

    // Create .root-file for temperature file
    TFile *f_temperature = TFile::Open( ( std::string( "root_files/file_TT_ambiance_Top.root" ) ).c_str(), "READ" );

    TTree* t_temperature = nullptr;

    std::vector<double> *temperature_val = 0;
    std::vector<double> *time_temperature_val = 0;

    TBranch *temperature_branch = 0;
    TBranch *time_temperature_branch = 0;

    int nentries_temperature;

    cout << "Temperature correction file information:" << endl;

    if ( f_temperature->GetListOfKeys()->Contains( "corrections" ) ) {
        cout << "TTree corrections exists in temperature correction file." << endl;
    f_temperature->GetObject( "corrections", t_temperature );

    // Get existing branches
    t_temperature->SetBranchAddress("strain_0avg_corr", &temperature_val, &temperature_branch);
    t_temperature->SetBranchAddress("time", &time_temperature_val, &time_temperature_branch);

    nentries_temperature = t_temperature->Draw( "strain_0avg_corr:time","", "goff" );
    }

    else { // When this script is run for the temperature sensor that is used for the temperature correction (currently set to sensor TT_ambiance_Top at column 63) before it is run for any other sensor every the not averaged corrected data is used for being able to run this script (see first part of run_it.sh)
    f_temperature->GetObject( "data", t_temperature );
        cout << "TTree corrections does not exist in temperature correction file." << endl;
    // Get existing branches
    t_temperature->SetBranchAddress("temperature_val", &temperature_val, &temperature_branch);
    t_temperature->SetBranchAddress("time_val", &time_temperature_val, &time_temperature_branch);

    nentries_temperature = t_temperature->Draw( "temperature_val:time_val","", "goff" );
    }

    cout << "Temperature check: There are " << nentries_temperature << " values and starting value is " << temperature_val->at(0) << "." <<  endl;

    const std::string input_filename("txt_files/root_path_file.txt");
    std::ifstream input_file( input_filename );

    int count = 0;
    std::string lines;

    while ( getline ( input_file, lines ) ) count++;

    input_file.close();
    input_file.open( input_filename );
    cout << "Number of root-files with sensors' readout: " << count << endl << endl;

    cout << "Start of going through files in root_path_file.txt:" << endl << endl;

    std::string data_filename;
    int file_num = 0;
    while ( input_file ) {
        if ( file_num == count ) break;
        input_file >> data_filename;
        ++file_num;

        cout << "Name of .root-file: " << data_filename.substr(data_filename.length() - 24, data_filename.length() ) << endl;

        // Open files with paths and names specified in 'txt_files/path_root_file.root'
        TFile *f = TFile::Open( ( data_filename.substr(data_filename.length() - 24, data_filename.length() ) ).c_str(), "UPDATE" );


        if (f->IsOpen()) {
            cout << "File '" << data_filename << "' has been opened." << endl;
        }
        else {
            cout << "NO FILE HAS BEEN OPENED" << endl;
        }

        // Get TBranch
        TTree* t = nullptr;
        f->GetObject( "data", t );

        std::vector<double> *strain_val = 0;
        std::vector<double> *displacement_val = 0;
        std::vector<double> *pressure_val = 0;
        std::vector<double> *time_val = 0;
        std::vector<double> *pressure_coeff = 0;
        std::vector<double> *temperature_coeff = 0;
        std::vector<int> *start_day = 0;
        std::vector<int> *day_length = 0;

        // Get existing branches
        TBranch *strain_branch = 0;
        t->SetBranchAddress("strain_val", &strain_val, &strain_branch);

        TBranch *pressure_val_branch = 0;
        t->SetBranchAddress("pressure_val", &pressure_val, &pressure_val_branch);

        TBranch *time_branch = 0;
        t->SetBranchAddress("time_val", &time_val, &time_branch);

        TBranch *pressure_coeff_branch = 0;
        t->SetBranchAddress("pressure_coefficient", &pressure_coeff, &pressure_coeff_branch);

        TBranch *temperature_coeff_branch = 0;
        t->SetBranchAddress("temperature_coefficient", &temperature_coeff, &temperature_coeff_branch);

        TBranch *startday_branch = 0;
        t->SetBranchAddress("start_day", &start_day, &startday_branch);

        TBranch *daylength_branch = 0;
        t->SetBranchAddress("day_length", &day_length, &daylength_branch);

        // Create new TTree to 'file_*.root'
        TTree* tree = new TTree("corrections", "Simple");
        f->GetObject("corrections", tree);

        // Make new branches to existing TTree
        std::vector<double> strain;
        std::vector<double> time;
        std::vector<double> strain0corr;
        std::vector<double> strain0avgcorr;
        std::vector<double> strain0avgpcorr;
        std::vector<double> strain0avgptcorr;

        tree->Branch("strain", &strain);
        tree->Branch("time", &time);
        tree->Branch("strain_0_corr", &strain0corr);
        tree->Branch("strain_0avg_corr", &strain0avgcorr);
        tree->Branch("strain_0avgp_corr", &strain0avgpcorr);
        tree->Branch("strain_0avgpt_corr", &strain0avgptcorr);

        double strain_value_ZERO_AVG_CORR;
        double strain_value_ZERO_AVG_P_CORR;

        int nentries4 = t->Draw( "sensor_ID", "", "goff" );
        int nentries = t->Draw( "strain_val:time_val","", "goff" );

        cout << "Entries: " << nentries << endl;

        // Find out what kind of sensor it is (labeled by plot_mode) and write respective value into strain_val (even though it might not be a strain value)
        if ( strain_val->at(0) == 999 ) { // if it is not strain, check displacement
        ++plot_mode;
        t->SetBranchAddress("displacement_val", &strain_val, &strain_branch);
        nentries = t->Draw( "displacement_val:time_val","", "goff" );
                if ( strain_val->at(0) == 999 ) { // if it is not displacement, check temperature
                ++plot_mode;
                cout << "Displacement value: " << strain_val->at(0) << endl;
                t->SetBranchAddress("temperature_val", &strain_val, &strain_branch);
                nentries = t->Draw( "temperature_val:time_val","", "goff" );
                }
                        if ( strain_val->at(0) == 999 ) { // if it is not temperature, check pressure
                        ++plot_mode;
                        cout << "Temperature value: " << strain_val->at(0) << endl;
                        t->SetBranchAddress("pressure_val2", &strain_val, &strain_branch);
                        nentries = t->Draw( "pressure_val2:time_val","", "goff" );
                        }
        }

        if ( plot_mode == 0 ) {
        cout << "This strain sensor's first value is " << strain_val->at(0) << " #mu m/m";
        }

        if ( plot_mode == 1 ) {
        cout << "This displacement sensor's first value is " << strain_val->at(0) << " mm";
        }

        if ( plot_mode == 2 ) {
        cout << "This temperature sensor's first value is " << strain_val->at(0) << " #circ C";
        }

        if ( plot_mode == 3 ) {
        cout << "This pressure sensor's first value is " << strain_val->at(0) << " mbar";
        }

        cout << " at time " << time_val->at(0) << "." << endl;

        int nentries3 = t->Draw( "pressure_coefficient:temperature_coefficient", "", "goff" ); // Get entries
        int nentries2 = t->Draw( "start_day:day_length", "", "goff" );

        cout << "Pressure coefficient: " << pressure_coeff->at(0) << endl;
        cout << "Temperature coefficient: " << temperature_coeff->at(0) << endl;

        // End of data readout from TBrach 'data', beginning of corrections -----------------------------------------------------------------------

        cout << "Start time of filling: ";
        start_of_filling_time_value.Print();
        cout << "Start time of filling in seconds: " << start_of_filling_time_value.Convert() << endl;

        // Find out which time value comes closest to the set filling start time
        for ( int uu = 0; uu < time_val->size(); uu++ ) {
            if ( time_difference > abs( start_of_filling_time_value.Convert() - time_val->at(uu) ) ) {
                time_difference = abs( start_of_filling_time_value.Convert() - time_val->at(uu) );
                filling_index = uu;
            }
        }
        cout << "Difference of closest time to filling time: " << time_difference << " s" <<  endl;
        time_difference = 999; // Reset

        cout << "Start of the filling value has index " << filling_index << " with strain_value " << strain_val->at( filling_index ) << " taken as the offset to zero the data." << endl;


        // Find out which time value comes closest to the set start time of plots
        for ( int uu = 0; uu < time_val->size(); uu++ ) {
            if ( time_difference > abs( start_of_plots_time_value.Convert() - time_val->at(uu) ) ) {
                time_difference = abs( start_of_plots_time_value.Convert() - time_val->at(uu) );
                plots_index = uu;
            }
        }
        cout << "Difference of closest time to set plots starting time: " << time_difference << " s" <<  endl;
        time_difference = 999; // Reset

        cout << "Start of the plots value has index " << plots_index << " with strain_value " << strain_val->at( plots_index ) << " taken as the start value for the plotting." << endl;


        // Set start_index (index of value closest to specified time), this is going to be the first value appearing in plot
        // Set manually to desired starting value for plot (now the plotting starts at the filling time), use option 
        // 'user_defined' in 'strain.C', otherwise the plotting will start from the beginning of the .csv-file
        start_index = filling_index; // Uncomment one
//      start_index = plots_index;

        // Beginning of zeroing -------------------------------------------------------------------------------------------------------------

        // Set offset as the strain value at the starting time of the filling
        offset = strain_val->at(filling_index);
        cout << "Offset: " << offset << endl;

        // Write data into new TBranch 'strain' of TTree 'corrections' to consider starting point of the plotting
        for (int jj = start_index ; jj < nentries; jj++) {
            strain.push_back( strain_val->at(jj) );
            time.push_back( time_val->at(jj) );
            strain0corr.push_back( strain_val->at(jj) - offset );
        }

        // End of zeroing --------------------------------------------------------------------------------------------------------------------

        // Calculate delta_T for temperature correction --------------------------------------------------------------------------------------
        // For temperature correction, put loops like later for strain here to get the delta_T and then implement after pressure correction push_backs
        std::vector<double> delta_T;

        for ( int uuu = 0; uuu < strain_val->size(); uuu++ ) {
            // Find out which delta_T matches the corresponding strain value most timewise
            for ( int uu = 0; uu < time_temperature_val->size(); uu++ ) {
                if ( time_difference > abs( time_val->at( uuu ) - time_temperature_val->at(uu) ) ) {
                    time_difference = abs( time_val->at( uuu ) - time_temperature_val->at(uu) );
                    temp_time_index = uu;
                }
            }
            time_difference = 999; // Reset
            delta_T.push_back( temperature_val->at(temp_time_index) - T_REF );
        }

        cout << "Number of temperature values (should be the same as entries): " << delta_T.size() << endl;
//      cout << "Delta_T size (should be the same as entries): " << delta_T.size() << endl;

        // Beginning of moving average (using convolution), pressure and temperature correction ----------------------------------------------
        // Array filled with ones for convolution
        double ones[day_length->at(0)];
        double ones_val = 1. /day_length->at(0);

        for ( int i = 0; i < day_length->at(0); i++ ) {
            ones[i] = ones_val;
        }

        cout << "Number of values of one day: " << day_length->at(0) << endl << endl;


        // Calculate moving average correction for three different time intervals 
        // 1) [start_index, half a day after start_index)
        // 2) [half a day after start_index, half a day before end date)
        // 3) [half a day before end date, end date]


        // Value loop
//      for ( r = start_index ; r < start_index + day_length->at(0)/2; r++ ) { // 1)
//      for ( r = start_index ; r <= strain_val->size() - day_length->at(0)/2 -1; r++ ) { // 1) and 2)
        for ( r = start_index ; r < strain_val->size(); r++ ) { // 1), 2) and 3)

            strain_value_ZERO_AVG_CORR = ( strain_val->at( r ) - offset ) * ones[rr];

            // Loop over values used for moving average correction (using half a day before to half a day after the current value with index r)
            for ( rr = 1 ; rr <= day_length->at(0)/2; rr++ ) {

            // 1) [start_index, half a day after start_index):
            // Mirror data used for calculating the moving average for the first day_length/2 data points
            if ( r < start_index + day_length->at(0)/2 ) strain_value_ZERO_AVG_CORR += ( strain_val->at( r + rr ) - offset ) * ones[rr]; // 1)
            if ( r < start_index + day_length->at(0)/2 && abs ( r - rr ) < start_index ) strain_value_ZERO_AVG_CORR += ( strain_val->at( r + rr ) - offset ) * ones[rr]; // 1)
            if ( r < start_index + day_length->at(0)/2 && abs ( r - rr ) > start_index ) strain_value_ZERO_AVG_CORR += ( strain_val->at( r - rr ) - offset ) * ones[rr]; // 1)

//          2) [half a day after start_index, half a day before end date)
            if ( r >= start_index + day_length->at(0)/2 && r < strain_val->size() - day_length->at(0)/2 ) {
                strain_value_ZERO_AVG_CORR += ( strain_val->at( r + rr ) - offset ) * ones[rr]; // 2)
                strain_value_ZERO_AVG_CORR += ( strain_val->at( r - rr ) - offset ) * ones[rr]; // 2)
            }

//          3) [half a day before end date, end date]
            if ( r >= ( strain_val->size() - day_length->at(0)/2 ) ) strain_value_ZERO_AVG_CORR += ( strain_val->at( r - rr ) - offset ) * ones[rr]; // 3)
            if ( r >= ( strain_val->size() - day_length->at(0)/2 ) && abs ( r + rr ) < strain_val->size() ) strain_value_ZERO_AVG_CORR += ( strain_val->at( r + rr ) - offset ) * ones[rr]; // 3)
            if ( r >= ( strain_val->size() - day_length->at(0)/2 ) && abs ( r + rr ) > strain_val->size() ) strain_value_ZERO_AVG_CORR += ( strain_val->at( r - rr ) - offset ) * ones[rr]; // 3)

            } // End of moving average loop

            // Fill vectors
        if ( plot_mode == 0 || plot_mode == 1 ) { // for strain and displacement sensors
            strain0avgcorr.push_back( strain_value_ZERO_AVG_CORR ); // Moving average correction
            strain0avgpcorr.push_back( strain_value_ZERO_AVG_CORR + pressure_coeff->at(0) * P_REF ); // Pressure correction with moving average smoothening before
            strain_value_ZERO_AVG_P_CORR = strain_value_ZERO_AVG_CORR + pressure_coeff->at(0) * P_REF;
            strain0avgptcorr.push_back( strain_value_ZERO_AVG_P_CORR + delta_T.at(r) * temperature_coeff->at(0) ); // Temperature correction on top of that
        }

        if ( plot_mode == 2 || plot_mode == 3 ) { // for temperature and pressure sensors do not apply pressure and temperature corrections
            strain0avgcorr.push_back( strain_value_ZERO_AVG_CORR + offset );
            strain0avgpcorr.push_back( 999 );
            strain0avgptcorr.push_back( 999 );
        }

//      cout << "strain0avg " << strain_value_ZERO_AVG_CORR << " at " << r << "." << endl; // Check
        strain_value_ZERO_AVG_CORR = 0.; // Reset
        } // End of value loop

        r = filling_index; // Reset loop indices
        rr = 0;

        plot_mode = 0; // Reset plot mode

        // End of moving average, pressure and temperature correction -----------------------------------------------------

        // Fill TTree 'corrections'
        tree->Fill();
        tree->Write("", TObject::kOverwrite);

        // Close files
        f->Close();
        f_temperature->Close();
    } // while ( input_file ), end of root-file readout from root_path_file.txt

    return 0;
}
