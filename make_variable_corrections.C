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

double get_column_length( std::string sensor_id_string, std::string date ) {

    double column_num;

    TFile tree_file( "root_files/column_count.root", "read" );
    TTree* sensor_tree = nullptr;
    tree_file.GetObject( (std::string("sensors_") + date).c_str(), sensor_tree );

    sensor_tree->Draw( "columns>>hcolumns", ( "sensor_ID == " + sensor_id_string ).c_str(), "");
//    cout << sensor_id_string << endl;
    TH1* hcolumns = (TH1*)gPad->GetPrimitive("hcolumns");
    column_num = hcolumns->GetMean();
//    cout << column_num << " column_num " << endl;
//    cout << hcolumns->GetMean() << endl;;

    delete hcolumns;
    delete sensor_tree;
 
    tree_file.Close();
    gDirectory->cd("root_files/merged_file.root:/");
 
   return column_num;
}


	int i = 0;
void get_sensor_names( std::string *sensor_name, std::string *sensor_ID) {
        
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
//        cout << sensor_name[i] << "," << sensor_ID[i] << endl;
        i = i + 1;
        }
    return 0;
}

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
    f_temperature->GetObject( "data", t_temperature );

    std::vector<double> *temperature_val = 0;
    std::vector<double> *time_temperature_val = 0;

    // Get existing branches
    TBranch *temperature_branch = 0;
    t_temperature->SetBranchAddress("temperature_val", &temperature_val, &temperature_branch);

    TBranch *time_temperature_branch = 0;
    t_temperature->SetBranchAddress("time_val", &time_temperature_val, &time_temperature_branch);
	
    int nentries_temperature = t_temperature->Draw( "temperature_val:time_val","", "goff" );
    cout << "Temperature check " << nentries_temperature << " and " << temperature_val->at(0) <<  endl;

    const std::string input_filename("txt_files/root_path_file.txt");
    std::ifstream input_file( input_filename );

    int count = 0;
    std::string lines;

    while ( getline ( input_file, lines ) ) count++;

    input_file.close();
    input_file.open( input_filename );
    cout << "Number of root-files with sensors' readout: " << count << endl;

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
            cout << "File '" << data_filename << "' has been opened" << endl;
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
		cout << "Displacement value " << strain_val->at(0) << endl;
	        t->SetBranchAddress("temperature_val", &strain_val, &strain_branch);
        	nentries = t->Draw( "temperature_val:time_val","", "goff" );
		}
                	if ( strain_val->at(0) == 999 ) { // if it is not temperature, check pressure
			++plot_mode;
                	cout << "Temperature value " << strain_val->at(0) << endl;
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
//	start_index = plots_index;

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

//        cout << "Delta_T size (should be the same as entries): " << delta_T.size() << endl;


        // Beginning of moving average, pressure and temperature correction (using convolution) ----------------------------------------------

        // Array filled with ones for convolution
        double ones[day_length->at(0)];
        double ones_val = 1. /day_length->at(0);

        for ( int i = 0; i < day_length->at(0); i++ ) {
            ones[i] = ones_val;
        }

	cout << "Number of values of one day: " << day_length->at(0) << endl;


        // Calculate moving average correction for three different time intervals 
	// 1) [start_index, half a day after start_index]
	// 2) [half a day after start_index, half a day before end date]
	// 3) [half a day before end date, end date]


        // 1) [start_index, half a day after start_index]:
	// Mirror data used for calculating the moving average for the first day_length/2 data points

	// Loop over values
        for ( r = start_index ; r < start_index + day_length->at(0)/2; r++ ) {
            if ( r == start_index ) {
	        if ( plot_mode == 0 || plot_mode == 1 ) {
		    strain0avgcorr.push_back( 0. );
		    strain0avgpcorr.push_back( strain0avgcorr.at(0) + pressure_coeff->at(0) * P_REF );
                    strain0avgptcorr.push_back( strain0avgpcorr.at(0) + delta_T.at(r) * temperature_coeff->at(0) );
		}
                if ( plot_mode == 2 || plot_mode == 3 ) strain0avgcorr.push_back( 999 );
                continue; // first value is zeroed anyways
            }  

	    // Loop over values used for moving average correction (using half a day before to half a day after the current value with index r)
            for ( rr = 0; rr < day_length->at(0); rr++ ) {
                strain_value_ZERO_AVG_CORR += ( strain_val->at( abs( r - day_length->at(0)/2 + rr + 1 ) ) - offset ) * ones[rr]; // Convolution is the moving average
//		cout << "Index " << abs( r - day_length->at(0)/2 + rr + 1 ) << " and " << strain_value_ZERO_AVG_CORR << endl; // Check
	    }

	    // Discard non-phyical values and set these to 100
	    if ( abs( strain_value_ZERO_AVG_CORR) > 1e+8 ) {
	        strain_value_ZERO_AVG_CORR = 100;
	    }

	    // Fill vectors
	    if ( plot_mode == 0 || plot_mode == 1 ) { // for strain and displacement sensors
	        strain0avgcorr.push_back( strain_value_ZERO_AVG_CORR ); // Moving average correction
		strain0avgpcorr.push_back( strain_value_ZERO_AVG_CORR + pressure_coeff->at(0) * P_REF ); // Pressure correction with moving average smoothening before
                strain_value_ZERO_AVG_P_CORR = strain_value_ZERO_AVG_CORR + pressure_coeff->at(0) * P_REF;
                strain0avgptcorr.push_back( strain_value_ZERO_AVG_P_CORR + delta_T.at(r) * temperature_coeff->at(0) ); // Temperature correction on top of that
	    }

            if ( plot_mode == 2 || plot_mode == 3 ) { // for temperature and pressure sensors do not apply these corrections
	        strain0avgcorr.push_back( 999 );
		strain0avgpcorr.push_back( 999 );
                strain0avgptcorr.push_back( 999 );
	    }

//	    cout << "strain0avg " << strain_value_ZERO_AVG_CORR << " at " << r << " and strain value at " << strain_val->at(  abs( r - day_length->at(0)/2 + 0 + 1 ) ) << endl; // Check
            strain_value_ZERO_AVG_CORR = 0.;
	} // End of 1) loop

	r = 0; // Reset loop indices
	rr = 0;


//	2) [half a day after start_index, half a day before end date]:
	// Use day_length/2 data points before and after data point that is moving average corrected	
        for ( r = start_index + day_length->at(0)/2; r < strain_val->size() - day_length->at(0)/2 -1; r++ ) {
	    for ( rr = 0; rr < day_length->at(0); rr++ ) {
	        strain_value_ZERO_AVG_CORR += ( strain_val->at( abs( r - day_length->at(0)/2 + rr + 1 ) ) - offset ) * ones[rr];
//		cout << "Index " << abs( r - day_length->at(0)/2 + rr + 1 ) << " and " << strain_value_ZERO_AVG_CORR << endl; // Check
	    }

            if ( abs( strain_value_ZERO_AVG_CORR) > 1e+8 ) {
                strain_value_ZERO_AVG_CORR = 100;
            }

            if ( plot_mode == 0 || plot_mode == 1 ) {
                strain0avgcorr.push_back( strain_value_ZERO_AVG_CORR ); // Moving average correction
                strain0avgpcorr.push_back( strain_value_ZERO_AVG_CORR + pressure_coeff->at(0) * P_REF ); // Pressure correction with moving average smoothening before
		strain_value_ZERO_AVG_P_CORR = strain_value_ZERO_AVG_CORR + pressure_coeff->at(0) * P_REF;
                strain0avgptcorr.push_back( strain_value_ZERO_AVG_P_CORR + delta_T.at(r) * temperature_coeff->at(0) ); // Temperature correction on top of that
	    }

	    if ( plot_mode == 2 || plot_mode == 3 ) {
                strain0avgcorr.push_back( 999 );
                strain0avgpcorr.push_back( 999 );
                strain0avgptcorr.push_back( 999 );
            }

//	    cout << "strain0avg2 " << strain_value_ZERO_AVG_CORR << endl; // Check
	    strain_value_ZERO_AVG_CORR = 0.;
	} // End of 2) loop

        r = 0; // Reset loop indices
        rr = 0;


//      3) [half a day before end date, end date]:
        // Mirror data used for calculating the moving average for the last day_length/2 data points
        for ( r = strain_val->size() - day_length->at(0)/2 - 1; r < strain_val->size(); r++ ) {
	    for ( rr = 0; rr < day_length->at(0); rr++ ) {
	        if ( abs( r - day_length->at(0)/2 + rr + 1 ) > strain_val->size() - 1 ) {
                    strain_value_ZERO_AVG_CORR += ( strain_val->at( abs( r - day_length->at(0) + rr + 1 ) ) - offset ) * ones[rr];		
		}
		else {
                    strain_value_ZERO_AVG_CORR += ( strain_val->at( abs( r - day_length->at(0)/2 + rr + 1 ) ) - offset ) * ones[rr]; // Convolution is the moving average
		}
//		cout << "Index " << abs( r - day_length->at(0)/2 + rr + 1 ) << " and " << strain_value_ZERO_AVG_CORR << endl; // Check
	    }

            if ( abs( strain_value_ZERO_AVG_CORR) > 1e+8 ) {
                strain_value_ZERO_AVG_CORR = 100;
            }

            if ( plot_mode == 0 || plot_mode == 1 ) {
                strain0avgcorr.push_back( strain_value_ZERO_AVG_CORR ); // Moving average correction
                strain0avgpcorr.push_back( strain_value_ZERO_AVG_CORR + pressure_coeff->at(0) * P_REF ); // Pressure correction with moving average smoothening before
                strain_value_ZERO_AVG_P_CORR = strain_value_ZERO_AVG_CORR + pressure_coeff->at(0) * P_REF;
                strain0avgptcorr.push_back( strain_value_ZERO_AVG_P_CORR + delta_T.at(r) * temperature_coeff->at(0) ); // Temperature correction on top of that
            }

            if ( plot_mode == 2 || plot_mode == 3 ) {
                strain0avgcorr.push_back( 999 );
                strain0avgpcorr.push_back( 999 );
                strain0avgptcorr.push_back( 999 );
            }
//        cout << "strain0avg3 " << strain_value_ZERO_AVG_CORR << endl; // Check
        strain_value_ZERO_AVG_CORR = 0.;
        }

        r = filling_index; // Reset loop indices
        rr = 0;

	plot_mode = 0; // Reset plot mode

        // End of moving average, pressure and temperature correction (using convolution) -----------------------------------------------------

	// Fill TTree 'corrections'
	tree->Fill();
	tree->Write("", TObject::kOverwrite);

	// Close files
	f->Close();
	f_temperature->Close();
    } // while ( input_file ), end of root-file readout from root_path_file.txt

    return 0;
}
