////////////////////////////////////////////////////////////////////////////////////////////////////

// File: create_variabe_file.C

// This script writes the data stored in the .csv-files to .root-files (one for each sensor).
// It takes 

//	* 'txt_files/data_input.txt'
//	* 'root_files/column_length.root'
//	* 'data_files/P_coeff_NP02.txt'
//	* 'data_files/T_coeff_NP02.txt'

// as inputs to create a .root-file called 'file_*.root' in the root_files/ directory.
// After starting root, this script can be run by typing: .x create_variable_file.C( int x )
// Hereby x is an integer labeling the x-th sensor that is readout (from first to last columns in .csv-files).

// Author: Lars Bathe-Peters <lars.bathe-peters@cern.ch>
// CERN Summer Student Programme 2019

////////////////////////////////////////////////////////////////////////////////////////////////////


std::string coeff_file_pressure("data_files/P_coeff_NP02.txt");
std::string coeff_file_temperature("data_files/T_coeff_NP02.txt");

int N_SENSORS_TOTAL = 200; // specifies the maximum number of sensors and has to greater than the actual number of sensors in .csv-files 
const int N_SENSORS_CONST_TOTAL = 200; // 
int N_SENSORS; // Number of sensors will be specified later automatically

static double r[N_SENSORS_CONST_TOTAL]; // array for pressure and temperature readout 
int day_length_file; // integer to determine file for setting day_length

double *pressure_coeff;
double *temperature_coeff;

int k = 0;

double* pressure_coefficients( std::string coeff_file, int lo ) {

    static double r[N_SENSORS_CONST_TOTAL];
    std::ifstream inFile;
    inFile.open( coeff_file );

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

	if ( atof( row[1].c_str() ) == 0  ) continue;
        r[k] = atof( row[1].c_str() );
        k = k + 1;
    }
    k = 0;
    inFile.close();
    return r;
}


double* temperature_coefficients( std::string coeff_file, int lo ) {

    std::ifstream inFile;
    inFile.open( coeff_file );

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

        r[k] = atof( row[1].c_str() );
        k = k + 1;
    }

    k = 0;
    inFile.close();
    return r;
}

// Functions to retrieve column_length from 'column_length.root'
double get_column_length( std::string sensor_id_string, std::string date, int lo, std::string insert_zero ) {

    TCanvas *c1 = new TCanvas("c1","c1");

    double column_num;

    TFile tree_file( "root_files/column_count.root", "read" );

    TTree* sensor_tree = nullptr;
    tree_file.GetObject( (std::string("sensors_") + date).c_str(), sensor_tree );
    sensor_tree->Draw( "columns>>hcolumns", ( "sensor_ID == " + sensor_id_string ).c_str(), "");

    TH1* hcolumns = (TH1*)gPad->GetPrimitive("hcolumns");
    column_num = hcolumns->GetMean();

    delete hcolumns;
    delete sensor_tree;
    delete c1;

    tree_file.Close();
    gDirectory->cd( ( std::string( "root_files/file_" ) + insert_zero + std::to_string(lo).c_str() + std::string( ".root:/" ) ).c_str() );
    return column_num;
}


void create_variable_file( int lo ) {

    // Variable declarations
    std::string time_value_string_year;
    std::string time_value_string_month;
    std::string time_value_string_day;
    std::string time_value_string_hour;
    std::string time_value_string_min;
    std::string time_value_string_sec;
    std::string strain_value_string;

    int year;
    int month;
    int day;
    int hour;
    int min;
    int sec;

    int i = 1; // count lines of csv-file
    int j = 0; // count number of values written to .root-file
    int column_count[4000]; // array for column lenghts
    std::fill_n( column_count, 4000, 111 );
    int max = *max_element( std::begin(column_count), std::end(column_count) );
    int breaking_m[64]; // needed to save values of sensor with fewest to largest number of values
    std::fill_n( breaking_m, 64, 1111 );
    bool exists;
    int k = 1;
    int ü = 0;
    int end_line = 0;
    int sensor_initial;
    std::string date;
    std::string sensor_name;
    TDatime time_value;
    double strain_value;
    double ones[9000];
    int column_num;

    const std::string input_filename("txt_files/data_input.txt"); // has to be same .txt file as in 'get_column_length.C'
    std::ifstream input_file( input_filename );
    std::string data_filename;

    pressure_coeff = pressure_coefficients( coeff_file_pressure, lo );
    temperature_coeff = temperature_coefficients( coeff_file_temperature, lo );

    cout << "Pressure " << *( pressure_coeff + lo ) << endl;
    cout << "Temperature " << *( temperature_coeff + lo ) << endl;

    ofstream outFile;
    outFile.open("txt_files/sensors_name_id.txt");

    std::string insert_zero; // for correct name of file

    // Consider digits when recreating .root-files
    if ( lo < 10 ) { // one digit
        insert_zero = "00";
    }
    if ( lo > 9 && lo < 100 ) { // two digits
        insert_zero = "0";
    }
    if ( lo > 99  ) {
        insert_zero = "";
    }

    // Create stream for writing path and filenames into .txt-file
    ofstream root_path_file;

    // Create root file
    TFile* f = new TFile( ( std::string( "root_files/file_" ) + insert_zero + std::to_string(lo).c_str() + std::string( ".root" ) ).c_str(), "RECREATE");

    // Create TTree to that file
    TTree* tree = new TTree("data", "Simple");

    // Create vectors for TTree
    std::vector<double> strainval;
    std::vector<int> timeval;
    std::vector<double> temperatureval;
    std::vector<long> sensorID;
    std::vector<double> displacementval;
    std::vector<double> pressurecoeff;
    std::vector<double> pressureval;
    std::vector<double> temperaturecoeff;
    std::vector<int> startday;
    std::vector<int> daylength;

    // Write branches to that TTree
    tree->Branch("strain_val", &strainval);
    tree->Branch("displacement_val", &displacementval);
    tree->Branch("temperature_val", &temperatureval);
    tree->Branch("pressure_val", &pressureval);
    tree->Branch("time_val", &timeval);
    tree->Branch("sensor_ID", &sensorID);
    tree->Branch("pressure_coefficient", &pressurecoeff);
    tree->Branch("temperature_coefficient", &temperaturecoeff);
    tree->Branch("start_day", &startday);
    tree->Branch("day_length", &daylength);

    int day_length = 0;
    int start_day;

    std::string sensor_id;
    long sensor_ID;

    int count = 0;
    std::string lines;

    // Find number of input csv-files in txt-file
    while ( getline( input_file, lines ) ) count++;

    input_file.close();
    input_file.open( input_filename );

    cout << "Number of input csv-files: " << count << endl;

    // When .csv-file is read, count it with file_num and stop when the last .csv-file has been read
    int file_num = 0;

    while ( input_file ) {

        if ( file_num == count ) break;
	input_file >> data_filename; 
	++file_num;

        date = data_filename.substr(19,16);
        cout << data_filename << endl;
        cout << "Date: [" << date  << "] (from to [dd_mm-dd_mm_year])" << endl;

        ifstream inFile;
 
        inFile.open(data_filename.c_str());

        if (inFile.is_open()) {
            cout << "File '" << data_filename << "' has been opened" << endl;
        }
        else {
            cout << "NO FILE HAS BEEN OPENED" << endl;
        }

        std::vector<std::string> row;
        std::string line, word;

        getline(inFile, line);

        while( inFile ) {
 
            row.clear();
            getline( inFile, line );
            stringstream s(line);

	    if ( i == column_num ) { // readin has reached last value of respective sensor at line i

                int cmax = 0; // counts columns of maximum length in csv-file
                int t = 0; // to change index of array that stores the last column lengths of sensor columns

                cout << "End of file '" << data_filename << "' at " << i << endl;

	        // Reset variables for reading next .csv-file
                i = 1;
                j = 0;
                std::fill_n( column_count, 4000, 111 );
                max = *max_element( std::begin(column_count), std::end(column_count) );
                std::fill_n( breaking_m, 64, 1111 );
                k = 1;
                ü = 0;
                end_line = 0;
                inFile.close();
	        break;
            }

            while (getline(s, word, ',') ) {
                row.push_back(word);
            }

            if (i%10000 == 0) cout << i << endl;
 
            // Sensor names
            if ( i == 1 && file_num == 1 ) {
                sensor_name = row[2*lo + 1];
                sensor_id =  row[2*lo].substr(14, row[2*lo].length() ); // Took not the whole sensor_ID -> Problems in TTree variable, because too long
                std::stringstream iss(sensor_id);
              	iss >> sensor_ID;
//              cout << sensor_ID << endl;
//	        cout << sensor_name << " at " << lo << endl;
//	        sensorID.push_back( sensor_ID[kk-k] );
	        N_SENSORS = row.size()/2;
	        const int N_SENSORS_CONST = row.size()/2;

                for ( int oo = 0; oo < 2 * N_SENSORS; oo++  ) {
	        outFile << row[oo + 1] << "," << row[oo].substr( 14, row[2*lo].length() ) << endl;
	        ++oo;
	        }
	    }
            column_num = get_column_length( sensor_id, date, lo, insert_zero );

            i = i + 1; // next line

            if ( i == 2 ) continue; // skip line after sensor names and IDs
  
            if ( i % (int)ceil(column_num/8640.) != 0  || abs( atof( ( row[2*lo + 1] ).c_str() ) ) > 1e+3 ) continue; // skip values, so that the resulting plot has more or less 8640 values (which corresponds to one value in 10 seconds when datataking with constant time steps in one day)

            // Skip data points according to columns length of sensor in data_input.csv file
            if ( row[2*lo].empty() && end_line == 0 ) {
	        breaking_m[lo] = 2*lo;
//	        cout << column_count[m] << endl;
	        end_line = i;
	        if (column_count[lo] == 111) { 	
	            column_count[lo] = i;
//		    cout << sensor_name[m/2] <<endl;
	        }
	    end_line = 0;
	    }

            sensor_initial = lo;
//          cout << sensor_initial << endl;
//	    cout << "Time: " << row[2*lo] << endl;
//	    cout << "Value: " << row[2*lo+1] << endl;

	    exists = std::find(std::begin( breaking_m ), std::end( breaking_m ), 2*lo) != std::end( breaking_m );

            if ( exists ) { 
                time_value_string_year = "1996";
                time_value_string_month = "01";
                time_value_string_day = "01";
                time_value_string_hour = "01";
                time_value_string_min = "01";
                time_value_string_sec = "01";

                strain_value_string = "999";
	    }

            else {
                time_value_string_year = row[2*lo].substr(0,4);
                time_value_string_month = row[2*lo].substr(5,7);
                time_value_string_day = row[2*lo].substr(8,10);
                time_value_string_hour = row[2*lo].substr(11,13);
                time_value_string_min = row[2*lo].substr(14,17);
                time_value_string_sec = row[2*lo].substr(17,20);
	
                strain_value_string = row[2*lo+1];
	    }

            // Count sensors with most values up to line that is currently readin
            if (max < *max_element( std::begin(column_count), std::end(column_count) ) ) {
                max =  *max_element( std::begin(column_count), std::end(column_count) );
                //cout << "change to " << max << endl;
                ü++; 
	    }
            if ( ü == 2 * N_SENSORS ) break; // break when all sensors have been found

            strain_value = atof(strain_value_string.c_str());

	    if (strain_value != 999 ) {  
                if ( sensor_name.find("TT") != std::string::npos ) {
                    temperatureval.push_back(strain_value);
                    strainval.push_back( 999 );
                    displacementval.push_back( 999 );
                    sensorID.push_back(sensor_ID);
	            pressurecoeff.push_back( *( pressure_coeff + lo ) );
	            temperaturecoeff.push_back( *( temperature_coeff + lo ) );
                    pressureval.push_back( 999 );
                }

                else if ( sensor_name.find("DS") != std::string::npos ) {
                    displacementval.push_back(strain_value);
                    strainval.push_back( 999 );
                    temperatureval.push_back( 999 );
                    sensorID.push_back(sensor_ID);
                    pressurecoeff.push_back( *( pressure_coeff + lo ) );
                    temperaturecoeff.push_back( *( temperature_coeff + lo ) );
                    pressureval.push_back( 999 );
                }

                else if ( sensor_name.find("ID") != std::string::npos ) {
                    displacementval.push_back( 999 );
                    strainval.push_back( 999 );
                    temperatureval.push_back( 999 );
                    sensorID.push_back(sensor_ID);
                    pressurecoeff.push_back( *( pressure_coeff + lo ) );
                    temperaturecoeff.push_back( *( temperature_coeff + lo ) );
                    pressureval.push_back( strain_value );
                }

                else {
                    strainval.push_back(strain_value);
                    displacementval.push_back( 999 );
                    temperatureval.push_back( 999 );
                    sensorID.push_back(sensor_ID);
	            pressurecoeff.push_back( *( pressure_coeff + lo ) );
                    temperaturecoeff.push_back( *( temperature_coeff + lo ) );
                    pressureval.push_back( 999 );
                }
            } // end of if-statement

            // Time value
            year = atof(time_value_string_year.c_str());
            month = atof(time_value_string_month.c_str());
            day = atof(time_value_string_day.c_str());
            hour = atof(time_value_string_hour.c_str());
            min = atof(time_value_string_min.c_str());
            sec = atof(time_value_string_sec.c_str());

            // If there is more that one file, choose the latter one for determining the length of one day
            if ( count == 1 ) day_length_file = 1;
            if ( count > 1 ) day_length_file = 2;

            // Find start day
            if ( j == 1 && file_num == day_length_file ) {
                start_day = day;
	        cout << "Start day " << start_day << endl;
                startday.push_back(start_day);
	    }

            // If day length (number of values in one day per sensor) has not been set yet, day length is set
            if (day == start_day + 1 && day_length == 0 ) {
	        day_length = j;
                daylength.push_back(day_length);
	    }

            time_value.Set( year, month, day, hour, min, sec );

            if (strain_value != 999) {
                 timeval.push_back(time_value.Convert());
	    }

//	    cout << "Start day " << start_day << endl;
//	    cout << "Day length " << day_length << endl;
//	    cout << "File number " << file_num << endl;
            j++; // count line that is readin
        } // while ( inFile ), end of csv-file readout
    } // while ( input_file ), end of loop over csv-files in .txt-file

    // Discard created file, if there is no reasonable readout
    if ( strainval.size() == 0) {
        cout << "Strain value vector size is " << strainval.size() << ". The corresponding file will therefore be discarded."  << endl;
        remove( ( std::string( "root_files/file_" ) + insert_zero + std::to_string(lo).c_str() + std::string( ".root" ) ).c_str() );
    }

    // Write path and name of file to 'root_path_file.txt'
    if ( strainval.size() != 0 ) {
        root_path_file.open("txt_files/root_path_file.txt", ios::app); // appends path_to_root_file to txt-file
        if ( lo < 10 ) { // one digiit
            root_path_file << gSystem->pwd() << "/" << ( std::string( "root_files/file_00" ) + std::to_string(lo).c_str() + std::string( ".root" ) ).c_str() << endl;
        }
        if ( lo > 9 && lo < 100 ) { // two digits
            root_path_file << gSystem->pwd() << "/" << ( std::string( "root_files/file_0" ) + std::to_string(lo).c_str() + std::string( ".root" ) ).c_str() << endl;
        }
        if ( lo > 99  ) {
            root_path_file << gSystem->pwd() << "/" << ( std::string( "root_files/file_" ) + std::to_string(lo).c_str() + std::string( ".root" ) ).c_str() << endl;
        }
        cout << gSystem->pwd() << endl;
        tree->Fill();
        tree->Write("", TObject::kOverwrite);
    }  

    // Save one temperature sensor in an extra .root-file for the temperature correction
    if ( sensor_name.find("TT_ambience_Top") != std::string::npos ) {
        TFile* f2 = new TFile( "root_files/file_TT_ambiance_Top.root", "RECREATE");
        TTree *newtree = tree->CloneTree();
        newtree->Write("", TObject::kOverwrite);
        f2->Close();
        cout << "File 'root_files/file_TT_ambiance_Top.root' has been created for the temperature correction." << endl;
    }
    // Close file
    f->Close();
    return 0;
}
