////////////////////////////////////////////////////////////////////////////////////////////////////

// File: get_column_length.C

// This script finds the sensor dependent and hence variable column length (i. e. the number of values readout by each sensor).
// The names of the .csv-files are specified in 

//      * 'txt_files/data_input.txt'

// as the input for this script.

// It will automatically find the number of x- and y-values, find the sensor ID in order to tag each sensor and write this information to a .root-file called 'column_length.root' in the root_files/ directory.
// After starting root, this script can be run by typing: .x get_column_length.C()

// Author: Lars Bathe-Peters <lars.bathe-peters@cern.ch>
// CERN Summer Student Programme 2019

////////////////////////////////////////////////////////////////////////////////////////////////////


const int N_SENSORS_TOTAL = 200; // With our 64 sensors, N_SENSORS_TOTAL = 64 would be enough (this has to be always greater than the number of sensors)
int N_SENSORS; // = 64; // Variable of how many sensors shall be plotted (will be automatically set to number of y-coordinate columns in csv-file)
const int BOGUS = 999;

int i = 1; // count lines of csv-file
int u = 0; // column index
int uu = 0; // count how many sensors' column lengths have been found 
int r; // index for loop over sensors with maximum column length

void get_column_length( ) {

    // To determine whether length of column (corresponding to sensor) has not been counted yet
    int column_count[4000];
    std::fill_n( column_count, 4000, 111 ); // fill it with 111 to use it later to find out whether it has been filled with value or not

    // Label for column
    int column_index[N_SENSORS_TOTAL];
    std::fill_n( column_index, N_SENSORS_TOTAL, BOGUS );

    // Labels sensors with longest column
    int new_number[N_SENSORS_TOTAL];

    bool exists;
    int k = 1;
    int end_line = 0;
    std::string sensor_name[N_SENSORS_TOTAL];

//    const std::string input_filename("data_input_short.txt"); // Be sure not to have additional blank lines in txt-file (will be counted)
    const std::string input_filename("txt_files/data_input.txt");
    std::ifstream input_file( input_filename );
    std::string data_filename;

    int count = 0;
    std::string lines;

    // Find number of input csv-files in txt-file
    while ( getline( input_file, lines ) ) count++;

    input_file.close();
    input_file.open( input_filename );
    cout << "Number of input csv-files: " << count << endl;

    std::vector<std::string> row;
    std::string line, word;

    // Create root file
    TFile* f = new TFile("root_files/column_count.root", "RECREATE");

    // Input stream to readout csv-file
    ifstream inFile;

    int file_num = 0; // count input csv-files (data_filename-files without input_file)
    while ( input_file ) {
        if ( file_num == count ) break;
    input_file >> data_filename;
    ++file_num;
    cout << "File number " << file_num << " is called: " << data_filename << endl;

    std::string date = data_filename.substr(19,16);
    cout << "Date: [" << date  << "] (from to [dd_mm-dd_mm-yyyy])" << endl;

    // Create TTree to that file
    TTree* tree = new TTree( ( std::string("sensors_") + date ).c_str(), "CREATE");

    // Write branches to that TTree
    std::vector<long> sensorID;
    std::vector<int> columns;

    tree->Branch("sensor_ID", &sensorID);
    tree->Branch("columns", &columns);

    std::string sensor_id[N_SENSORS_TOTAL];
    long sensor_ID[N_SENSORS_TOTAL];

    // Open file
    inFile.open(data_filename.c_str());
    if (inFile.is_open()) {
        cout << "File '" << data_filename << "' has been opened" << endl;
    }
    else {
        cout << "NO FILE HAS BEEN OPENED" << endl;
    }

    cout << "--------------------------SENSORS------------------------------------------------" << endl;

    getline(inFile, line);

    // Start of reading file 
    while(inFile) {
        row.clear();
        getline( inFile, line );
        stringstream s(line);


        // If end of file is reached... (sensor columns of maximum length have not been readout yet)
        if ( inFile.eof() )  {
            int cmax = 0; // counts columns of maximum length in csv-file
            int t = 0; // to change index of array that stores the last column lengths of sensor columns

            cout << "End of file '" << data_filename << "' at " << i + 1 << endl;

            // Loop over sensors 
            for (u = 0; u < N_SENSORS; u++) {
                exists = std::find(std::begin(column_index), std::end(column_index), u) != std::end(column_index); // find the columns that have not been readout yet
                    if ( column_index[u] == 999 ) {
                        cmax = cmax + 1;
                    }

                    if ( !exists ) {
                        new_number[t] = u; // store column index in array
                        t = t + 1;
                    }
            }

            // Loop over number of maximum columns
            for ( r = 0; r < cmax; r++) {
                column_index[N_SENSORS - cmax + r] = new_number[r]; // Fill remaining array elements
                cout << sensor_name[ new_number[r] ] << " (" << sensor_ID[ new_number[r] ] << "):";
                cout << " change to " << column_index[N_SENSORS - cmax + r] << " at " << i << ", sensors' columns found is " << N_SENSORS - cmax + r + 1 <<  endl;

                columns.push_back(i);
               sensorID.push_back( sensor_ID[ new_number[r] ] );
            }

        cout << "--------------------------SENSORS------------------------------------------------" << endl;
        cout << endl;

        // End of first file readout. Reset variables etc.
        cmax = 0;
        t = 0;
        uu = 0;
        i = 1;
        std::fill_n( column_count, 4000, 111 ); // fill it with 111 to use it later to find out whether it has been filled with value or not
        std::fill_n( column_index, 64, 999 );
        std::fill_n( new_number, 64, 999 );
        k = 1;
        end_line = 0;
        inFile.close();
        break;
        }

        while (getline(s, word, ',') ) {
            row.push_back(word);
        }

        if (i%10000 == 0) cout << "Reading out line: " << i << endl;

        // Get sensor names and IDs in first line of file
        if (i == 1) {
            cout << "Number of sensors in .csv-file is " << row.size()/2 << endl;
            N_SENSORS = row.size()/2;
            for ( int kk = 1; kk < N_SENSORS_TOTAL + 2; kk++ ) {
                if ( kk  > row.size() ) break; // stop when reaching end of file
                sensor_name[kk-k] = row[kk];
                sensor_id[kk-k] =  row[kk - 1].substr(14, row[kk-1].length() ); // Took not the whole sensor_ID -> Problems in TTree variable, because too long
                std::stringstream iss(sensor_id[kk-k]);
                iss >> sensor_ID[kk-k];
//              cout << sensor_ID[kk-k] << endl;
//              cout << sensor_name[kk-k] << endl;
                kk = kk + 1;
                k = k + 1;
            }
        }

        i = i + 1; // Counter for line that is currently read in csv-file

        // Loop over sensors and go through every line of csv-file
        for ( int m = 0; m <= 2 * N_SENSORS; m++ ) {
            if ( row[m].empty() && end_line == 0 ) { // if empty cell...
                end_line = i;
                if (column_count[m/2] == 111) { // if column has not been counted yet (tagged by 111)...
                    column_count[m/2] = i; // ...store respective line
                  column_index[uu] = m/2;     // ...store column index

                    cout << sensor_name[m/2] << " (" << sensor_ID[m/2] << "):";
                    cout << " change to sensor " << column_index[uu] << " at " << i -1 << ", sensors' columns found is " << uu +1 <<  "." << endl;

                    columns.push_back(i-1);
                    sensorID.push_back(sensor_ID[ m/2 ]);

                    uu = uu + 1;
                }
                end_line = 0;
            }
            m = m + 1;
        }
    }

    // Fill and write TTree
    tree->Fill();
    tree->Write("", TObject::kOverwrite);
    delete tree;
    }

f->Close();
return 0;
}

