////////////////////////////////////////////////////////////////////////////////////////////////////

// File: merge_root_files.C

// Script to merge .root-files in 'txt_files/root_path_file.root'.
// It takes 

//      * 'txt_files/root_path_file.txt'
//	* 'root_files/file_*.root'

// as inputs to create a new .root-file called 'merged_file.root' in the 'root_files/' directory.
// After starting root, this script can be run by typing: .x merge_root_files.C()

// Author: Lars Bathe-Peters <lars.bathe-peters@cern.ch>
// CERN Summer Student Programme 2019

////////////////////////////////////////////////////////////////////////////////////////////////////


void merge_root_files() {

	TChain chain("data");
	TChain chain2("corrections");
	std::ifstream input_root_files("txt_files/root_path_file.txt");

	if (input_root_files.is_open()) {
		cout << "File" << " has been opened" << endl;
	}
	else {
		cout << "NO FILE HAS BEEN OPENED" << endl;
	}

	std::string line;
	
	while ( std::getline( input_root_files, line ) ) {
		chain.Add( line.c_str() );
		chain2.Add( line.c_str() );
		cout << "line " <<  line << endl;
	}

	TFile* f = new TFile( "root_files/merged_file.root", "RECREATE" );

	chain.CloneTree( -1, "fast" );
	chain2.CloneTree( -1, "fast" );
	f->Write();	

	return 0;
}
