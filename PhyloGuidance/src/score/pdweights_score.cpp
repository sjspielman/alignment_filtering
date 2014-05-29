// Written by COW. Compile with either of the following lines (depending on your system):
// g++ -std=c++0x -O3 pdweights_score.cpp -o pdweights_score
// clang++ -O3 pdweights_score.cpp -o pdweights_score

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;

void readMSAfromFasta( vector<string>& MSA, const char* filename )
{
	string line, sequence;
	ifstream infile( filename );
	if ( infile.is_open() )
	{
		bool have_sequence = false;
		while ( infile.good() )
		{
			getline( infile, line );
			if ( line.empty() ) continue;
			if ( line[0]=='>' )
			{
				if ( have_sequence )
					MSA.push_back( sequence );
				have_sequence = false;
				sequence = "";
			}
			else
			{
				sequence += line;
				have_sequence = true;
			}
		}
		if ( have_sequence )
			MSA.push_back( sequence );
		infile.close();
	}
	else
	{
		cerr << "Unable to open file " << filename << endl;
		exit( -1 );
	}
	
	
}

void readDistances( vector<vector<double> >& distances, const char* filename )
{
	string line;
	ifstream infile( filename );
	if ( infile.is_open() )
	{
		while ( infile.good() )
		{
			getline( infile, line );
			//cout << line << endl;
			if ( !line.empty() )
			{
				vector<double> row;
				double distance;
				stringstream ls(line);
				while( (ls >> distance) )
					//cout << distance << endl;
					row.push_back( distance ) ;
				distances.push_back( row );
			}
		}
		infile.close();
	}
	else
	{
		cerr << "Unable to open file " << filename << endl;
		exit( -1 );
	}
	
	/*
	for ( unsigned int i=0; i<distances.size(); i++ )
	{
		for ( unsigned int j=0; j<distances[i].size(); j++ )
			//cout << i << " " << j << " " << distances[i][j] << endl;
			cout << " " << distances[i][j];
		cout << endl;
	}
	cout << "==============" << endl;
	*/
}


void readMSA( vector<string>& MSA, const char* filename )
{
	string line;
	ifstream infile( filename );
	if ( infile.is_open() )
	{
		while ( infile.good() )
		{
			getline( infile, line );
			if ( !line.empty() )
				MSA.push_back( line );
		}
		infile.close();
	}
	else
	{
		cerr << "Unable to open file " << filename << endl;
		exit( -1 );
	}
}

void printMSA( const vector<string>& MSA )
{
//	cout << "----------------------------------------------------------" << endl;
	for ( auto itr = MSA.begin(); itr != MSA.end(); ++itr )
	{
		cout << *itr << endl;
	}
//	cout << "----------------------------------------------------------" << endl;
}

template <typename T> void initMatrix( vector<vector<T> >& matrix, const unsigned int x, const unsigned int y, const T init_value = 0 )
{
	matrix.resize( x );
	for ( auto itr=matrix.begin(); itr!=matrix.end(); ++itr )
	{
		(*itr).resize( y );
		for ( auto itr2=(*itr).begin(); itr2!=(*itr).end(); ++itr2 )
			(*itr2) = init_value;
	}
}

void printMatrix( const vector<vector<double> >& matrix, ostream& out=cout )
{
	for ( auto itr=matrix.begin(); itr!=matrix.end(); ++itr )
	{
		auto itr2=(*itr).begin();
		out << *itr2;
		++itr2;
		for ( ; itr2!=(*itr).end(); ++itr2 )
			out << " " << *itr2;
		out << endl;
	}
}

int skipGaps( const string& str, int i )
{
	while ( i < str.size() )
	{
		if ( str[i] != '-' )
			return i;
		else ++i;
	}
	return -1;
}


double calcPosScore( const vector<vector<int> >& map, const vector<vector<double> >& distances, int column, int row, int x0 )
{
	// sum up all the pairs that point to the same location in base MSA
	double sum = 0;
	for ( unsigned int i=0; i<map.size(); i++ )
		if ( map[i][column]==x0 )
			sum += distances[i][row];
	return sum;
}

void calcScores( const vector<string>& base_MSA, const vector<string>& test_MSA, const vector<vector<double> >& distances, vector<vector<double> >& scores )
{
	//This function calculates a matrix of scores for the base MSA, relative to the test MSA. Caution: no sanity checking is performed. If the test MSA does not contain the exact same sequences as the base MSA, the function will fail.


	auto bc = base_MSA[0].size(); // number of columns in base alignment
	auto tc = test_MSA[0].size(); // number of columns in test alignment
	auto n = base_MSA.size(); // number of sequences

	//printMSA( base_MSA );
	//printMSA( test_MSA );
	//cout << bc << " " << tc << " " << n << endl;
	
	// set up the map of positions from test to base array
	// initialize map to -1's
	vector<vector<int> > tb_map;
	initMatrix<int>( tb_map, n, tc, -1 );
	
	//printMatrix( tb_map );
	
	auto bs_itr = base_MSA.begin(); // iterator over base MSA sequences
	auto ts_itr = test_MSA.begin(); // iterator over test MSA sequences
	int j = 0; // counts current sequence
	
	while( bs_itr != base_MSA.end() )
	{
		int bi = 0; // pointer to characters in sequence
		int ti = 0; // pointer to characters in sequence
		
		while( 1 )
		{
			bi = skipGaps( *bs_itr, bi );
			ti = skipGaps( *ts_itr, ti );
			if ( bi < 0 || ti < 0 ) break;
			//cout << bi << " " << (*bs_itr)[bi] << endl;
			//cout << ti << " " << (*ts_itr)[ti] << endl;
			tb_map[j][ti] = bi;
			++bi; ++ti;
		}
		//cout << *bs_itr << endl << *ts_itr << endl << j << endl;
		++bs_itr; ++ts_itr; ++j;
	}
	//printMatrix( tb_map );
	
	
	// now that we have the map, we can calculate the scores in the original MSA
	initMatrix<double>( scores, n, bc, 0 ); // start with matrix of zeros

	for ( int ti=0; ti<tc; ti++ ) // go over all columns in map
	{
		//cout << "Column" << ti << endl;
		for ( j=0; j<n; j++ ) // go over all rows
		{
			//cout << "Row" << j << endl;
			int x0 = tb_map[j][ti];
			
			if ( x0 >= 0 )
				scores[j][x0] += calcPosScore( tb_map, distances, ti, j, x0 );
			//cout << j << endl;
		}
	}
	//printMatrix( scores );
}

int main( int argc, char* argv[] )
{
	if ( argc != 5 )
	{
		cout << "Usage is: " << argv[0] << " <base-MSA file> <test-MSA file> <distance file> <outfile>" << endl;
	}
	else
	{
		vector<string> base_MSA;
		vector<string> test_MSA;
		vector<vector<double> > distances;
		vector<vector<double> > scores;
		
		readMSAfromFasta( base_MSA, argv[1] );
		readMSAfromFasta( test_MSA, argv[2] );
		readDistances( distances, argv[3] );
		calcScores( base_MSA, test_MSA, distances, scores );
		
		ofstream outfile( argv[4] );
		printMatrix( scores, outfile );
	}
}