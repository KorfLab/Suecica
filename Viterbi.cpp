//Viterbi.cpp : Runs an HMM on integer observations given a parameter file
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

class Viterbi
{
private:
	string param_filename;
	string obs_filename;
	bool control_use;
	vector<int> obs;							// Create vector to contain HMM observations
	vector<char> states;						// Create vector to contain HMM states
	vector<double> start_p;						// Create vector to contain HMM start probabilities
	vector<vector<double> > trans_p;			// Create vector to contain HMM transition probabilities
	vector<vector<double> > emit_p;				// Create vector to contain HMM emission probabilities
	vector<vector<int> > path;					// Create vector to contain HMM paths. Each nested vector contains a path for each starting state
	vector<vector<double> > V;					// Create vector to contain HMM path probabilities. Each entry corresponds to the cumulative probability for each starting state's path

	void ParseParams();
	bool is_long(const string &str);

public:
	string GetOptimalPath();
	Viterbi(string t_paramfilename, string t_obsfilename);
};

string Viterbi::GetOptimalPath()
{
	// Now that parameters have been stored in the object, calculate the most probable path and return it
	const int emit_num = emit_p[0].size();	// Number of possible outputs from each state (with Prob != 0)
	const int obs_num = obs.size();
	const int state_num = states.size();
	path.resize(state_num);
	for (int i = 0; i < state_num; i++) path[i].reserve(obs_num);
	V.resize(state_num);
	for (int i = 0; i < state_num; i++) V[i].resize(2);
	
	// Start by initializing Viterbi paths for t=0 for each state

	for (int i = 0; i < state_num; i++)
	{
		if (obs[0] < emit_num) {
			V[i][0] = start_p[i] + emit_p[i][obs[0]];
		}
		else {V[i][0] = 1;} // For when probability is so low for one of the states, floating point error says it is 0. This leads to log(0) errors, so just say log(prob)=1 aka prob=0 for such observations.
		path[i].push_back(i);
	}

	// Continue Viterbi for remaining observations

	vector <double> temp_prob;
	temp_prob.resize(state_num);
	for (int t = 1; t < obs_num; t++) // Cycle through each observation
	{
		for (int y = 0; y < state_num; y++) // Find the probability of seeing an observation for each state
		{
			bool non_zero_value = true; // For when probability is so low it's 0 (or left out due to a user p-value cutoff)
			for (int j = 0; j < state_num; j++)
			{
				if (obs[t] > emit_num) {non_zero_value = false;} // Probability is too low, so we will keep each path's new state the same as the previous state
			}

			if (non_zero_value == true)	// Calculate path normally
			{
				for (int y0 = 0; y0 < state_num; y0++)	{temp_prob[y0] = V[y0][0] + trans_p[y0][y] + emit_p[y][obs[t]];}
			}
			else // Can't calculate log(0) so keep previous path position's state
			{
				for (int y0 = 0; y0 < state_num; y0++)	{temp_prob[y0] = V[y0][0];}
			}

			double max_prob = temp_prob[0];
			int max_state = 0;
			for (int y0 = 1; y0 < state_num; y0++)
			{
				if (temp_prob[y0] > max_prob)
				{
					max_prob = temp_prob[y0];
					max_state = y0;
				}
			}
			V[y][1] = temp_prob[max_state];
			path[y].push_back(max_state);
		}
		for (int y0 = 0; y0 < state_num; y0++)
			V[y0][0] = V[y0][1];
	}

	// Determine path of maximum likelihood and return it

	double max_prob = V[0][0];		// V[0][0] = Final prob for first state, V[1][0] = Final prob for second state, etc.
	int max_state = 0;
	for (int i = 1; i < state_num; i++)
	{
		if (V[i][0] > max_prob)
		{
			max_prob = V[i][0];
			max_state = i;
		}
	}
	string finalpath = "";
	for (int i = 0; i < obs_num; i++)
	{
		if (!path[max_state][i] == max_state)
			max_state = path[max_state][i];
		finalpath += states[path[max_state][i]];
	}
	return finalpath;
}

bool Viterbi::is_long(const string &str)
{
	return str.find_first_not_of("0123456789.e-") == string::npos;
}

void Viterbi::ParseParams()
{
	ifstream inFile;
	const int NULL_VAL = 0, STATES = 1, START_PROB = 2, TRANS_PROB = 3, EMIT_PROB = 4, OBSERVATIONS = 5;
	int mode = NULL_VAL;
	inFile.open(param_filename.c_str());
	if (!inFile.is_open())
	{
		cout << "Error opening " << param_filename << ", program terminating" << endl;
		exit(EXIT_FAILURE);
	}
	string line;
	int cur_state = -1, prev_state = -1, next_state = 0; // cur_state used in emit_p, prev_state and next_state used in trans_p
	while (inFile.good())
	{
		inFile >> line;
		if (line == "States:") {mode = STATES;}
		else if (line == "Start") {inFile >> line; mode = START_PROB;}
		else if (line == "Transition")	// Prepare incoming transition probabilities to enter a nested vector
		{
			inFile >> line;
			mode = TRANS_PROB;
			trans_p.resize(states.size());
			for (unsigned int i = 0; i < states.size(); i++)
				trans_p[i].resize(states.size());
		}
		else if (line == "Emission")	// Prepare incoming emission probabilities to enter a vector
		{
			inFile >> line;
			mode = EMIT_PROB;
			emit_p.resize(states.size());
		}
		else if (mode == STATES)
		{
			char temp_state = line[0];
			states.push_back(temp_state);
		}
		else if (mode == START_PROB)
		{
			char state;
			inFile >> state;	// State
			inFile >> line;		// Start probability
			double start_prob = log(strtod(line.c_str(), NULL));
			start_p.push_back(start_prob);
		}
		else if (mode == TRANS_PROB)
		{
			//	S	S	0.99998
			//		D	1.0000000000010001e-05
			//		T	1.0000000000010001e-05
			inFile >> line;			// Next state or Transition probability
			if (!is_long(line))		// Line was next state, so get transition probability, as well as increment prev_state and reset next_state
			{
				prev_state++;
				next_state = 0;
				inFile >> line;		// Transition probability being grabbed
			}
			else {next_state++;}	// Line was transition probability
			double trans_prob = log(strtod(line.c_str(), NULL));
			trans_p[prev_state][next_state] = trans_prob;
		}
		else if (mode == EMIT_PROB)
		{
			//	S	0	0.555019350543
			//		1	0.326768919177
			string line2;
			int tab_pos = 0;
			getline(inFile, line2);
			for (unsigned int i = 1; i < line2.length(); i++) {
				if (line2.substr(i,1) == "\t") {
					tab_pos = i;
					break;
				}
			}
			if (tab_pos != 0) {
				// New state's integer observation and observation probability are present on line
				cur_state++;
				line2 = line2.substr(tab_pos+1, line2.length() - (tab_pos + 1));
			}
			double emit_prob = log(strtod(line2.c_str(), NULL));
			emit_p[cur_state].push_back(emit_prob); // (Integer observations assumed to be ordered such that they start from 0 and increment by 1 in the parameter file)
		}
	}
	inFile.close();
	vector<char>(states).swap(states);
	vector<double>(start_p).swap(start_p);
	// Read in observations
	inFile.open(obs_filename.c_str());
	inFile.seekg(0,ios::end);			// Used to find end of file, determine file size
	int obs_num = inFile.tellg();		// Maximum number of possible observations present is equal to length of observation string, minus 1 = obs_num - current - 1
	obs.reserve(obs_num);				// Reserve the HMM observation vector to have enough elements to hold the maximum possible number of observations
	inFile.clear();						// Clear EOF flag
	inFile.seekg(0,ios::beg);			// Return to start of observations
	int t_obs;
	while (inFile >> t_obs) {obs.push_back(t_obs);}	// Putting each observation into an array element
	if (inFile.eof()) {}
	else if (inFile.fail())
	{
		cout << "An error occurred while reading " << param_filename << ": data type mismatch." << endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		cout << "An unkown error occurred while reading " << param_filename << endl;
		exit(EXIT_FAILURE);
	}
}

Viterbi::Viterbi(string t_paramfilename, string t_obsfilename)
{
	param_filename = t_paramfilename;
	obs_filename = t_obsfilename;
	ParseParams();
}

int main(int argc, char* argv[])
{
	string param_filename = "";
	string obs_filename = "";
	bool control_use = false;
	const int NULL_VAL = 0, GRAB_FILENAME = 1, GRAB_OBS = 2;
	int mode = NULL_VAL;
	for (int i = 0; i < argc; i++)
	{
		string temp = argv[i];
		if (temp == "-p") {mode = GRAB_FILENAME;}
		else if (temp == "-o") {mode = GRAB_OBS;}
		else if (mode == GRAB_FILENAME)
		{
			param_filename += argv[i];
		}
		else if (mode == GRAB_OBS)
		{
			obs_filename += argv[i];
		}
	}
	Viterbi myHMM(param_filename, obs_filename);
	string path = myHMM.GetOptimalPath();
	cout << path;
	
	return 0;
}