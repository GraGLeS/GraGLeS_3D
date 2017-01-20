/*
	GraGLeS 2D A grain growth simulation utilizing level set approaches
    Copyright (C) 2015  Christian Miessen, Nikola Velinov

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <exception>
#include <iostream>
#include <sys/time.h>
#include "grainhdl.h"
#include "ggLS.h"
#include "Settings.h"

using namespace std;

int main(int argc, char *argv[]) {

try
{
	if (argc > 1)
		Settings::initializeParameters(argv[1]);
	else
		Settings::initializeParameters();
}
catch(exception& e)
{
	cout<<"Unable to parse parameters file! Details:\n";
	cout<<e.what()<<endl;
	cout<<"Simulation will now halt."<<endl;
	return 0;
}
	grainhdl* my_sim = NULL;

	my_sim = new grainhdl();
try
{
	my_sim->initializeSimulation();
}
catch(exception& e)
{
	cout<<"Unable to initialize simulation! Reason:\n";
	cout<<e.what()<<endl;
	cout<<"Simulation will now halt."<<endl;
	return 0;
}

	if (Settings::MicrostructureGenMode == E_GENERATE_WITH_VORONOY)
		my_sim->saveMicrostructure();

	timeval time_start, time_end;
	gettimeofday(&time_start, NULL);
try
{
	my_sim->run_sim();
}
catch(exception& e)
{
	cout<<"Simulation encountered a fatal error:"<<endl;
	cout<<e.what()<<endl;
	cout<<"Will now stop the simulation."<<endl;
}
	gettimeofday(&time_end, NULL);
	double elapsed_secs = (time_end.tv_sec - time_start.tv_sec)
			+ (time_end.tv_usec - time_start.tv_usec) / 1000000.0;
	cout << "Elapsed seconds for simulation:" << elapsed_secs << endl;

	//my_sim->save_NrGrainsStats();

	my_sim->clear_mem();

	delete my_sim;

}
