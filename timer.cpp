#include <ctime>
//#include <chrono>
#include <cstdlib>
#include <iostream>

using namespace std;

int main(int argc, char** argv)
{
	/* get time */
	time_t startTime = time(0);

	/* run command */
	system(argv[1]);

	/* get time */
	time_t endTime = time(0);

	/* calculate elapsed time */
	time_t elapsedTime = endTime - startTime;

	/* print data */
	cout << "The program began: " << asctime(localtime(&startTime)) << endl;

	cout << "The program ended: " << asctime(localtime(&endTime)) << endl;

	cout << "Elapsed time: " << elapsedTime << " seconds." << endl;

	cout << "Done." << endl;

    return 0;
}
