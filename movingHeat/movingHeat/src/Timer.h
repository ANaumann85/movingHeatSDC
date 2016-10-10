#ifndef TIMER_H
#define TIMER_H

#include <string>

using namespace std;
struct MyTimer 
{
	//time_t beg, end;
	timespec beg,end;
	std::string name;
	MyTimer(string s):
		name(s)
	{
		//time(&beg);
		clock_gettime(CLOCK_MONOTONIC, &beg);
	}

	~MyTimer() 
	{ 
		//time(&end); 
		clock_gettime(CLOCK_MONOTONIC, &end);
		if( end.tv_nsec < beg.tv_nsec) {
			double tms = (end.tv_sec-beg.tv_sec-1)*1000+(1000000000+end.tv_nsec-beg.tv_nsec)/1000000.0;
			std::cout << name << " needed :" << tms << " ms\n";
		}else {
			double tms = (end.tv_sec-beg.tv_sec)*1000.0+(end.tv_nsec-beg.tv_nsec)/1000000.0;
			std::cout << name << " needed :" <<  tms << " ms\n";
		}
	}

};
#endif
