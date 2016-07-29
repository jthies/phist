#ifndef __CP_OPTIONS_H__
#define __CP_OPTIONS_H__

#include <string>

class Cp_Options{
private:
	int cpFreq;
	char * cpPath;
	bool isRestarted;
	bool withSCR;
public:
	Cp_Options();
	~Cp_Options();
		void setCpPath ( const char * cpPath_); 

	void setCpFreq ( const int cpFreq_);
	int getCpFreq();
	void getCpPath ( char * cpPath_); 
	
	void setRestartStatus(bool status);
	bool getRestartStatus();



};

#endif
