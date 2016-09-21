#ifndef __CP_OPTIONS_H__
#define __CP_OPTIONS_H__

#include <string>

class Cp_Options{
private:
	int cpFreq;
	int nIter;
	std::string cpPath;
	bool isRestarted;
	bool withSCR;

public:
	Cp_Options();
	~Cp_Options();
		void setCpPath ( const std::string cpPath_); 

	void setCpFreq ( const int cpFreq_);
	int getCpFreq();
	void setnIter(const int nIter_);
	int getnIter();
	std::string getCpPath(); 
	
	void setRestartStatus(bool status);
	bool getRestartStatus();

};

#endif
