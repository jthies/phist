#ifndef __CP_OPTIONS_H__
#define __CP_OPTIONS_H__


typedef struct {

	int cpFreq;
	char * cpPath;
	bool isRestarted;
	bool withSCR;
//	cpOpt->cpFreq = cpFreq_;
/*	CP_Options();
	~CP_Options();
	void setCpPath ( const char * cpPath_); 

	void getCpFreq( int * cpFreq_);
	void getCpPath ( char * cpPath_); 
	
	void setRestartStatus(bool status);
	bool getRestartStatus();
*/
}CP_Options;
// TODO: build the functions for cpFreq and cpPath 
//
//void setCpFreq( struct CP_Options * cpOpt, const int cpFreq_);


#endif
