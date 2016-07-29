#include "cp_options.h"
#include <stdio.h>


Cp_Options::Cp_Options(){
	cpFreq = 0;
	cpPath = new char[256];
	isRestarted = false;	
}


Cp_Options::~Cp_Options(){
	
}

void Cp_Options::setCpFreq(const int cpFreq_){
	cpFreq = cpFreq_;
}

void Cp_Options::setCpPath(const char * cpPath_){
	sprintf(cpPath, "%s", cpPath_);
}

void Cp_Options::setRestartStatus( bool status_){
	isRestarted = status_;
}

int Cp_Options::getCpFreq (){
	return cpFreq;
}

void Cp_Options::getCpPath ( char * cpPath_){
	sprintf(cpPath_, "%s\n", cpPath);
}

bool  Cp_Options::getRestartStatus(){
	return isRestarted;
}


