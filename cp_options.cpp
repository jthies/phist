#include "cp_options.h"
#include <stdio.h>


Cp_Options::Cp_Options(){
	cpFreq = 0;
	nIter = 0;
	cpPath = new char[256];
	isRestarted = false;	
}


Cp_Options::~Cp_Options(){
	
}

void Cp_Options::setCpFreq(const int cpFreq_){
	cpFreq = cpFreq_;
}

void Cp_Options::setnIter(const int nIter_){
	nIter = nIter_;
}

void Cp_Options::setCpPath(const std::string cpPath_){
	cpPath = cpPath_;
}

void Cp_Options::setRestartStatus( bool status_){
	isRestarted = status_;
}

int Cp_Options::getCpFreq (){
	return cpFreq;
}

std::string Cp_Options::getCpPath (){
	return cpPath;
}

bool  Cp_Options::getRestartStatus(){
	return isRestarted;
}

int Cp_Options::getnIter(){
	return nIter;
}

