#ifndef __CPBASE_HPP__
#define __CPBASE_HPP__
#include <string>

class CpBase
{
protected: 

public:
	CpBase(){}
	~CpBase(){}	
	virtual int read(const std::string* filename) = 0;
	virtual int write(const std::string* filename) = 0;
	virtual int update() = 0;
};

#endif
