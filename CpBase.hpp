#ifndef CPBASE_HPP
#define CPBASE_HPP

class CpBase
{
protected: 

public:
	//CpBase(){}
	//~CpBase(){}	
	virtual void read(const std::string* filename) = 0;
	virtual void write(const std::string* filename) = 0;
	virtual void update() = 0;
};

#endif
