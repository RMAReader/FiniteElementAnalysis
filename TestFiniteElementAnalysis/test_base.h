#ifndef _TEST_BASE_H_
#define _TEST_BASE_H_

#include <vector>
#include <string>


class test_base
{
public:
	test_base();
	~test_base();

	virtual void run();

	std::string name;
	std::vector < std::string > errors;

};


class test_geometry : public test_base
{
public:
	void run();
};


#endif _TEST_BASE_H_