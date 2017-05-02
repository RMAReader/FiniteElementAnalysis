#ifndef _TEST_BASE_H_
#define _TEST_BASE_H_

#include <vector>
#include <string>
#include <sstream>

class test_base
{
public:
	test_base();
	~test_base();

	virtual void run();

	std::string name;
	std::vector < std::string > errors;

	void AssertIsTrue(bool);
	void AssertIsFalse(bool);
	void AssertAreEqual(float,float,float);

};


class geometry_min_height_sphere_on_point : public test_base { public:void run(); };
class geometry_min_height_sphere_on_line : public test_base { public:void run(); };
class geometry_min_height_sphere_on_triangle : public test_base { public:void run(); };


#endif _TEST_BASE_H_