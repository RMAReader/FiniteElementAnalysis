#include "test_base.h"


test_base::test_base()
{
}


test_base::~test_base()
{
}

void test_base::run()
{
}


void test_base::AssertIsTrue(bool value)
{
	try{
		if (value == false){ throw std::string("Expected True but was False"); }
	}
	catch (std::string e)
	{
		errors.push_back(e);
	}
}
void test_base::AssertIsFalse(bool value)
{
	try{
		if (value == true){ throw std::string("Expected False but was True"); }
	}
	catch (std::string e)
	{
		errors.push_back(e);
	}
}
void test_base::AssertAreEqual(float actual, float expected, float tolerance)
{
	try
	{
		if (abs(actual - expected) > abs(tolerance))
		{
			std::stringstream err;
			err << "Error: actual = " << actual;
			err << ", expected = " << expected;
			err << ", tolerance = " << tolerance;
			err << std::endl;
			throw err.str();
		}
	}
	catch (std::string e)
	{
		errors.push_back(e);
	}
}
