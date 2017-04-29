#include "test_base.h"
#include "geometry.h"

void test_geometry::run()
{
	name = "test_geometry";

	geoVEC3F p1(std::array < float, 3 > {{0, 0, 0}});
	geoVEC3F p2(std::array < float, 3 > {{10, 0, 0}});
	geoVEC3F p3(std::array < float, 3 > {{0, 10, 0}});
	
	float r = 5;

	for (int i = -2 * r; i <= 10 + 2 * r; i++)
	{
		float h;
		float x = i;
		float y = -3;
		if (geometry::min_height_sphere_on_triangle(x, y, r, p1, p2, p3, h))
		{
			std::cout << "point (" << x << "," << y << "," << h << ")" << std::endl;
		}
		else
		{
			std::cout << "point (" << x << "," << y << ", - )" << std::endl;
		}
	}

}