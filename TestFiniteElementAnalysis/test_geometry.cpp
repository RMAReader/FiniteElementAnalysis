#include "test_base.h"
#include "geometry.h"

void geometry_min_height_sphere_on_point::run()
{
	name = "geometry_min_height_sphere_on_point";
	float tolerance = 0.000001;

	geoVEC3F p1(std::array < float, 3 > {{0, 0, 0}});
	
	float r = 5;
	float x, y, h;

	y = -3;

	x = 5;
	AssertIsFalse(geometry::min_height_sphere_on_point(x, y, r, p1, h));

	x = 4;
	AssertIsTrue(geometry::min_height_sphere_on_point(x, y, r, p1, h));
	AssertAreEqual(h, -5, tolerance);

	x = 0;
	AssertIsTrue(geometry::min_height_sphere_on_point(x, y, r, p1, h));
	AssertAreEqual(h, -1, tolerance);
}


void geometry_min_height_sphere_on_line::run()
{
	name = "geometry_min_height_sphere_on_line";
	float tolerance = 0.000001;

	geoVEC3F p1(std::array < float, 3 > {{0, 0, 0}});
	geoVEC3F p2(std::array < float, 3 > {{10, 0, 0}});

	float r = 5;
	float x, y, h;

	y = -3;

	x = -0.1;
	AssertIsFalse(geometry::min_height_sphere_on_line(x, y, r, p1, p2, h));

	x = 0;
	AssertIsTrue(geometry::min_height_sphere_on_line(x, y, r, p1, p2, h));
	AssertAreEqual(h, -1, tolerance);
	AssertIsTrue(geometry::min_height_sphere_on_line(x, y, r, p2, p1, h));
	AssertAreEqual(h, -1, tolerance);

	x = 5;
	AssertIsTrue(geometry::min_height_sphere_on_line(x, y, r, p1, p2, h));
	AssertAreEqual(h, -1, tolerance);

	x = 10;
	AssertIsTrue(geometry::min_height_sphere_on_line(x, y, r, p1, p2, h));
	AssertAreEqual(h, -1, tolerance);

	x = 10.1;
	AssertIsFalse(geometry::min_height_sphere_on_line(x, y, r, p1, p2, h));
	
	y = 5.1;
	x = 5;
	AssertIsFalse(geometry::min_height_sphere_on_line(x, y, r, p1, p2, h));
}



void geometry_min_height_sphere_on_triangle::run()
{
	name = "geometry_min_height_sphere_on_triangle";
	float tolerance = 0.000001;

	geoVEC3F p1(std::array < float, 3 > {{0, 0, 0}});
	geoVEC3F p2(std::array < float, 3 > {{10, 0, 0}});
	geoVEC3F p3(std::array < float, 3 > {{0, 10, 0}});

	float r = 5;
	float x, y, h;

	y = 1;

	x = -0.1;
	AssertIsFalse(geometry::min_height_sphere_on_triangle(x, y, r, p1, p2, p3, h));

	x = 0;
	AssertIsTrue(geometry::min_height_sphere_on_triangle(x, y, r, p1, p2, p3, h));
	AssertAreEqual(h, 0, tolerance);
	AssertIsTrue(geometry::min_height_sphere_on_triangle(x, y, r, p2, p3, p1, h));
	AssertAreEqual(h, 0, tolerance);
	AssertIsTrue(geometry::min_height_sphere_on_triangle(x, y, r, p3, p1, p2, h));
	AssertAreEqual(h, 0, tolerance);
	AssertIsTrue(geometry::min_height_sphere_on_triangle(x, y, r, p3, p2, p1, h));
	AssertAreEqual(h, 0, tolerance);
	AssertIsTrue(geometry::min_height_sphere_on_triangle(x, y, r, p2, p1, p3, h));
	AssertAreEqual(h, 0, tolerance);
	AssertIsTrue(geometry::min_height_sphere_on_triangle(x, y, r, p1, p3, p2, h));
	AssertAreEqual(h, 0, tolerance);

	x = 5;
	AssertIsTrue(geometry::min_height_sphere_on_triangle(x, y, r, p1, p2, p3, h));
	AssertAreEqual(h, 0, tolerance);

	x = 9;
	AssertIsTrue(geometry::min_height_sphere_on_triangle(x, y, r, p1, p2, p3, h));
	AssertAreEqual(h, 0, tolerance);

	x = 9.1;
	AssertIsFalse(geometry::min_height_sphere_on_triangle(x, y, r, p1, p2, p3, h));

}