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


void geometry_mesh3f_range::run()
{
	name = "geometry_mesh3f_range";
	float tolerance = 0.000001;

	geometry::mesh3f mesh;
	geoVEC3F p1(std::array < float, 3 > {{0, 0, 0}});
	geoVEC3F p2(std::array < float, 3 > {{1, 0, 0}});
	geoVEC3F p3(std::array < float, 3 > {{2, 0, 0}});
	geoVEC3F p4(std::array < float, 3 > {{0, 1, 0}});
	geoVEC3F p5(std::array < float, 3 > {{1, 1, 0}});
	geoVEC3F p6(std::array < float, 3 > {{2, 1, 0}});
	geometry::triangle<float, 3> t1(p1, p2, p4);
	geometry::triangle<float, 3> t2(p4, p5, p2);
	geometry::triangle<float, 3> t3(p2, p3, p5);
	geometry::triangle<float, 3> t4(p5, p6, p3);
	mesh.add_triangle(t1);
	mesh.add_triangle(t2);
	mesh.add_triangle(t3);
	mesh.add_triangle(t4);

	int nx = 1;
	int ny = 1;
	geometry::mesh3f_region region(&mesh, nx, ny);
	
	region.set_region(-10, 10, -10, 10);
	region.set_region(-10,-5, -10, -5);
	region.set_region(5, 10, 5, 10);
}
