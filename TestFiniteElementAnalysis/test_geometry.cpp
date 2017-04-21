#include "test_base.h"
#include "geometry.h"

void test_geometry::run()
{
	name = "test_geometry";

	geoVEC3F p1(std::array < float, 3 > {{0, 0, 0}});
	geoVEC3F p2(std::array < float, 3 > {{0, 2, 0}});
	geoVEC3F p3(std::array < float, 3 > {{2, 1, 2}});
	
	float cx = 0;
	float cy = 1;
	float r = 1;

	float h1, h2, h3, h4, h5, h6;

	bool cz1 = geometry::min_height_sphere_on_triangle(cx, cy, r, p1, p2, p3, h1);
	bool cz2 = geometry::min_height_sphere_on_triangle(cx, cy, r, p2, p3, p1, h2);
	bool cz3 = geometry::min_height_sphere_on_triangle(cx, cy, r, p3, p1, p2, h3);
	bool cz4 = geometry::min_height_sphere_on_triangle(cx, cy, r, p1, p3, p2, h4);
	bool cz5 = geometry::min_height_sphere_on_triangle(cx, cy, r, p3, p2, p1, h5);
	bool cz6 = geometry::min_height_sphere_on_triangle(cx, cy, r, p2, p1, p3, h6);


	cx = 0;
	cy = 0;
	cz1 = geometry::min_height_sphere_on_triangle(cx, cy, r, p1, p2, p3, h1);


	p1 = geoVEC3F(std::array < float, 3 > {{0, 0, 0}});
	p2 = geoVEC3F(std::array < float, 3 > {{2, 0, 2}});
	
	cx = 0;
	cy = 0;
	r = 1;

	geometry::mesh_edge<float> e;
	e.initialise(p1, p2);
	cz3 = e.min_height_sphere(cx, cy, r, h3);

	e.initialise(p2, p1);
	cz4 = e.min_height_sphere(cx, cy, r, h4);

}