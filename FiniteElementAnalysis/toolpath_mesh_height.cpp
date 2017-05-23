#include "toolpath_base.h";


toolpath_mesh_height::toolpath_mesh_height(geometry::mesh3f* mesh, int nx, int ny, float r, float default_z)
{
	this->region = geometry::mesh3f_region(mesh, nx, ny);
	this->r = r;
	this->default_z = default_z;
}


float toolpath_mesh_height::get_height(float x, float y)
{
	
	region.set_region(x-r, x+r, y-r, y+r);

	float z = -999;
	float h;

	for each(auto i in region.distinct_triangles)
	{
		//find height at which tool tip touches triangle on surface.  If no touches, false is returned
		//geometry::triangle<float, 3> &t = region.get_triangle(i);
		//if (geometry::min_height_sphere_on_triangle(x, y, r, t.p1, t.p2, t.p3, h))
		//{
		//	z = fmaxf(z, h);
		//}

		if (region.get_triangle(i).min_height_sphere(x, y, r, h))
		{
			z = fmaxf(z, h);
		}
	}
	for each(auto i in region.distinct_edges)
	{
		//find height at which tool tip touches triangle on edge.  If no touches, false is returned
		//geometry::line<float, 3> &e = region.get_edge(i);
		//if (geometry::min_height_sphere_on_line(x, y, r, e.p1, e.p2, h))
		//{
		//	z = fmaxf(z, h);
		//}
		if (region.get_edge(i).min_height_sphere(x, y, r, h))
		{
			z = fmaxf(z, h);
		}
	}
	for each(auto i in  region.distinct_points)
	{
		//find height at which tool tip touches triangle on vertices.  If no touches, false is returned
		if (geometry::min_height_sphere_on_point(x, y, r, region.get_point(i), h))
		{
			z = fmaxf(z, h);
		}
	}

	if (z > -999)
	{
		return z;
	}
	else
	{
		return default_height();
	}
}


float toolpath_mesh_height::default_height()
{
	return default_z;
}
