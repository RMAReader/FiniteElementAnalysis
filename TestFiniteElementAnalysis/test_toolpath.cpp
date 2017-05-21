#include "test_base.h"
#include "geometry.h"
#include "toolpath_base.h"



void toolpath_build_scanlines::run()
{
	name = "toolpath_build_scanlines";
	float tolerance = 0.0001;
	std::vector<geometry::line<float, 2>> border;

	//outline
	geoVEC2F p1(std::array < float, 2 > {{0, 0}});
	geoVEC2F p2(std::array < float, 2 > {{0, 10}});
	geoVEC2F p3(std::array < float, 2 > {{10, 20}});
	geoVEC2F p4(std::array < float, 2 > {{20, 10}});
	geoVEC2F p5(std::array < float, 2 > {{20, 0}});
	geoVEC2F p6(std::array < float, 2 > {{5, 5}});
	border.push_back(geometry::line< float, 2 >(p1, p2));
	border.push_back(geometry::line< float, 2 >(p2, p3));
	border.push_back(geometry::line< float, 2 >(p3, p4));
	border.push_back(geometry::line< float, 2 >(p4, p5));
	border.push_back(geometry::line< float, 2 >(p5, p6));
	border.push_back(geometry::line< float, 2 >(p6, p1));

	//island
	geoVEC2F p7(std::array < float, 2 > {{10, 10}});
	geoVEC2F p8(std::array < float, 2 > {{10, 15}});
	geoVEC2F p9(std::array < float, 2 > {{15, 10}});
	border.push_back(geometry::line< float, 2 >(p7, p8));
	border.push_back(geometry::line< float, 2 >(p8, p9));
	border.push_back(geometry::line< float, 2 >(p9, p7));

	float step_y = 3;
	float tool_diameter = 6;
	std::vector<scan_line<float>> scan_lines;

	toolpath_base::build_scanlines(scan_lines, border, step_y);
	
	std::vector<scan_line<float>> expected_scan_lines
	{
		scan_line<float>(0.00000f, std::vector < float > { 0, 0, 20, 20 }),
		scan_line<float>(2.85714f, std::vector < float > { 0, 2.85714f, 11.4286f, 20 }),
		scan_line<float>(5.71429f, std::vector < float > { 0, 20 }),
		scan_line<float>(8.57143f, std::vector < float > { 0, 20 }),
		scan_line<float>(11.4286f, std::vector < float > { 1.42857f, 10, 13.5714f, 18.5714f }),
		scan_line<float>(14.2857f, std::vector < float > { 4.28572f, 10, 10.7143f, 15.7143f }),
		scan_line<float>(17.1429f, std::vector < float > { 7.14286f, 12.8571f }),
		scan_line<float>(20.0000f, std::vector < float > { 10, 10 }),
	};
	
	for (int i = 0; i < scan_lines.size(); i++)
	{
		AssertAreEqual(scan_lines[i].y, expected_scan_lines[i].y, tolerance);
		for (int j = 0; j < scan_lines[i].x.size(); j++)
		{
			AssertAreEqual(scan_lines[i].x[j], expected_scan_lines[i].x[j], tolerance);
		}
	}

}



void scanline_path_2D::run()
{
	name = "toolpath_scanline_path_2D";
	float tolerance = 0.000001;

	float tool_diameter = 6;
	std::vector<scan_line<float>> scan_lines
	{
		scan_line<float>(0, std::vector < float > { 0, 20 }),
		scan_line<float>(5, std::vector < float > { 1, 21 }),
		scan_line<float>(10, std::vector < float > { 2, 8, 9, 22 }),
		scan_line<float>(15, std::vector < float > { 3, 8, 9, 23 }),
		scan_line<float>(20, std::vector < float > { 4, 24 }),
		scan_line<float>(25, std::vector < float > { 5, 25, 40,50 }),
		scan_line<float>(30, std::vector < float > { 6, 26, 46, 59 }),
	};

	std::vector<std::vector<geoVEC2F>> path;
	toolpath_base::scanline_path_2D(&path, scan_lines, tool_diameter);


	//expected results
	std::vector<std::vector<geoVEC2F>> expected_path
	{
		std::vector<geoVEC2F>{
			geoVEC2F(std::array < float, 2 > {{0, 0}}),
			geoVEC2F(std::array < float, 2 > {{20, 0}}),
			geoVEC2F(std::array < float, 2 > {{21, 5}}),
			geoVEC2F(std::array < float, 2 > {{1, 5}}),
			geoVEC2F(std::array < float, 2 > {{2, 10}}),
			geoVEC2F(std::array < float, 2 > {{8, 10}}),
			geoVEC2F(std::array < float, 2 > {{8, 15}}),
			geoVEC2F(std::array < float, 2 > {{3, 15}}),
			geoVEC2F(std::array < float, 2 > {{4, 20}}),
			geoVEC2F(std::array < float, 2 > {{24, 20}}),
			geoVEC2F(std::array < float, 2 > {{25, 25}}),
			geoVEC2F(std::array < float, 2 > {{5, 25}}),
			geoVEC2F(std::array < float, 2 > {{6, 30}}),
			geoVEC2F(std::array < float, 2 > {{26, 30}}),
		},
		std::vector<geoVEC2F>{
			geoVEC2F(std::array < float, 2 > {{9, 10}}),
			geoVEC2F(std::array < float, 2 > {{22, 10}}),
			geoVEC2F(std::array < float, 2 > {{23, 15}}),
			geoVEC2F(std::array < float, 2 > {{9, 15}}),
		},
		std::vector<geoVEC2F>{
			geoVEC2F(std::array < float, 2 > {{40, 25}}),
			geoVEC2F(std::array < float, 2 > {{50, 25}}),
			geoVEC2F(std::array < float, 2 > {{59, 30}}),
			geoVEC2F(std::array < float, 2 > {{46, 30}}),
		},
	};
	
	for (int r = 0; r < path.size(); r++)
	{
		for (int p = 0; p < path[r].size(); p++)
		{
			AssertAreEqual(path[r][p][0], expected_path[r][p][0], tolerance);
			AssertAreEqual(path[r][p][1], expected_path[r][p][1], tolerance);
		}
	}
}


void toolpath_scanning_path_3D::run()
{
	name = "toolpath_scanning_path_3D";
	float tolerance = 0.000001;

	geometry::mesh3f mesh;
	geoVEC3F p1(std::array < float, 3 > {{0, 0, 0}});
	geoVEC3F p2(std::array < float, 3 > {{1, 0, 0}});
	geoVEC3F p3(std::array < float, 3 > {{0, 2, 0}});
	mesh.add_triangle(geometry::triangle<float, 3>(p1, p2, p3));
	
	int nx = 1;
	int ny = 1;
	float radius = 0.2;
	float default_z = 30;
	float step = 0.1;
	toolpath_mesh_height h(&mesh, nx, ny, radius, default_z);

	std::vector<std::vector<geoVEC2F>> path_2D
	{
		std::vector < geoVEC2F > {
			geoVEC2F(std::array < float, 2 > {{-0.1, 0.1}}),
			geoVEC2F(std::array < float, 2 > {{1.1, 0.1}}),
		},
		std::vector < geoVEC2F > {
			geoVEC2F(std::array < float, 2 > {{-0.1, 1.0}}),
			geoVEC2F(std::array < float, 2 > {{1.1, 1.0}}),
		}
	};

	std::vector<geoVEC3F> path_3D;
	toolpath_base::scanning_path_3D(&path_3D, path_2D, step, &h);

	for each (auto p in path_3D)
	{
		std::cout << p[0] << "," << p[1] << "," << p[2] << std::endl;
	}
	
}