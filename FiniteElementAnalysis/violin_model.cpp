#include "violin_model.h"


violin_model::violin_model()
{
}


violin_model::~violin_model()
{
}




void violin_ribs::scale_model(double ratio){

	for (auto it = curves.begin(); it != curves.end(); ++it){
		for (int i = 0; i < it->second->nPoints(); i++){
			it->second->item(i)->operator*=(ratio);
		}
	}
	for (auto it = floats.begin(); it != floats.end(); ++it){
		it->second *=ratio;
	}
}
void violin_model::scale_model(double ratio){
	ribs->scale_model(ratio);
}



base2f::base2f()
{

}
base2f::~base2f()
{
}




point2f::point2f() : vec2<float>()
{
}
point2f::point2f(float u, float v) : vec2<float>(u, v)
{
}
point2f::point2f(base2f* parent, float u, float v) : vec2<float>(u, v)
{
	this->parent = parent;
}

void point2f::Translate(float u, float v)
{
	x += u;
	y += v;
}
void point2f::Translate(point2f p)
{
	x += p.x;
	y += p.y;
}
void point2f::SetParent(base2f* parent)
{
	this->parent = parent;
}
base2f* point2f::GetParent()
{
	return parent;
}
float point2f::GetDistance(point2f& point)
{
	return sqrtf((point.x - x) * (point.x - x) + (point.y - y) * (point.y - y));
}
//float point2f::X()
//{
//	return x;
//}
//float point2f::Y()
//{
//	return y;
//}



point3f::point3f() : vec3<float>()
{
}
point3f::point3f(float u, float v,float w) : vec3<float>(u, v, w)
{
}
point3f::point3f(const point3f &v) : vec3<float>(v)
{
}



circle2f::circle2f(){}
circle2f::circle2f(point2f p, float d)
{
	centre = p;
	diameter = d;
}
circle2f::circle2f(float x, float y, float d)
{
	centre = point2f(x,y);
	diameter = d;
}
circle2f::circle2f(const circle2f &c)
{
	centre = point2f(c.centre);
	diameter = c.diameter;

}



curve2f::curve2f(int order, point2f point) : curve<point2f>(order, point)
{
	for (int i = 0; i < nPoints(); i++)
	{
		item(i)->SetParent(this);
	}
}

curve2f::curve2f(int p, int n, int lknot, double* knot, point2f* points) : curve<point2f>(p,n,lknot,knot,points)
{
	for (int i = 0; i < nPoints(); i++)
	{
		item(i)->SetParent(this);
	}
}

void curve2f::append(point2f point)
{
	bspline::curve < point2f >::append(point);
	item(nPoints() - 1)->SetParent(this);
}



