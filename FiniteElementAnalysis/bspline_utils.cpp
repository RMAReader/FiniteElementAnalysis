//#include "bspline_utils.h"
//
//
//bspline::vec2<float> operator+(bspline::vec2<float>& u, bspline::vec2<float>& v)
//{
//	bspline::vec2<float> r(u.x + v.x, u.y + v.y);
//	return u;
//}
//bspline::vec2<float> operator+(float d, bspline::vec2<float>& v)
//{
//	bspline::vec2<float> r(d + v.x, d + v.y);
//	return r;
//}
//bspline::vec2<float> operator+(bspline::vec2<float>& v, float d)
//{
//	bspline::vec2<float> r(d + v.x, d + v.y);
//	return r;
//}
//bspline::vec2<float> operator*(bspline::vec2<float>& u, bspline::vec2<float>& v)
//{
//	bspline::vec2<float> r(u.x * v.x, u.y * v.y);
//	return u;
//}
//bspline::vec2<float> operator*(float d, bspline::vec2<float>& v)
//{
//	bspline::vec2<float> r(d * v.x, d * v.y);
//	return r;
//}
//bspline::vec2<float> operator*(bspline::vec2<float>& v, float d)
//{
//	bspline::vec2<float> r(d * v.x, d * v.y);
//	return r;
//}
//
//
//bspline::offsetCurve2f::offsetCurve2f(bspline::curve<bspline::vec2<float>>* curve, double x)
//{
//	this->base = curve;
//	this->offset = x;
//}
//void bspline::offsetCurve2f::SetBase(bspline::curve<bspline::vec2<float>>* curve)
//{
//	this->base = curve;
//}
//void bspline::offsetCurve2f::SetOffset(double x)
//{
//	this->offset = x;
//}
//bspline::curve<bspline::vec2<float>>* bspline::offsetCurve2f::GetBase()
//{
//	return this->base;
//}
//double bspline::offsetCurve2f::GetOffset()
//{
//	return this->offset;
//}
//bspline::vec2<float>  bspline::offsetCurve2f::evaluate(double t)
//{
//	double dt = 0.000001;
//
//	if (t < base->minParam() || t > base->maxParam()){ return point2f(0, 0); }
//
//	double u = t + dt;
//	double v = t - dt;
//
//	if (u > base->maxParam()){ u = base->maxParam(); }
//	if (v < base->minParam()){ v = base->minParam(); }
//
//	bspline::vec2<float> p1 = this->base->evaluate(u);
//	bspline::vec2<float> p2 = this->base->evaluate(v);
//	bspline::vec2<float> n(p1.y - p2.y, p2.x - p1.x);
//	if (n.x == 0 && n.y == 0)
//	{
//		return base->evaluate(t);
//	}
//	else
//	{
//		n *= offset / sqrtf(n.x*n.x + n.y*n.y);
//		bspline::vec2<float> p = base->evaluate(t);
//		p += n;
//		return p;
//	}
//}
