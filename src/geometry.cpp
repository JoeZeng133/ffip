#include <geometry.hpp>

namespace ffip {
	/* Box Geometry*/
	Box::Box(const fVec3& _center, const real _lenx, const real _leny, const real _lenz): lenx(_lenx), leny(_leny), lenz(_lenz), center(_center) {
		auto temp = fVec3{lenx, leny, lenz} / 2;
		p1 = center - temp;
		p2 = center + temp;
	}
	
	Box::Box(const fVec3& _p1, const fVec3& _p2): p1(_p1), p2(_p2) {
		center = {(p1.x + p2.x) / 2, (p1.y + p2.y) / 2, (p1.z + p2.z) / 2};
		lenx = p2.x - p1.x;
		leny = p2.y - p1.y;
		lenz = p2.z - p1.z;
	}
	
	void Box::init(const fVec3 &_center, const real _lenx, const real _leny, const real _lenz) {
		center = _center;
		lenx = _lenx;
		leny = _leny;
		lenz = _lenz;
		p1 = center - Vec3<>{lenx / 2, leny / 2, lenz / 2};
		p2 = center + Vec3<>{lenx / 2, leny / 2, lenz / 2};
	}
	
	void Box::init(const fVec3 &_p1, const fVec3 &_p2) {
		p1 = _p1;
		p2 = _p2;
		center = {(p1.x + p2.x) / 2, (p1.y + p2.y) / 2, (p1.z + p2.z) / 2};
		lenx = p2.x - p1.x;
		leny = p2.y - p1.y;
		lenz = p2.z - p1.z;
	}
	
	bool Box::is_in_closure(const fVec3 &p) const {
		return	(p.x >= p1.x) && (p.x <= p2.x) &&
				(p.y >= p1.y) && (p.y <= p2.y) &&
				(p.z >= p1.z) && (p.z <= p2.z);
	}
	
	bool Box::is_in_interior(const fVec3 &p) const {
		return	(p.x > p1.x) && (p.x < p2.x) &&
				(p.y > p1.y) && (p.y < p2.y) &&
				(p.z > p1.z) && (p.z < p2.z);
	}
	
	bool Box::is_in_exterior(const fVec3 &p) const {
		return !is_in_closure(p);
	}
	
	fVec3 Box::get_p1() {
		return p1;
	}
	
	fVec3 Box::get_p2() {
		return p2;
	}
	
	fVec3 Box::get_center() {
		return center;
	}
	/* sphere geometry*/
	
	Sphere::Sphere(const fVec3& _center, const real _radius):radius(_radius), center(_center) {}
	
	void Sphere::init(const fVec3 &_center, const real _radius) {
		center = _center;
		radius = _radius;
	}
	
	bool Sphere::is_in_closure(const fVec3 &p) const {
		return	(p.x - center.x) * (p.x - center.x) +
				(p.y - center.y) * (p.y - center.y) +
				(p.z - center.z) * (p.z - center.z)
				<= radius * radius;
	}
	
	bool Sphere::is_in_interior(const fVec3 &p) const {
		return	(p.x - center.x) * (p.x - center.x) +
				(p.y - center.y) * (p.y - center.y) +
				(p.z - center.z) * (p.z - center.z)
				< radius * radius;
	}
	
	bool Sphere::is_in_exterior(const fVec3 &p) const {
		return !is_in_closure(p);
	}
	
	Geometry_Node::Geometry_Node(Primitive const* _val): val(_val) {}
	
	Geometry_Node::Geometry_Node(const operation _op, Geometry_Node const& _left, Geometry_Node const& _right): op(_op), left(new Geometry_Node{_left}), right(new Geometry_Node{_right}) {}
	
	bool Geometry_Node::is_in_interior(const fVec3 &point) const {
		if(val)
			return val->is_in_interior(point);
		else {
			switch (op) {
				case operation::plus:
					return left->is_in_interior(point) || right->is_in_interior(point);
					break;
					
				case operation::minus:
					return left->is_in_interior(point) && right->is_in_exterior(point);
					break;
					
				case operation::mult:
					return left->is_in_interior(point) && right->is_in_interior(point);
					break;
					
				default:
					return false;
					break;
			}
		}
	}
	
	bool Geometry_Node::is_in_closure(const fVec3 &point) const {
		if(val)
			return val->is_in_closure(point);
		else {
			switch (op) {
				case operation::plus:
					return left->is_in_closure(point) || right->is_in_closure(point);
					break;
					
				case operation::minus:
					return left->is_in_closure(point) && !right->is_in_interior(point);
					break;
					
				case operation::mult:
					return left->is_in_closure(point) && right->is_in_closure(point);
					break;
					
				default:
					return false;
					break;
			}
		}
	}
	
	bool Geometry_Node::is_in_exterior(const fVec3 &point) const {
		return !is_in_closure(point);
	}

	Geometry_Node operator*(Geometry_Node const& l, Geometry_Node const& r) {
		return Geometry_Node(Geometry_Node::operation::mult, new Geometry_Node{l}, new Geometry_Node{r});
	}
	
	Geometry_Node operator+(Geometry_Node const& l, Geometry_Node const& r) {
		return Geometry_Node(Geometry_Node::operation::plus, new Geometry_Node{l}, new Geometry_Node{r});
	}
	
	Geometry_Node operator-(Geometry_Node const& l, Geometry_Node const& r) {
		return Geometry_Node(Geometry_Node::operation::minus, new Geometry_Node{l}, new Geometry_Node{r});
	}
	
	/* Inhomogeneous box*/
	Inhomogeneous_Box::Inhomogeneous_Box(Medium const*  m1, Medium const*  m2,  const std::string& filename): medium1(m1), medium2(m2), interp{filename} {
		Box::init(interp.get_p1(), interp.get_p2());
	}
	
	bool Inhomogeneous_Box::update_weights(const fVec3& p, std::vector<real> &weights) const {
		if (!is_in_interior(p))
			return 0;
		else {
			real density = get_density(p);
			weights[medium1->index] += density;
			weights[medium2->index] += 1 - density;
			return 1;
		}
	}
	
	real Inhomogeneous_Box::get_density(const fVec3 &p) const {
		return interp.request_value(p);
	}
	
	Homogeneous_Object::Homogeneous_Object(Medium const*  m, const Geometry_Node& _base): Geometry_Node(_base), medium(m) {}
	
	bool Homogeneous_Object::update_weights(const fVec3 &p, std::vector<real> &weights) const {
		if (!is_in_interior(p))
			return 0;
		else {
			weights[medium->index] += 1;
			return 1;
		}
	}
	
	/* geometry factory */
	std::vector<Inhomogeneous_Box*> inhom_box_holder;
	std::vector<Homogeneous_Object*> hom_box_holder;
	std::vector<Primitive*> primitive_holder;
	
	Solid const* make_solid(Medium const* m1, Medium const*  m2, const std::string& filename) {
		inhom_box_holder.push_back(new Inhomogeneous_Box{m1, m2, filename});
		return inhom_box_holder.back();
	}
	
	Solid const* make_solid(Medium const* m, const Geometry_Node& geom) {
		hom_box_holder.push_back(new Homogeneous_Object{m, geom});
		return hom_box_holder.back();
	}
}

