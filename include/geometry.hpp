#pragma once

#include <utility.hpp>
#include <medium.hpp>
#include <memory>

namespace ffip {
	class Primitive {
	public:
		virtual ~Primitive() {}
		Primitive() = default;
		virtual bool is_in_interior(const fVec3& p) const = 0;							//return whether the given point is inside the geometry
		virtual bool is_in_closure(const fVec3& p) const = 0;
		virtual bool is_in_exterior(const fVec3& p) const = 0;
	};
	
	/* box geometry */
	class Box : public Primitive {
	private:
		real lenx, leny, lenz;
		fVec3 center;
		fVec3 p1, p2;
		
	public:
		Box(const fVec3& _center, const real _lenx, const real _leny, const real _lenz);
		Box(const fVec3& _p1, const fVec3& _p2);
		Box() = default;
		Box(const Box&) = default;					//copy
		Box& operator=(const Box&) = default;
		Box(Box&&) = default;							//move
		Box& operator=(Box&&) = default;
		
		/* re-construction */
		void init(const fVec3& _center, const real _lenx, const real _leny, const real _lenz);
		void init(const fVec3& _p1, const fVec3& _p2);
		
		/* override functions*/
		bool is_in_interior(const fVec3& p) const override;
		bool is_in_closure(const fVec3& p) const override;
		bool is_in_exterior(const fVec3& p) const override;
		
		/* other functions*/
		fVec3 get_p1();
		fVec3 get_p2();
		fVec3 get_center();
	};
	
	/* sphere geometry*/
	class Sphere : public Primitive {
	private:
		real radius;
		fVec3 center;
		
	public:
		Sphere() = default;
		Sphere(const fVec3& _center, const real _radius);
		Sphere(const Sphere&) = default;							//copy
		Sphere& operator=(const Sphere&) = default;
		Sphere(Sphere&&) = default;								//move
		Sphere& operator=(Sphere&&) = default;
		
		/* re-construction */
		void init(const fVec3& _center, const real _radius);
		
		/* override functions*/
		bool is_in_interior(const fVec3& p) const override;
		bool is_in_closure(const fVec3& p) const override;
		bool is_in_exterior(const fVec3& p) const override;
		
		/* other functions*/
		fVec3 get_center();
		real get_radius();
	};
	
	/* Constructive Solid Geometry structure*/
	struct Geometry_Node : public Primitive{
		enum class operation {plus, minus, mult};
		
		operation op;
		Geometry_Node* left{nullptr};
		Geometry_Node* right{nullptr};
		Primitive* val{nullptr};
		
		Geometry_Node(Primitive* const _val);
		Geometry_Node(const operation op, Geometry_Node const& left, Geometry_Node const& right);	//make copy of left, right nodes
		
		Geometry_Node(const Geometry_Node&) = default;				//copy
		Geometry_Node& operator=(const Geometry_Node&) = default;
		Geometry_Node(Geometry_Node&&) = default;					//move
		Geometry_Node& operator=(Geometry_Node&&) = default;
		
		bool is_in_interior(const fVec3& p) const override;
		bool is_in_closure(const fVec3& p) const override;
		bool is_in_exterior(const fVec3& p) const override;
	};
	
	void reset(Geometry_Node* x);			//release memory
	
	Geometry_Node operator*(Geometry_Node const& l, Geometry_Node const& r);	//intersection
	Geometry_Node operator+(Geometry_Node const& l, Geometry_Node const& r);	//union
	Geometry_Node operator-(Geometry_Node const& l, Geometry_Node const& r);	//difference
	
	/* a region with materials */
	class Solid {
	public:
		virtual ~Solid() {}
		virtual bool get_weights(const fVec3& p, std::vector<real>& weights) const;			//get weights (R^n) of each material at each point, used to calculate mixing materials, if it is outside, weights are 0^n, return 1 if it is inside
	};
	
	/* a box region with inhomogeneous material property specified by e = rho * e1 + (1 - rho) * e2*/
	class Inhomogeneous_Box : public Solid, public Box {
	private:
		Medium_Type medium1;
		Medium_Type medium2;
		GriddedInterp interp;
		
	public:
		Inhomogeneous_Box(const Medium_Type& _medium1, const Medium_Type& _medium2, const std::string& filename);	//given medium1, give medium2 and filename of the interpolation data of rho
		
		real get_density(const fVec3& p) const;			//return rho at a given point
		/* override functions*/
		bool get_weights(const fVec3& p, std::vector<real>& weights) const override;
	};
	
	class Homogeneous_Object : public Solid, public Geometry_Node {
	private:
		Medium_Type medium;
		
	public:
		Homogeneous_Object(const Geometry_Node& _base, const Medium_Type& _medium);
		bool get_weights(const fVec3& p, std::vector<real>& weights) const override;
	};
	
	/* type factory */
	template<typename... Args>
	Geometry_Node make_sphere(Args&&... args) {
		return Geometry_Node{new Sphere{std::forward<Args>(args)...}};
	}
	
	template<typename... Args>
	Geometry_Node make_box(Args&&... args) {
		return Geometry_Node{new Box{std::forward<Args>(args)...}};
	}
	
	Solid* make_solid(const Medium_Type& medium1, const Medium_Type& medium2, const std::string& filename);	//make inhomogeneous regions
	Solid* make_solid(const Geometry_Node& geom, const Medium_Type& medium);											//make homogeneous regions
}
