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

	/* disk geometry*/
	class Disk : public Primitive {
	private:
		fVec3 center;
		real radius;
		real height;
		Direction dir;

	public:
		Disk(const fVec3& center, const real radius, const real height, const Direction dir);
		Disk() = default;
		Disk(const Disk&) = default;
		Disk& operator=(const Disk&) = default;
		Disk(Disk&&) = default;
		Disk& operator=(Disk&&) = default;

		bool is_in_interior(const fVec3& p) const override;
		bool is_in_closure(const fVec3& p) const override;
		bool is_in_exterior(const fVec3& p) const override;
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
		fVec3 get_p1() const;
		fVec3 get_p2() const;
		fVec3 get_center() const;
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
	
	/* Constructive Solid Geometry structure
	   it does not own any resources
	 */
	struct Geometry_Node : public Primitive{
		enum class operation {plus, minus, mult};
		
		operation op;
		Geometry_Node const* left{nullptr};
		Geometry_Node const* right{nullptr};
		Primitive const* val{nullptr};
		
		Geometry_Node(Primitive const* _val);
		Geometry_Node(const operation op, Geometry_Node const& left, Geometry_Node const& right);	//make copy of left, right nodes
		
		Geometry_Node(const Geometry_Node&) = default;				//copy
		Geometry_Node& operator=(const Geometry_Node&) = default;
		Geometry_Node(Geometry_Node&&) = default;					//move
		Geometry_Node& operator=(Geometry_Node&&) = default;
		
		bool is_in_interior(const fVec3& p) const override;
		bool is_in_closure(const fVec3& p) const override;
		bool is_in_exterior(const fVec3& p) const override;
	};
	
	
	Geometry_Node operator*(Geometry_Node const& l, Geometry_Node const& r);	//intersection
	Geometry_Node operator+(Geometry_Node const& l, Geometry_Node const& r);	//union
	Geometry_Node operator-(Geometry_Node const& l, Geometry_Node const& r);	//difference
	
	/* a region with materials */
	class Solid {
	public:
		virtual ~Solid() {}
		virtual bool update_weights(const fVec3& p, Medium_Voxel& weights) const = 0;			//get weights (R^n) of each material at each point, used to calculate mixing materials, if it is outside, weights are 0^n, return 1 if it is inside
	};
	
	/* a box region with inhomogeneous material property specified by e = rho * e1 + (1 - rho) * e2*/
	class Inhomogeneous_Box : public Box, public Solid {
	private:
		Medium const* medium1;
		Medium const* medium2;
		
		GriddedInterp interp;
	public:
		Inhomogeneous_Box(Medium const* m1, Medium const*  m2,  const std::string& filename);	//given medium1, give medium2 and filename of the interpolation data of rho
		
		real get_density(const fVec3& p) const;			//return rho at a given point
		/* override functions*/
		bool update_weights(const fVec3& p, Medium_Voxel& weights) const override;
	};
	
	class Homogeneous_Object : public Solid, public Geometry_Node {
	private:
		Medium const* medium;
		
	public:
		Homogeneous_Object(Medium const*  m, const Geometry_Node& _base);
		bool update_weights(const fVec3& p, Medium_Voxel& weights) const override;
	};
	
}
