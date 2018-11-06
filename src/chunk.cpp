#include <chunk.hpp>

namespace ffip {
	Chunk::Chunk(const iVec3& _sim_p1, const iVec3& _sim_p2, const iVec3& _ch_p1, const iVec3& _ch_p2, real _dx, real _dt):  sim_p1(_sim_p1), sim_p2(_sim_p2), ch_p1(_ch_p1), ch_p2(_ch_p2), dx(_dx), dt(_dt) {
		
		ch_origin = ch_p1 - iVec3{1, 1, 1};
		ch_dim = ch_p2 - ch_p1 + iVec3{3, 3, 3};
		
		ch_jump_x = ch_dim.y * ch_dim.z;
		ch_jump_y = ch_dim.z;
		ch_jump_z = 1;
		
		//jumps for use in averaging fields around a point
//		for(int i = 0; i < 8; ++i) {
//			jump[i].push_back(0);
//			if(i & 1) {
//				int tmp = jump[i].size();
//				for(int j = 0; j < tmp; j++) {
//					jump[i].push_back(jump[i][j] + ch_jump_x);
//					jump[i][j] -= jump_x;
//				}
//			}
//
//			if(i & 1) {
//				int tmp = jump[i].size();
//				for(int j = 0; j < tmp; j++) {
//					jump[i].push_back(jump[i][j] + jump_y);
//					jump[i][j] -= jump_y;
//				}
//			}
//
//			if(i & 1) {
//				int tmp = jump[i].size();
//				for(int j = 0; j < tmp; j++) {
//					jump[i].push_back(jump[i][j] + jump_z);
//					jump[i][j] -= jump_z;
//				}
//			}
//		}
		
	}
	
	iVec3 Chunk::get_dim() const{
		return ch_dim;
	}
	
	iVec3 Chunk::get_origin() const {
		return ch_origin;
	}
	
	iVec3 Chunk::get_p1() const {
		return ch_p1;
	}
	
	iVec3 Chunk::get_p2() const {
		return ch_p2;
	}
	
	int Chunk::get_index_ch(const iVec3& p) const {
		return (p.x - ch_origin.x) * ch_jump_x + (p.y - ch_origin.y) * ch_jump_y + (p.z - ch_origin.z);
	}
	
	void Chunk::set_medium_point(const iVec3 &point, Medium_Internal * const medium_internal) {
		int index = get_index_ch(point);
		medium_chunk[index] = medium_internal;
	}
	
	void Chunk::add_source_internal(Source_Internal *const source) {
		source_list.push_back(source);
	}
	
	void Chunk::PML_init(const real_arr &kx, const real_arr &ky, const real_arr &kz, const real_arr &bx, const real_arr &by, const real_arr &bz, const real_arr &cx, const real_arr &cy, const real_arr &cz) {
		
		auto ch_p1_ch = ch_p1 - ch_origin;
		auto ch_p2_ch = ch_p2 - ch_origin;
		auto ch_p1_sim = ch_p1 - sim_p1;
		
		for(int x = ch_p1_ch.x, x_sim = ch_p1_sim.x; x <= ch_p2_ch.x; x++, x_sim++)
			for(int y = ch_p1_ch.y, y_sim = ch_p1_sim.y; y <= ch_p2_ch.y; y++, y_sim++)
				for(int z = ch_p1_ch.z, z_sim = ch_p1_sim.z; z <= ch_p2_ch.z; z++, z_sim++) {
					
					int index = x * ch_jump_x + y * ch_jump_y + z * ch_jump_z;
					auto tmp = iVec3{x_sim, y_sim, z_sim};
					
					Coord_Type ctype = tmp.get_type(ch_origin.get_type());
					
					switch (ctype) {
						case Ex :
							e_PML.push_back(PML_Point(index, ch_jump_y, ch_jump_z, by[y_sim], bz[z_sim], cy[y_sim], cz[z_sim]));
							break;
							
						case Hx :
							h_PML.push_back(PML_Point(index, ch_jump_y, ch_jump_z, by[y_sim], bz[z_sim], cy[y_sim], cz[z_sim]));
							break;
						
						case Ey :
							e_PML.push_back(PML_Point(index, ch_jump_z, ch_jump_x, bz[z_sim], bx[x_sim], cz[z_sim], cx[x_sim]));
							break;
							
						case Hy :
							h_PML.push_back(PML_Point(index, ch_jump_z, ch_jump_x, bz[z_sim], bx[x_sim], cz[z_sim], cx[x_sim]));
							break;
							
						case Ez :
							e_PML.push_back(PML_Point(index, ch_jump_x, ch_jump_y, bx[x_sim], by[y_sim], cx[y_sim], cy[z_sim]));
							break;
							
						case Hz :
							h_PML.push_back(PML_Point(index, ch_jump_x, ch_jump_y, bx[x_sim], by[y_sim], cx[y_sim], cy[z_sim]));
							break;
							
						default:
							break;
					}
				}
		
		
		this->kx = kx;
		this->kx.insert(this->kx.begin(), 1);
		this->kx.insert(this->kx.end(), 1);
		
		this->ky = ky;
		this->ky.insert(this->ky.begin(), 1);
		this->ky.insert(this->ky.end(), 1);
		

	}
	
	void Chunk::update_Jd(real time) {
		update_JMd_helper(Ex, Ey, Ez, e_PML);
		for(auto& x : source_list) {
			x->get_Jd(time);
			x->update_Md(jmd);
		}
	}
	
	void Chunk::update_Md(real time) {
		update_JMd_helper(Hx, Hy, Hz, h_PML);
		for(auto& x : source_list) {
			x->get_Md(time);
			x->update_Jd(jmd);
		}
	}
	
	void Chunk::update_JMd_helper(const Coord_Type Fx, const Coord_Type Fy, const Coord_Type Fz, std::vector<PML_Point>& PML) {
		/* Curl updates without PML */
		auto tmp = get_component_interior(ch_p1, ch_p2, Fx);
		auto p1_ch = tmp.first - ch_origin;
		auto p2_ch = tmp.second - ch_origin;
		
		for(int x = p1_ch.x; x <= p2_ch.x; x += 2)
			for(int y = p1_ch.y; y <= p2_ch.y; y += 2)
				for(int z = p1_ch.z; z <= p2_ch.z; z += 2) {
					int index = x * ch_jump_x + y * ch_jump_y + z * ch_jump_z;
					
					eh[index] = (eh[index + ch_jump_y] - eh[index - ch_jump_y]) / dx / ky[y] - (eh[index + ch_jump_z] - eh[index - ch_jump_z]) / dx / kz[z];
				}
		
		tmp = get_component_interior(ch_p1, ch_p2, Fy);
		p1_ch = tmp.first - ch_origin;
		p2_ch = tmp.second - ch_origin;
		
		for(int x = p1_ch.x; x <= p2_ch.x; x += 2)
			for(int y = p1_ch.y; y <= p2_ch.y; y += 2)
				for(int z = p1_ch.z; z <= p2_ch.z; z += 2) {
					int index = x * ch_jump_x + y * ch_jump_y + z * ch_jump_z;
					
					eh[index] = (eh[index + ch_jump_z] - eh[index - ch_jump_z]) / dx / kz[z] - (eh[index + ch_jump_x] - eh[index - ch_jump_z]) / dx / kx[x];
				}
		
		tmp = get_component_interior(ch_p1, ch_p2, Fz);
		p1_ch = tmp.first - ch_origin;
		p2_ch = tmp.second - ch_origin;
		
		for(int x = p1_ch.x; x <= p2_ch.x; x += 2)
			for(int y = p1_ch.y; y <= p2_ch.y; y += 2)
				for(int z = p1_ch.z; z <= p2_ch.z; z += 2) {
					int index = x * ch_jump_x + y * ch_jump_y + z * ch_jump_z;
					
					eh[index] = (eh[index + ch_jump_x] - eh[index - ch_jump_x]) / dx / kx[x] - (eh[index + ch_jump_y] - eh[index - ch_jump_y]) / dx / ky[y];
				}
		
		for(auto& x : PML) {
			if(x.c_pos > 0)
				x.psi_pos = x.b_pos * x.psi_pos + x.c_pos * (eh[x.index + x.jump_pos] - eh[x.index - ch_jump_y]) / dx;
			if(x.c_neg > 0)
				x.psi_neg = x.b_neg * x.psi_neg + x.c_neg * (eh[x.index + x.jump_neg] - eh[x.index - ch_jump_z]) / dx;
			jmd[x.index] += x.psi_pos - x.psi_neg;
		}
	}
	
	void Chunk::update_H2B(real time) {
		update_DEHB_helper(Ex);
		update_DEHB_helper(Ey);
		update_DEHB_helper(Ez);
	}
	
	void Chunk::update_D2E(real time) {
		update_DEHB_helper(Hx);
		update_DEHB_helper(Hy);
		update_DEHB_helper(Hz);
	}
	
	void Chunk::update_DEHB_helper(const Coord_Type F) {
		auto tmp = get_component_interior(ch_p1, ch_p2, F);
		auto p1_ch = tmp.first - ch_origin;
		auto p2_ch = tmp.second - ch_origin;
		
		for(int x = p1_ch.x; x <= p2_ch.x; x += 2)
			for(int y = p1_ch.y; y <= p2_ch.y; y += 2)
				for(int z = p1_ch.z; z <= p2_ch.z; z += 2) {
					int index = x * ch_jump_x + y * ch_jump_y + z * ch_jump_z;
					real tmp = eh[index];
					
					eh[index] = medium_chunk[index]->update_field(eh[index], eh1[index], jmd[index], dispersive_field_chunk[index]);
					eh1[index] = tmp;
				}
	}
	
	void Chunk::update_padded_E(real time) {}
	
	void Chunk::update_padded_H(real time) {}
	
	real Chunk::at(const fVec3& p, const Coord_Type ctype) const{
		return operator()(p * (2 / dx), ctype);
	}

	real Chunk::operator()(const iVec3 &p) const {
		if (!ElementWise_Less_Eq(ch_p1, p) || !ElementWise_Less_Eq(p, ch_p2))
			return 0;

		return eh[(p.x - ch_origin.x) * ch_jump_x + (p.y - ch_origin.y) * ch_jump_y + (p.z - ch_origin.z) * ch_jump_z];
	}

	real Chunk::operator()(const iVec3& p, const Coord_Type ctype) const {
		if (!ElementWise_Less_Eq(ch_p1, p) || !ElementWise_Less_Eq(p, ch_p2))
			return 0;

		return ave(ctype ^ p.get_type(), get_index_ch(p));
	}

	real Chunk::interp_helper(const real* data, const real w) const{
		return data[0] * (1 - w) + data[1] * w;
	}

	real Chunk::operator()(const fVec3& p, const Coord_Type ctype) const {
		if (!ElementWise_Less_Eq(ch_p1, p) || !ElementWise_Less_Eq(p, ch_p2))
			return 0;

		static real f[8];
		
		auto tmp = get_component_closure(p, p, ctype);
		int base_index_ch = get_index_ch(tmp.first);

		for (int i = 0; i < 8; ++i)
			f[i] = eh[base_index_ch + !(i & 1) * ch_jump_x + !(i & 2) * ch_jump_y + !(i & 4) * ch_jump_z];

		return interp_helper(f, 
			(p.z - tmp.first.z) / 2,
			(p.y - tmp.first.y) / 2,
			(p.x - tmp.first.x) / 2
			);
	}

	template<>
	int Chunk::ave_helper(const int bit, const int index, const int jump) const{
		if(bit & 1) {
			return eh[index - jump] + eh[index + jump];
		}else{
			return eh[index];
		}
	}
	
	int Chunk::ave(const int bit, const int index) const{
		return ave_helper(bit, index, ch_jump_x, ch_jump_y, ch_jump_z) / jump[bit].size();
	}
	
}
