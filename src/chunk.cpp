#include <chunk.hpp>

namespace ffip {
	Chunk::Chunk(const iVec3& _sim_p1, const iVec3& _sim_p2, const iVec3& _ch_p1, const iVec3& _ch_p2, real _dx, real _dt):  sim_p1(_sim_p1), sim_p2(_sim_p2), ch_p1(_ch_p1), ch_p2(_ch_p2), dx(_dx), dt(_dt) {
		
		ch_origin = ch_p1 - iVec3{1, 1, 1};
		ch_dim = ch_p2 - ch_p1 + iVec3{3, 3, 3};
		
		ch_jump_x = 1;
		ch_jump_y = ch_dim.x;
		ch_jump_z = ch_dim.x * ch_dim.y;
		
		eh.resize(ch_dim.x * ch_dim.y * ch_dim.z, 0);
		eh1.resize(eh.size(), 0);
		jmd.resize(eh.size(), 0);
		medium_chunk.resize(eh.size(), nullptr);
		dispersive_field_chunk.resize(eh.size(), nullptr);
		
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
		return (p.x - ch_origin.x) * ch_jump_x + (p.y - ch_origin.y) * ch_jump_y + (p.z - ch_origin.z) * ch_jump_z;
	}
	
	void Chunk::set_medium_point(const iVec3 &point, Medium_Internal * const medium_internal) {
		int index = get_index_ch(point);
		medium_chunk[index] = medium_internal;
	}
	
	void Chunk::add_source_internal(Source_Internal *const source) {
		source_list.push_back(source);
	}
	
	template<>
	const int Chunk::get_ch_jump<dir_x_tag>() const {
		return ch_jump_x;
	}
	
	template<>
	const int Chunk::get_ch_jump<dir_y_tag>() const {
		return ch_jump_y;
	}
	
	template<>
	const int Chunk::get_ch_jump<dir_z_tag>() const {
		return ch_jump_z;
	}
	
	template<>
	const real Chunk::get_k<dir_x_tag>(const int id) const {
		return kx[id];
	}
	
	template<>
	const real Chunk::get_k<dir_y_tag>(const int id) const {
		return ky[id];
	}
	
	template<>
	const real Chunk::get_k<dir_z_tag>(const int id) const {
		return kz[id];
	}
	
	void Chunk::PML_init(const real_arr &kx, const real_arr &ky, const real_arr &kz, const real_arr &bx, const real_arr &by, const real_arr &bz, const real_arr &cx, const real_arr &cy, const real_arr &cz) {
		
		auto p1_ch = ch_p1 - ch_origin;
		auto p2_ch = ch_p2 - ch_origin;
		auto p1_sim = ch_p1 - sim_p1;
		auto p2_sim = ch_p2 - sim_p1;
		
		{
			auto itr_ch = my_iterator(p1_ch, p2_ch, Null);
			auto itr_sim = my_iterator(p1_sim, p2_sim, Null);
			
			for(;!itr_ch.is_end();itr_ch.advance(),
				itr_sim.advance()) {
				int index = itr_ch.x * ch_jump_x + itr_ch.y * ch_jump_y + itr_ch.z * ch_jump_z;
				auto tmp = itr_sim.get_vec();
				
				Coord_Type ctype = tmp.get_type(sim_p1.get_type());
				
				if (cx[itr_sim.x] == 0 && cy[itr_sim.y] == 0 && cz[itr_sim.z] == 0) continue;
				
				switch (ctype) {
					case Ex :
						e_PML.push_back(PML_Point(index, ch_jump_y, ch_jump_z, by[itr_sim.y], bz[itr_sim.z], cy[itr_sim.y], cz[itr_sim.z]));
						break;
						
					case Hx :
						h_PML.push_back(PML_Point(index, ch_jump_y, ch_jump_z, by[itr_sim.y], bz[itr_sim.z], cy[itr_sim.y], cz[itr_sim.z]));
						break;
						
					case Ey :
						e_PML.push_back(PML_Point(index, ch_jump_z, ch_jump_x, bz[itr_sim.z], bx[itr_sim.x], cz[itr_sim.z], cx[itr_sim.x]));
						break;
						
					case Hy :
						h_PML.push_back(PML_Point(index, ch_jump_z, ch_jump_x, bz[itr_sim.z], bx[itr_sim.x], cz[itr_sim.z], cx[itr_sim.x]));
						break;
						
					case Ez :
						e_PML.push_back(PML_Point(index, ch_jump_x, ch_jump_y, bx[itr_sim.x], by[itr_sim.y], cx[itr_sim.x], cy[itr_sim.y]));
						break;
						
					case Hz :
						h_PML.push_back(PML_Point(index, ch_jump_x, ch_jump_y, bx[itr_sim.x], by[itr_sim.y], cx[itr_sim.x], cy[itr_sim.y]));
						break;
						
					default:
						break;
				}
			}
		}
		
		this->kx.insert(this->kx.end(), kx.begin() + p1_sim.x, kx.begin() + p2_sim.x + 1);
		this->kx.insert(this->kx.begin(), 1);
		this->kx.insert(this->kx.end(), 1);
		
		this->ky.insert(this->ky.end(), ky.begin() + p1_sim.y, ky.begin() + p2_sim.y + 1);
		this->ky.insert(this->ky.begin(), 1);
		this->ky.insert(this->ky.end(), 1);
		
		this->kz.insert(this->kz.end(), kz.begin() + p1_sim.z, kz.begin() + p2_sim.z + 1);
		this->kz.insert(this->kz.begin(), 1);
		this->kz.insert(this->kz.end(), 1);
	}
	
	
	void Chunk::update_Jd(const real time) {
		update_JMd_helper<ex_tag>();
		update_JMd_helper<ey_tag>();
		update_JMd_helper<ez_tag>();
		
		for(auto& item : e_PML) {
			if(item.c_pos != 0)
				item.psi_pos = item.b_pos * item.psi_pos + item.c_pos * (eh[item.index + item.jump_pos] - eh[item.index - item.jump_pos]) / dx;
			
			if(item.c_neg != 0)
				item.psi_neg = item.b_neg * item.psi_neg + item.c_neg * (eh[item.index + item.jump_neg] - eh[item.index - item.jump_neg]) / dx;
			
			jmd[item.index] += item.psi_pos - item.psi_neg;
		}
		
		for(auto& x : source_list) {
			x->get_Jd(time);
			x->update_Jd(jmd);
		}
	}
	
	void Chunk::update_Md(const real time) {
		update_JMd_helper<hx_tag>();
		update_JMd_helper<hy_tag>();
		update_JMd_helper<hz_tag>();
		
		for(auto& item : h_PML) {
			if(item.c_pos != 0)
				item.psi_pos = item.b_pos * item.psi_pos + item.c_pos * (eh[item.index + item.jump_pos] - eh[item.index - item.jump_pos]) / dx;
			
			if(item.c_neg != 0)
				item.psi_neg = item.b_neg * item.psi_neg + item.c_neg * (eh[item.index + item.jump_neg] - eh[item.index - item.jump_neg]) / dx;
			
			jmd[item.index] += item.psi_pos - item.psi_neg;
		}
		
		for(auto& x : source_list) {
			x->get_Md(time);
			x->update_Md(jmd);
		}
	}
	
	void Chunk::update_B2H(const real time) {
		update_DEHB_helper(Hx);
		update_DEHB_helper(Hy);
		update_DEHB_helper(Hz);
	}
	
	void Chunk::update_D2E(const real time) {
		update_DEHB_helper(Ex);
		update_DEHB_helper(Ey);
		update_DEHB_helper(Ez);
	}
	
	void Chunk::update_DEHB_helper(const Coord_Type F) {
		auto tmp = get_component_interior(ch_p1, ch_p2, F);
		auto p1_ch = tmp.first - ch_origin;
		auto p2_ch = tmp.second - ch_origin;
		
		for(auto itr = my_iterator(p1_ch, p2_ch, p1_ch.get_type()); !itr.is_end(); itr.advance()) {
					int index = itr.x * ch_jump_x + itr.y * ch_jump_y + itr.z * ch_jump_z;
					real tmp = eh[index];
					
					eh[index] = medium_chunk[index]->update_field(eh[index], eh1[index], jmd[index], dispersive_field_chunk[index]);
					eh1[index] = tmp;
				}
	}
	
	void Chunk::update_padded_E(const real time) {}
	
	void Chunk::update_padded_H(const real time) {}
	
	real Chunk::at(const fVec3& p, const Coord_Type ctype) const{
		return operator()(p / (dx / 2), ctype);
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
		
		auto tmp = get_nearest_point<side_low_tag>(p, ctype);
		int base_index_ch = get_index_ch(tmp);

		for (int i = 0; i < 8; ++i)
			f[i] = eh[base_index_ch + (i & 1? ch_jump_x << 1 : 0) + (i & 2? ch_jump_y << 1 : 0) + (i & 4? ch_jump_z << 1 : 0)];

		return interp_helper(f, 
			(p.z - tmp.z) / 2,
			(p.y - tmp.y) / 2,
			(p.x - tmp.x) / 2
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
