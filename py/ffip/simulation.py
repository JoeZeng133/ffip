from ffip import Vector3, Medium

class Simulation:
    def __init__(self,
                 cell_size,
                 resolution,
                 geometry=[],
                 sources=[],
                 dimensions=3,
                 boundary_layers=[],
                 symmetries=[],
                 default_material=Medium(),
                 k_point=False,
                 dispersive_materials=[],
                 progress_interval=4,
                 Courant=0.5,
                 output_volume=None,
                 load_structure='',
                 geometry_center = Vector3()):
                 pass

    
class PML:
    def __init__(self):
        pass


    



