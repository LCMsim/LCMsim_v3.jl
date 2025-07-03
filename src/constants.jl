
const NO_NEIGHBOUR = -9

const KEY_ALPHA = "alpha"
const KEY_PERMEABILITY = "permeability"
const KEY_PERMEABILITY_NOISE = "permeability_noise"
const KEY_POROSITY = "porosity"
const KEY_POROSITY_NOISE = "porosity_noise"
const KEY_REFERENCE_DIRECTION = "reference_direction"
const KEY_THICKNESS = "thickness"
const KEY_TYPE = "type"
const KEY_POROSITY_1 = "porosity_1"
const KEY_P_1 = "p_1"
const KEY_POROSITY_2 = "porosity_2"
const KEY_P_2 = "p_2"
const KEY_PERMEABILITY_Z = "permeability_Z"
const KEY_THICKNESS_DM = "thickness_DM"
const KEY_POROSITY_DM = "porosity_DM"
const KEY_PERMEABILITY_DM = "permeability_DM"
const KEY_PERMEABILITY_DM_Z = "permeability_DM_Z"
const KEY_N_CK = "n_CK"
const KEY_PART_ID = "part_id"
const KEY_NAME = "name"


# 
const PERMEABILITY_MAX = 1.0
const PERMEABILITY_MIN = 0.0

const GAMMA_A = 1.0
const GAMMA_INIT = 0.0
const V_A = 0.
const V_INIT = 0.
const U_A = 0.
const U_INIT = 0.

# HDF group path constants
const GROUP_SIMULATION = "/simulation"
const GROUP_MESH = "/mesh"
const GROUP_CELLS = "/mesh/cells"
const GROUP_AUX = "/mesh/aux"
const GROUP_PARTS = "/mesh/parts"
const GROUP_NODES = "/mesh/nodes"
const GROUP_PARAMS = "/mesh/parameters"

# hdf name constants
const HDF_P = "p"
const HDF_RHO = "rho"
const HDF_GAMMA = "gamma"
const HDF_U = "u"
const HDF_V = "v"
const HDF_CELLPOROSITYTIMESCELLPOROSITY_FACTOR = "cellporositytimesporosityfactor"
const HDF_VISCOSITY = "viscosity"
const HDF_ITER = "iter"
const HDF_DELTAT = "deltat"
const HDF_T = "t"
const HDF_H_OUT = "h_out"
const HDF_POROSITY_OUT = "porosity_out"
const HDF_P_DM = "p_DM"
const HDF_RHO_DM = "rho_DM"
const HDF_GAMMA_DM = "gamma_DM"
const HDF_U_DM = "u_DM"
const HDF_V_DM = "v_DM"
const HDF_GAMMA_OUT = "gamma_out"
const HDF_W = "w"
const HDF_FILLED = "filled"
const HDF_FILLED_DM = "filled_DM"