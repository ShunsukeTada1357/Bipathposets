module Bipathposets

export interval_decomposition
export interval_decomposition_cubical
export plot_bipath_diagram
export bipath_records_cubical
export interactive_bipath_viewer

#Varidation
include("IsFSC.jl")
#Core
include("FSCUtils.jl")
include("reduction_and_basechange_sparse.jl")
include("AtypeMethod.jl")
include("BipathMatMethod.jl")
#Image / cubical preprocessing
include("Images2FSC.jl")
#Visualization
include("BipathPlaneUtils.jl")
include("BipathPD.jl")

include("BipathRecords.jl")
include("ImageBipathGUI.jl")
end # module Bipathposets
