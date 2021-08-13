module BWALGS
    filepath = "BWansatz/"

    include(join([filepath, "BWDMALGS.jl"]))
    export BW
end

using .BWALGS
