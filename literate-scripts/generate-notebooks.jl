using GlacioTeaching
pushfirst!(LOAD_PATH, "..")

# this file-types are just copied over
asset_files = [".png"]

# 
for dir in readdir(joinpath(@__DIR__, "."), join=true)
    !isdir(dir) && continue
    cd(dir)
    for typ in [:assignment, :sol]
        outdir = if typ==:assignment
            joinpath(@__DIR__, "notebooks", splitdir(dir)[2])
        else
            joinpath(@__DIR__, "notebooks", splitdir(dir)[2]* "-solution")
        end
        GlacioTeaching.process_folder(dir, outdir, :nb;
                                      make_outputs=typ,
                                      execute=:sol,
                                      path_nbinclude=nothing, # hmm, probably remove this hack.  Related to https://github.com/stevengj/NBInclude.jl/issues/28
                                      asset_files
                        )
    end
end
# copy Project.toml
cp(joinpath(@__DIR__, "scripts/Project.toml"), joinpath(@__DIR__, "notebooks/Project.toml"), force=true)

# # make homework
# for dir in readdir(joinpath(@__DIR__, "scripts/homework"), join=true)
#     !isdir(dir) && continue
#     for typ in [:assignment, :sol]
#         outdir = if typ==:assignment
#             joinpath(@__DIR__, "notebooks", "homework", splitdir(dir)[2])
#         else
#             joinpath(@__DIR__, "notebooks", "homework-solution", splitdir(dir)[2])
#         end
#         GlacioTeaching.process_folder(dir, outdir, :nb;
#                                       make_outputs=typ,
#                                       execute=false #:sol
#                                       )
#     end
# end
