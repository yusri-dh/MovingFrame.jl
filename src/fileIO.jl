

load_mesh = load
save_mesh = save
function create_group_df(data_folder)
    fname1= joinpath(data_folder,"cor_gc_gda.csv")
    fname2= joinpath(data_folder,"gda.tmdf")

    cor_gc_gda=CSV.read(fname1)
    gda = CSV.read(fname2)[!,[:cellid,:trackid,:x,:y,:z,:t]]
    n_sample =size(gda,1)
    rename!(cor_gc_gda,[:cellid,:t,:gc_cellid])
    df = join(gda,cor_gc_gda,on=[:t,:cellid],kind=:inner)
    dropmissing!(df)
    df[!,:t] = convert.(Float64,df[!,:t])
    return df
end
