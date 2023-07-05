import h_transport_materials as htm


# These are the materials we'll plot
materials = [
    htm.TUNGSTEN,
    htm.ALUMINIUM,
    htm.INCONEL_600,
    htm.TITANIUM,
    htm.VANADIUM,
    htm.IRON,
    htm.NICKEL,
    htm.COPPER,
    htm.CUCRZR,
    htm.ZIRCONIUM,
    htm.MOLYBDENUM,
    htm.NIOBIUM,
    htm.PALLADIUM,
    htm.SILVER,
    htm.TANTALUM,
    htm.GOLD,
    htm.STEEL_RAFM,
    htm.STEEL_SERIES_300,
    htm.STEEL_316L,
    htm.INCONEL_625,
    htm.INCOLOY_800,
    htm.V4CR4TI,
    htm.PDAG,
    htm.LIPB,
]

# create 3 empty groups of properties
diffusivities, solubilities, permeabilities = (
    htm.PropertiesGroup(),
    htm.PropertiesGroup(),
    htm.PropertiesGroup(),
)

# iterate through the materials
for mat in materials:
    # get the first diffusivity and add it to the group
    diff = htm.diffusivities.filter(material=mat)[0]
    diffusivities.append(diff)

    # get the first solubility and add it to the group
    sol = htm.solubilities.filter(material=mat)[0]
    solubilities.append(sol)

    # some materials don't have a permeability in the database

    # filter the permeabilities for this material
    filtered_perms = htm.permeabilities.filter(material=mat)
    if len(filtered_perms) == 0:  # if no property was found

        # create a permeability from the product of diffusivity and solubility
        perm = htm.Permeability(
            pre_exp=diff.pre_exp * sol.pre_exp,
            act_energy=diff.act_energy + sol.act_energy,
            range=diff.range,
            material=mat,
        )
    else:
        # else, use the first entry
        perm = filtered_perms[0]
    
    # if the range is None (see issue #227 of HTM) us that of diff
    if perm.range is None: 
        perm.range = diff.range

    # add the permeability to the group
    permeabilities.append(perm)
