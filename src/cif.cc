#include "cif.hh"
#include <gemmi/cif.hpp>
#include <iostream>

freesasa_structure *
freesasa_structure_from_cif(std::FILE *input,
                            const freesasa_classifier *classifier,
                            int structure_options)
{
    auto doc = gemmi::cif::read_cstream(input, 8192, "cif-input");
    auto columns = std::vector<std::string>({
        "group_PDB",
        "auth_asym_id",
        "auth_seq_id",
        "pdbx_PDB_ins_code",
        "auth_comp_id",
        "auth_atom_id",
        "label_alt_id",
        "type_symbol",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
    });

    freesasa_structure *structure = freesasa_structure_new();
    for (auto block : doc.blocks) {
        for (auto site : block.find("_atom_site.", columns)) {
            if (site[0] != "ATOM" && !(structure_options & FREESASA_INCLUDE_HETATM)) {
                continue;
            }

            freesasa_cif_atom atom = {
                .group_PDB = site[0].c_str(),
                .auth_asym_id = site[1][0],
                .auth_seq_id = site[2].c_str(),
                .pdbx_PDB_ins_code = site[3].c_str(),
                .auth_comp_id = site[4].c_str(),
                .auth_atom_id = site[5].c_str(),
                .label_alt_id = site[6].c_str(),
                .type_symbol = site[7].c_str(),
                .Cartn_x = atof(site[8].c_str()),
                .Cartn_y = atof(site[9].c_str()),
                .Cartn_z = atof(site[10].c_str())};

            freesasa_structure_add_cif_atom(structure, &atom, classifier, structure_options);
        }
    }

    return structure;
}
