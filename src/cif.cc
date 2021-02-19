#include <gemmi/cif.hpp>
#include <iostream>
#include <set>
#include <memory>

#include "cif.hh"

static std::unique_ptr<std::set<int>> 
get_models(const gemmi::cif::Document &doc)
{
    auto models = std::make_unique<std::set<int>> ();
    for (auto block : doc.blocks) {
        for (auto site : block.find("_atom_site.", {"pdbx_PDB_model_num"})) {
            models->insert(gemmi::cif::as_int(site[0]));
        }
    }
    return models;
}

static std::unique_ptr<std::set<std::string>>
get_chains(const gemmi::cif::Document &doc)
{
    auto chains = std::make_unique<std::set<std::string>> ();
    for (auto block : doc.blocks) {
        for (auto site : block.find("_atom_site.", {"auth_asym_id"})) {
            chains->insert(site[0]);
        }
    }
    return chains;
}

static const auto atom_site_columns = std::vector<std::string>({
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
    "pdbx_PDB_model_num",
});

static freesasa_structure *
structure_from_doc(const gemmi::cif::Document &doc,
                   const std::set<int> &models,
                   const freesasa_classifier *classifier,
                   int structure_options)
{
    freesasa_structure *structure = freesasa_structure_new();
    
    unsigned skip_count {0};
    for (auto block : doc.blocks) {
        auto prevAltId = '?';

        for (auto site : block.find("_atom_site.", atom_site_columns)) {

            // TODO figure out what this does
            if (site[0] != "ATOM" && !(structure_options & FREESASA_INCLUDE_HETATM)) {
                continue;
            }

            auto model = atoi(site[11].c_str());

            if (models.count(model) == 0) {
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

            auto currentAltId = site[6][0];

            if (!(structure_options & FREESASA_INCLUDE_HYDROGEN) && std::string(atom.type_symbol) == "H") {
                continue;
            }

            // Pick the first alternative conformation for an atom
            if (currentAltId != '.' && currentAltId != 'A'){
                skip_count++;
                std::cout << "Skipped Atom: ";
                for (auto idx : site) {std::cout << idx.c_str() << " ";}
                std::cout << std::endl;
                continue;
            }
            // if ((currentAltId != '.' && prevAltId == '.') || (currentAltId == '.')) {
            //     prevAltId = currentAltId;
            // } else if (currentAltId != '.' && currentAltId != prevAltId) {
            //     skip_count++;
            //     std::cout << "Skipped Atom: ";
            //     for (auto idx : site) {std::cout << idx.c_str() << " ";}
            //     std::cout << std::endl;
            //     continue;
            // }

            freesasa_structure_add_cif_atom(structure, &atom, classifier, structure_options);
        }
    }
    std::cout << "HETATM? " << (structure_options & FREESASA_INCLUDE_HETATM) << std::endl;
    std::cout << "Skipped " << skip_count << " atoms from cif due to Alt ID." << std::endl;
    return structure;
}

freesasa_structure *
freesasa_structure_from_cif(std::FILE *input,
                            const freesasa_classifier *classifier,
                            int structure_options)
{
    const auto doc = gemmi::cif::read_cstream(input, 8192, "cif-input");
    const auto models = get_models(doc);
    freesasa_structure *structure;

    if (structure_options & FREESASA_JOIN_MODELS) {
        structure = structure_from_doc(doc, *models, classifier, structure_options);
    } else {
        auto firstModel = models->begin();
        auto singleModel = std::set<int>();
        singleModel.insert(*firstModel);

        structure = structure_from_doc(doc, singleModel, classifier, structure_options);
    }

    return structure;
}
