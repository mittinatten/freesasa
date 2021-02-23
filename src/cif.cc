#include <iostream>
#include <set>
#include <memory>

#include <gemmi/cif.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmcif.hpp>

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

static std::unique_ptr<std::set<std::string>>
get_chains(const gemmi::Model& model)
{
    auto chains = std::make_unique<std::set<std::string>> ();

    for (auto& chain : model.chains)
    {
        chains->insert(chain.name);
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
    
    for (auto block : doc.blocks) {
        auto prevAltId = '?';

        for (auto site : block.find("_atom_site.", atom_site_columns)) {
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
            if (currentAltId != '.' && currentAltId != 'A') {
                continue;
            }

            freesasa_structure_add_cif_atom(structure, &atom, classifier, structure_options);
        }
    }
    return structure;
}

freesasa_structure *
freesasa_structure_from_cif(std::FILE *input,
                            const freesasa_classifier *classifier,
                            int structure_options)
{
    const auto doc = gemmi::cif::read_cstream(input, 8192, "cif-input");
    const auto models = get_models(doc);

    if (structure_options & FREESASA_JOIN_MODELS) {
        return structure_from_doc(doc, *models, classifier, structure_options);
    } else {
        auto firstModel = models->begin();
        auto singleModel = std::set<int>();
        singleModel.insert(*firstModel);

        return structure_from_doc(doc, singleModel, classifier, structure_options);
    }
}


std::vector<freesasa_structure*>
get_array_of_all_models(const gemmi::cif::Document &doc, int *n, const freesasa_classifier *classifier)
{
    const auto models = get_models(doc);

    std::vector<freesasa_structure*> ss;

    for (int i=0; i < models->size(); ++i)
    {
        ss.emplace_back(freesasa_structure_new());
    }  
    return ss; 
}


freesasa_structure*
chain_structure_from_doc(const gemmi::cif::Document doc, 
                         const std::string &chName,
                         const freesasa_classifier *classifier,
                         int structure_options)
{
    std::cout << "Creating a freesasa structure from chain!" << std::endl;
    freesasa_structure *structure = freesasa_structure_new();

    //TODO fill in structure with chain data

    return structure;
}

freesasa_structure*
model_structure_from_doc(const gemmi::cif::Document &doc,
                         const int model, 
                         const freesasa_classifier *classifier,
                         int structure_options)
{
    std::cout << "Creating a freesasa structure from model!" << std::endl;
    freesasa_structure *structure = freesasa_structure_new();

    //TODO fill in structure with model data

    return structure;
}

std::vector<freesasa_structure*>
freesasa_cif_structure_array(std::FILE *input,
                         int *n,
                         const freesasa_classifier *classifier,
                         int options)
{
    int n_models = 0, n_chains = 0;

    std::vector<freesasa_structure*> ss;

    const auto doc = gemmi::cif::read_cstream(input, 8192, "cif-input");

    gemmi::Structure gemmi_struct = gemmi::make_structure_from_block(doc.blocks[0]);

    const auto models = gemmi_struct.models;

    n_models = models.size();

    std::cout << "Number of models: " << n_models << std::endl;

    /* only keep first model if option not provided */
    if (!(options & FREESASA_SEPARATE_MODELS)) n_models = 1;

    /* for each model read chains if requested */
    if (options & FREESASA_SEPARATE_CHAINS) 
    {
        for (auto& model : models) 
        {
            auto chain_names  = get_chains(model);
            int n_new_chains  = chain_names->size();
            n_chains         +=  n_new_chains;

            if (n_new_chains == FREESASA_FAIL) gemmi::fail("No chains in protein");
            if (n_new_chains == 0) {
                // TODO Cant get this to link with freesasa.a for some reason.
                //freesasa_warn("in %s(): no chains found (in model %s)", __func__, model.name.c_str());
                continue;
            }

            ss.reserve(n_new_chains);
            for (auto& chain_name : *chain_names)
            {
                ss.emplace_back(
                    // TODO add model to freesasa_structure: implement this line from structure.c (ss[j0 + j]->model = i + 1;)
                    chain_structure_from_doc(doc, chain_name, classifier, options)
                );
            }
        }
        *n = n_chains;
    }

    else 
    {
        ss.reserve(n_models);
        for (auto& model : models) 
        {
            ss.emplace_back(
                // TODO add model to freesasa_structure: implement this line from structure.c (ss[j0 + j]->model = i + 1;)
                model_structure_from_doc(doc, std::stoi(model.name), classifier, options)
            );
        }
        *n = n_models;
    }
    return ss;
}
