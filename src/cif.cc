#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <algorithm>

#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>

#define GEMMI_WRITE_IMPLEMENTATION
#include <gemmi/to_cif.hpp>
#undef GEMMI_WRITE_IMPLEMENTATION

#include "cif.hh"
#include "freesasa.h"


static std::vector<gemmi::cif::Document> docs;

struct ModelDiscriminator {
    ModelDiscriminator(const std::string &model_name,
                       const int model_col = 11)
        : _model_name(model_name), _model_col(model_col)
    {
    }

    bool operator()(const gemmi::cif::Table::Row &site) const
    {
        return _model_name != site[_model_col];
    }

private:
    const std::string _model_name;
    int _model_col;
};

struct ModelSetDiscriminator {
    ModelSetDiscriminator(const std::set<int> models,
                          const int model_col = 11)
        : _models(models), _model_col(model_col)
    {
    }

    bool operator()(const gemmi::cif::Table::Row &site) const
    {
        return _models.count(std::stoi(site[_model_col])) == 0;
    }

private:
    const std::set<int> _models;
    int _model_col;
};

struct ChainDiscriminator {
    ChainDiscriminator(const std::string &model_name, const std::string &chain_name,
                       const int model_col = 11, const int chain_col = 1)
        : _model_name(model_name), _chain_name(chain_name),
          _model_col(model_col), _chain_col(chain_col)
    {
    }

    bool operator()(const gemmi::cif::Table::Row &site) const
    {
        return _model_name != site[_model_col] || _chain_name != site[_chain_col];
    }

private:
    const std::string _model_name, _chain_name;
    int _model_col, _chain_col;
};

static std::unique_ptr<std::set<int>>
get_models(const gemmi::cif::Document &doc)
{
    auto models = std::make_unique<std::set<int>>();
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
    auto chains = std::make_unique<std::set<std::string>>();
    for (auto block : doc.blocks) {
        for (auto site : block.find("_atom_site.", {"auth_asym_id"})) {
            chains->insert(site[0]);
        }
    }
    return chains;
}

static std::unique_ptr<std::set<std::string>>
get_chains(const gemmi::Model &model)
{
    auto chains = std::make_unique<std::set<std::string>>();

    for (auto &chain : model.chains) {
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

static freesasa_cif_atom
freesasa_atom_from_site(const gemmi::cif::Table::Row &site)
{

    std::unique_ptr<std::string> auth_atom_id;
    // remove quotation marks if necessary
    if (site[5][0] == '"') {
        auth_atom_id = std::make_unique<std::string>(site[5].substr(1, site[5].size() - 2));
    }

    return {
        .group_PDB = site[0].c_str(),
        .auth_asym_id = site[1][0],
        .auth_seq_id = site[2].c_str(),
        .pdbx_PDB_ins_code = site[3].c_str(),
        .auth_comp_id = site[4].c_str(),
        .auth_atom_id = auth_atom_id ? std::move(*auth_atom_id).c_str() : site[5].c_str(),
        .label_alt_id = site[6].c_str(),
        .type_symbol = site[7].c_str(),
        .Cartn_x = atof(site[8].c_str()),
        .Cartn_y = atof(site[9].c_str()),
        .Cartn_z = atof(site[10].c_str())};
}

template <typename T>
static freesasa_structure *
freesasa_structure_from_pred(const gemmi::cif::Document &doc,
                             const T &discriminator,
                             const freesasa_classifier *classifier,
                             int structure_options)
{
    freesasa_structure *structure = freesasa_structure_new();
    std::string auth_atom_id;
    char prevAltId = '.';

    for (auto block : doc.blocks) {
        for (auto site : block.find("_atom_site.", atom_site_columns)) {
            if (site[0] != "ATOM" && !(structure_options & FREESASA_INCLUDE_HETATM)) {
                continue;
            }

            if (discriminator(site)) continue;

            freesasa_cif_atom atom = freesasa_atom_from_site(site);

            if (!(structure_options & FREESASA_INCLUDE_HYDROGEN) && std::string(atom.type_symbol) == "H") {
                continue;
            }

            // Pick the first alternative conformation for an atom
            auto currentAltId = site[6][0];
            if ((currentAltId != '.' && prevAltId == '.') || currentAltId == '.') {
                prevAltId = currentAltId;
            } else if (currentAltId != '.' && currentAltId != prevAltId) {
                continue;
            }

            freesasa_structure_add_cif_atom(structure, &atom, classifier, structure_options);
        }
    }
    return structure;
}


static gemmi::cif::Document&
generate_gemmi_doc(std::FILE * input)
{
    docs.emplace_back(gemmi::cif::read_cstream(input, 8192, "cif-input"));
    auto &doc = docs.back();

    gemmi::Structure gemmi_struct = gemmi::make_structure_from_block(doc.blocks[0]);

    if (gemmi_struct.name.find(".cif") != std::string::npos) {
        doc.source = gemmi_struct.name;
    }  else {
        doc.source = gemmi_struct.name + ".cif";
    }
    transform(doc.source.begin(), doc.source.end(), doc.source.begin(), tolower);

    return doc;
}

freesasa_structure *
freesasa_structure_from_cif(std::FILE *input,
                            const freesasa_classifier *classifier,
                            int structure_options)
{
    auto &doc = generate_gemmi_doc(input);

    std::cout << "Renamed doc.source to " << doc.source << std::endl;

    const auto models = get_models(doc);

    std::unique_ptr<const ModelSetDiscriminator> discriminator;
    if (structure_options & FREESASA_JOIN_MODELS) {
        discriminator = std::make_unique<const ModelSetDiscriminator>(std::move(*models));
    } else {
        auto firstModel = models->begin();
        auto singleModel = std::set<int>{*firstModel};
        discriminator = std::make_unique<const ModelSetDiscriminator>(singleModel);
    }
    return freesasa_structure_from_pred(doc, *discriminator, classifier, structure_options);
}

freesasa_structure *
freesasa_structure_from_model(const gemmi::cif::Document &doc,
                              const std::string &model_name,
                              const freesasa_classifier *classifier,
                              int structure_options)
{
    const ModelDiscriminator discriminator(model_name);
    return freesasa_structure_from_pred(doc, discriminator, classifier, structure_options);
    freesasa_structure *structure = freesasa_structure_new();
}

freesasa_structure *
freesasa_structure_from_chain(const gemmi::cif::Document doc,
                              const std::string &model_name,
                              const std::string &chain_name,
                              const freesasa_classifier *classifier,
                              int structure_options)
{
    const ChainDiscriminator discriminator(model_name, chain_name);
    return freesasa_structure_from_pred(doc, discriminator, classifier, structure_options);
}

std::vector<freesasa_structure *>
freesasa_cif_structure_array(std::FILE *input,
                             int *n,
                             const freesasa_classifier *classifier,
                             int options)
{
    int n_models = 0, n_chains = 0;

    std::vector<freesasa_structure *> ss;

    auto &doc = generate_gemmi_doc(input);

    gemmi::Structure gemmi_struct = gemmi::make_structure_from_block(doc.blocks[0]);

    const auto models = gemmi_struct.models;

    n_models = models.size();

    /* only keep first model if option not provided */
    if (!(options & FREESASA_SEPARATE_MODELS)) n_models = 1;

    /* for each model read chains if requested */
    if (options & FREESASA_SEPARATE_CHAINS) {
        for (int i = 0; i < n_models; ++i) {
            auto chain_names = get_chains(models[i]);
            int n_new_chains = chain_names->size();
            n_chains += n_new_chains;

            if (n_new_chains == 0) {
                freesasa_warn("in %s(): no chains found (in model %s)", __func__, models[i].name.c_str());
                continue;
            }

            ss.reserve(n_new_chains);
            for (auto &chain_name : *chain_names) {
                ss.emplace_back(
                    freesasa_structure_from_chain(doc, models[i].name, chain_name, classifier, options));
                freesasa_structure_set_model(ss.back(), i + 1);
            }
        }
        if (n_chains == 0) freesasa_fail("In %s(): No chains in any model in protein: %s.", __func__, gemmi_struct.name.c_str());
        *n = n_chains;
    } else {
        ss.reserve(n_models);
        for (int i = 0; i < n_models; ++i) {
            ss.emplace_back(
                freesasa_structure_from_model(doc, models[i].name, classifier, options));
            freesasa_structure_set_model(ss.back(), i + 1);
        }
        *n = n_models;
    }
    return ss;
}

struct freesasa_MCRA {

    freesasa_MCRA(const int model, const std::string &chain, const std::string &residue, const std::string &res_num, const std::string &atom)
        : _model(model), _chain(chain), _residue(residue), _res_num(res_num), _atom(atom)
    {
    }

    int find_row(gemmi::cif::Table &table, int start_idx = 0) const
    {
        int idx = start_idx, total_rows = table.length();

        assert(idx < total_rows);

        while (idx < total_rows) {
            if (this->is_row(table[idx])) {
                return idx;
            }
            ++idx;
        }

        std::cout << "Did not find row for: " << *this << std::endl;
        std::cout << "Looping through entire table... " << std::endl;

        idx = 0;
        for (const auto site : table) {
            if (this->is_row(site)) {
                return idx;
            }
            ++idx;
        }
        std::cout << "Still couldnt find: " << *this << std::endl;
        return -1;
    }

    friend std::ostream &operator<<(std::ostream &os, const freesasa_MCRA &mcra)
    {
        os << "Atom(" << mcra._model << " " << mcra._chain << " " <<  mcra._res_num << " [" << mcra._residue <<  "] " << mcra._atom << ")";
        return os;
    }

private:
    bool is_row(const gemmi::cif::Table::Row &row) const
    {
        int model = std::stoul(row[11]);
        const std::string &chain = row[1];
        const std::string &residue = row[4];
        const std::string &res_num = row[2];
        const std::string atom = row[5][0] != '"' ? row[5] : row[5].substr(1, row[5].size() - 2);

        if (_model == model) {
            if (_chain == chain) {
                if (_res_num == res_num && _residue == residue) {
                    if (_atom == atom) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    const int _model;
    const std::string &_chain, &_residue, &_res_num, &_atom;
};


static int find_doc_idx(std::string filename){

    std::transform(filename.begin(), filename.end(), filename.begin(), tolower);    
    for (int i = 0; i != docs.size(); ++i){
        if (docs[i].source == filename){
            return i;
        }
    }
    return -1;
}


static void 
append_freesasa_params_to_cif(gemmi::cif::Block& block, freesasa_node *result)
{
    assert(freesasa_node_type(result) == FREESASA_NODE_RESULT);

    const freesasa_parameters *params = freesasa_node_result_parameters(result);

    if (params == NULL) params = &freesasa_default_parameters;

    std::string params_prefix {"_freeSASA_parameters."};

    std::vector<std::string> params_tags {"algorithm", "probe-radius"};
    std::vector<std::string> params_data {
        std::string{freesasa_alg_name(params->alg)}, std::to_string(params->probe_radius)
    };

    #if USE_THREADS
        params_tags.push_back("threads");
        params_data.push_back(std::to_string(params->n_threads));
    #endif

    switch (params->alg) {
        case FREESASA_SHRAKE_RUPLEY:
            params_tags.push_back("testpoints");
            params_data.push_back(std::to_string(params->shrake_rupley_n_points));
            break;
        case FREESASA_LEE_RICHARDS:
            params_tags.push_back("slices");
            params_data.push_back(std::to_string(params->lee_richards_n_slices));
            break;
        default:
            assert(0);
            break;
    }

    for (int i = 0; i != params_tags.size(); ++i)
    {
        // place quotes around the algorithm tag value
        if (i == 0) block.set_pair(params_prefix + params_tags[i], gemmi::cif::quote(params_data[i]));
        else block.set_pair(params_prefix + params_tags[i], params_data[i]);
    }
}


static void 
append_freesasa_result_summary_to_cif(gemmi::cif::Block& block, freesasa_node *result)
{
    assert(freesasa_node_type(result) == FREESASA_NODE_RESULT);

    freesasa_node *structure = NULL, *chain = NULL;
    const freesasa_nodearea *area = NULL;

    std::string name{freesasa_node_name(result)};

    structure = freesasa_node_children(result);
    assert(structure);

    area = freesasa_node_area(structure);
    assert(area);

    std::string results_prefix{"_freeSASA_results."};
    std::vector<std::string> result_tags {"source", "chains", "model", "atoms", "type", "surface_area"};
    std::vector<std::string> template_data(result_tags.size()); 
    std::vector<std::vector<std::string>>result_data;

    if (name.empty()) template_data[0] = "unknown";
    else template_data[0] = name;

    template_data[1] = std::string{freesasa_node_structure_chain_labels(structure)};
    template_data[2] = std::to_string(freesasa_node_structure_model(structure));
    template_data[3] = std::to_string(freesasa_node_structure_n_atoms(structure));
    
    template_data[4] = "Total";
    template_data[5] = std::to_string(area->total);
    result_data.push_back(template_data);

    template_data[4] = "Apolar";
    template_data[5] = std::to_string(area->apolar);
    result_data.push_back(template_data);

    template_data[4] = "Polar";
    template_data[5] = std::to_string(area->polar);
    result_data.push_back(template_data);

    if (area->unknown > 0)
    {
        template_data[4] = "Unknown";
        template_data[5] = std::to_string(area->unknown);
        result_data.push_back(template_data);
    }

    chain = freesasa_node_children(structure);
    while (chain) 
    {
        area = freesasa_node_area(chain);
        assert(area);

        template_data[4] = gemmi::cif::quote(std::string{"CHAIN "} + std::string{freesasa_node_name(chain)});
        template_data[5] = std::to_string(area->total);
        result_data.push_back(template_data);

        chain = freesasa_node_next(chain);
    }

    gemmi::cif::Loop *result_loop;
    if (block.find(results_prefix, result_tags).ok())
    {
        result_loop = block.find(results_prefix, result_tags).get_loop();
    } 
    else 
    {
        gemmi::cif::Loop& temp_loop = block.init_loop(results_prefix, result_tags);
        result_loop = &temp_loop;
    }

    for (auto& row : result_data) 
    { 
        result_loop->add_row(row); 
        
    }
}


int freesasa_write_cif(std::FILE *output,
                       freesasa_node *root,
                       int options)
{
    freesasa_node *result{freesasa_node_children(root)};
    freesasa_node *structure, *chain, *residue, *atom;

    assert(output);
    assert(root);
    assert(freesasa_node_type(root) == FREESASA_NODE_ROOT);

    std::string prev_file = "";
    bool write = false;
    std::vector<std::string> sasa_vals, sasa_radii;
    while (result) {
        std::string inp_file{freesasa_node_name(result)};
        std::cout << "New Result with Input file: " << inp_file << std::endl;

        int idx = find_doc_idx(inp_file);
        if (idx == -1){
            std::cout << "Unable to find gemmi doc for result: " << inp_file << " Skipping..." << std::endl;
            result = freesasa_node_next(result);
            continue;
        } 
        auto &doc = docs[idx];
        std::cout << "Doc name: " << doc.source << std::endl;
        bool equal = doc.source == inp_file;
        std::cout << "file == doc? " << equal << std::endl;

        auto& block = doc.sole_block();
        append_freesasa_params_to_cif(block, result);
        append_freesasa_result_summary_to_cif(block, result);

        auto table = block.find("_atom_site.", atom_site_columns);
        if (prev_file != inp_file) {
            std::cout << "New file new vectors!" << std::endl;
            sasa_vals  = std::vector<std::string>{table.length(), "?"};
            sasa_radii = std::vector<std::string>{table.length(), "?"};
        } 
        structure = freesasa_node_children(result);
        while (structure) {
            int rowNum = 0, idx = 0;
            auto model = freesasa_node_structure_model(structure);
            chain = freesasa_node_children(structure);
            std::cout << "New Structure with model: " << model << std::endl;
            while (chain) {
                auto cName = freesasa_node_name(chain);
                residue = freesasa_node_children(chain);
                std::cout << "New Chain: " << cName << std::endl;
                while (residue) {
                    auto rName = freesasa_node_name(residue);
                    auto rNum = freesasa_node_residue_number(residue);
                    atom = freesasa_node_children(residue);
                    while (atom) {
                        auto aName = freesasa_node_name(atom);
                        auto area = freesasa_node_area(atom);
                        auto radius = freesasa_node_atom_radius(atom);

                        idx = freesasa_MCRA{model, cName, rName, rNum, aName}.find_row(table, rowNum);
                        if (idx != -1) {
                            rowNum = idx;
                            sasa_vals[rowNum] = std::to_string(area->total);
                            sasa_radii[rowNum] = std::to_string(radius);
                        } else {
                            exit(1);
                        }
                        atom = freesasa_node_next(atom);
                    }
                    auto last_res_name = freesasa_node_name(residue);
                    auto last_res_number = freesasa_node_residue_number(residue);
                    residue = freesasa_node_next(residue);
                }
                auto last_chain = freesasa_node_name(chain);
                chain = freesasa_node_next(chain);
                std::cout << "Finished chain: " << cName << std::endl;
            }
            structure = freesasa_node_next(structure);
        }
        prev_file = freesasa_node_name(result);
        result = freesasa_node_next(result);

        if (!result) {
            // There is no next result so write out current (last) file. 
            write = true;
        } else if (freesasa_node_name(result) != doc.source){
            // Next result node is from a new file so write out current file.
            write = true;
        } else {
            // Next result node is from the some doc so not ready for writing. 
            write = false;
        }
        
        if (write){
            std::cout << "Writing: " << doc.source << std::endl;

            auto &loop = *table.get_loop();

            unsigned long orig_tag_size = loop.tags.size();
            unsigned long new_tag_size = orig_tag_size + 2;

            // Creates a new table full of empty strings with the correct number of dimensions
            // Outside vector size is the # of columns, inside vector size is the # of rows. 
            std::vector<std::vector<std::string>> newCols(new_tag_size, {loop.length(), {"Empty"}});

            // Copies data from original columns to their respecitve column in the new table filled with empty strings.
            // Leaving only the new appended columns as empty strings
            for (unsigned int i = 0; i != orig_tag_size; ++i)
            {
                auto iCol = block.find_loop(loop.tags[i]);
                std::copy(iCol.begin(), iCol.end(), newCols[i].begin());
            }

            newCols[new_tag_size - 2] = std::move(sasa_vals);
            newCols[new_tag_size - 1] = std::move(sasa_radii);

            std::vector<std::string> new_tags{"_atom_site.FreeSASA_value", "_atom_site.FreeSASA_radius"};
            for (auto tag : new_tags) loop.tags.push_back(tag); 

            loop.set_all_values(newCols);

            if (output == stdout){
                gemmi::cif::write_cif_to_stream(std::cout, doc);
            } else {
                std::string out_file = inp_file.substr(0, inp_file.find_last_of(".")) + ".sasa.cif";
                std::ofstream newCif;
                newCif.open(out_file);
                gemmi::cif::write_cif_to_stream(newCif, doc);
                newCif.close();
            }
            unsigned count{0};
            for (auto& value : newCols[new_tag_size - 2]){
                if (value == "?") ++count;
            }
            std::cout << "Number of rows with no FREESASA value: " << count << std::endl;
            std::cout << "Number of rows x columns in the out file: " << loop.length() << " x " << loop.width() << std::endl;

        }
    }
    return FREESASA_SUCCESS;
}