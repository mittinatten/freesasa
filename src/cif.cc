#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>

#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>

#define GEMMI_WRITE_IMPLEMENTATION
#include <gemmi/to_cif.hpp>
#undef GEMMI_WRITE_IMPLEMENTATION

#include "cif.hh"
#include "freesasa.h"
#include "freesasa_internal.h"

static std::map<size_t, gemmi::cif::Document> docs;

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

/**
 * Returns atom where .auth_atom_id is dynamically allocated and needs to be freed.
 *
 * TODO: Better solution needed. A previous version with std::move didn't work as expected,
 * and sometimes caused seg-faults.
 */
static freesasa_cif_atom
freesasa_atom_from_site(const gemmi::cif::Table::Row &site)
{
    const char *auth_atom_id;
    // remove quotation marks if necessary
    if (site[5][0] == '"') {
        auth_atom_id = std::string(site[5].substr(1, site[5].size() - 2)).c_str();
    } else {
        auth_atom_id = site[5].c_str();
    }

    return {
        .group_PDB = site[0].c_str(),
        .auth_asym_id = site[1][0],
        .auth_seq_id = site[2].c_str(),
        .pdbx_PDB_ins_code = site[3].c_str(),
        .auth_comp_id = site[4].c_str(),
        .auth_atom_id = strdup(auth_atom_id),
        .label_alt_id = site[6].c_str(),
        .type_symbol = site[7].c_str(),
        .Cartn_x = atof(site[8].c_str()),
        .Cartn_y = atof(site[9].c_str()),
        .Cartn_z = atof(site[10].c_str())};
}

template <typename T>
static freesasa_structure *
structure_from_pred(const gemmi::cif::Document &doc,
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

            // since this is in the interface between C and C++ code, some hackery is needed
            free((void *)atom.auth_atom_id);
        }
    }
    return structure;
}

static std::pair<gemmi::cif::Document &, size_t>
generate_gemmi_doc(std::FILE *input)
{
    static size_t index = 1;

    // theoretically there could be thread safety issues here, not sure they're relevant.
    size_t my_idx = index++;

    docs.emplace(my_idx, gemmi::cif::read_cstream(input, 8192, "cif-input"));
    return std::make_pair(std::ref(docs[my_idx]), my_idx);
}

static void release_gemmi_doc(size_t doc_ref)
{
    docs.erase(doc_ref);
}

freesasa_structure *
freesasa_structure_from_cif(std::FILE *input,
                            const freesasa_classifier *classifier,
                            int structure_options)
{
    auto doc_idx_pair = generate_gemmi_doc(input);

    const auto models = get_models(doc_idx_pair.first);

    std::unique_ptr<const ModelSetDiscriminator> discriminator;
    if (structure_options & FREESASA_JOIN_MODELS) {
        discriminator = std::make_unique<const ModelSetDiscriminator>(std::move(*models));
    } else {
        auto firstModel = models->begin();
        auto singleModel = std::set<int>{*firstModel};
        discriminator = std::make_unique<const ModelSetDiscriminator>(singleModel);
    }

    auto structure = structure_from_pred(doc_idx_pair.first, *discriminator, classifier, structure_options);
    freesasa_structure_set_cif_ref(structure, doc_idx_pair.second, &release_gemmi_doc);

    return structure;
}

static freesasa_structure *
structure_from_model(const gemmi::cif::Document &doc,
                     const std::string &model_name,
                     const freesasa_classifier *classifier,
                     int structure_options)
{
    const ModelDiscriminator discriminator(model_name);
    return structure_from_pred(doc, discriminator, classifier, structure_options);
}

static freesasa_structure *
structure_from_chain(const gemmi::cif::Document doc,
                     const std::string &model_name,
                     const std::string &chain_name,
                     const freesasa_classifier *classifier,
                     int structure_options)
{
    const ChainDiscriminator discriminator(model_name, chain_name);
    return structure_from_pred(doc, discriminator, classifier, structure_options);
}

std::vector<freesasa_structure *>
freesasa_cif_structure_array(std::FILE *input,
                             int *n,
                             const freesasa_classifier *classifier,
                             int options)
{
    int n_models = 0, n_chains = 0;

    std::vector<freesasa_structure *> ss;

    auto doc_idx_pair = generate_gemmi_doc(input);
    auto &doc = doc_idx_pair.first;

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
                freesasa_warn("in %s(): no chains found (in model %s)",
                              __func__, models[i].name.c_str());
                continue;
            }

            ss.reserve(n_new_chains);
            for (auto &chain_name : *chain_names) {
                freesasa_structure *structure = structure_from_chain(doc, models[i].name, chain_name, classifier, options);
                if (freesasa_structure_n(structure) == 0) {
                    --n_chains;
                    free(structure);
                    continue;
                }

                ss.push_back(std::move(structure));

                freesasa_structure_set_model(ss.back(), i + 1);
                freesasa_structure_set_cif_ref(ss.back(), doc_idx_pair.second, NULL);
            }
        }
        if (n_chains == 0)
            freesasa_fail("In %s(): No chains in any model in protein: %s.",
                          __func__, gemmi_struct.name.c_str());
        *n = n_chains;
    } else {
        ss.reserve(n_models);
        for (int i = 0; i < n_models; ++i) {
            ss.emplace_back(
                structure_from_model(doc, models[i].name, classifier, options));
            freesasa_structure_set_model(ss.back(), i + 1);
            freesasa_structure_set_cif_ref(ss.back(), doc_idx_pair.second, NULL);
        }
        *n = n_models;
    }

    // this is a hack, we only want to release docs once per input
    freesasa_structure_set_cif_ref(ss.back(), doc_idx_pair.second, &release_gemmi_doc);
    return ss;
}

struct freesasa_MCRA {
    freesasa_MCRA(const int model,
                  const std::string &chain,
                  const std::string &res_num,
                  const std::string &residue,
                  const std::string &atom)
        : _model(model), _chain(chain), _res_num(res_num), _residue(residue), _atom(atom)
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

        freesasa_warn(
            "In %s(), unable to find row in _atom_site for atom (%d, %s, %s, %s, %s). "
            "Looping through entire table to double check...",
            __func__, this->_model, this->_chain.c_str(),
            this->_res_num.c_str(), this->_residue.c_str(),
            this->_atom.c_str());

        idx = 0;
        for (const auto site : table) {
            if (this->is_row(site)) {
                return idx;
            }
            ++idx;
        }
        return FREESASA_FAIL;
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

static void
append_freesasa_params_to_block(gemmi::cif::Block &block, freesasa_node *result)
{
    assert(freesasa_node_type(result) == FREESASA_NODE_RESULT);

    const freesasa_parameters *params = freesasa_node_result_parameters(result);

    std::string version{PACKAGE_STRING};
    version = version.substr(version.find(" "));

    if (params == NULL) params = &freesasa_default_parameters;

    std::string params_prefix{"_freeSASA_parameters."};

    std::vector<std::string> params_tags{"version", "algorithm", "probe-radius"};
    std::vector<std::string> params_data{
        version,
        std::string{freesasa_alg_name(params->alg)},
        std::to_string(params->probe_radius)};

    switch (params->alg) {
    case FREESASA_SHRAKE_RUPLEY:
        params_tags.emplace_back("testpoints");
        params_data.emplace_back(std::to_string(params->shrake_rupley_n_points));
        break;
    case FREESASA_LEE_RICHARDS:
        params_tags.emplace_back("slices");
        params_data.emplace_back(std::to_string(params->lee_richards_n_slices));
        break;
    default:
        assert(0);
        break;
    }

    for (int i = 0; i != params_tags.size(); ++i) {
        // place quotes around the algorithm tag value
        if (i == 1)
            block.set_pair(params_prefix + params_tags[i], gemmi::cif::quote(params_data[i]));
        else
            block.set_pair(params_prefix + params_tags[i], params_data[i]);
    }
}

static std::string
correct_inf_nan_values(const double value)
{
    if (std::isnan(value)) return ".";
    if (std::isinf(value)) return "?";
    return std::to_string(value);
}

static void
append_freesasa_rsa_residue_to_block(gemmi::cif::Block &block, freesasa_node *residue)
{
    assert(freesasa_node_type(residue) == FREESASA_NODE_RESIDUE);

    const freesasa_nodearea *abs, *reference;
    freesasa_nodearea rel;

    abs = freesasa_node_area(residue);
    reference = freesasa_node_residue_reference(residue);

    if (reference)
        freesasa_residue_rel_nodearea(&rel, abs, reference);
    else
        rel = freesasa_nodearea_null;

    std::string rsa_prefix{"_freeSASA_rsa."};
    std::vector<std::string> rsa_tags{
        "asym_id",
        "seq_id",
        "comp_id",
        "abs_total",
        "rel_total",
        "abs_side_chain",
        "rel_side_chain",
        "abs_main_chain",
        "rel_main_chain",
        "abs_apolar",
        "rel_apolar",
        "abs_polar",
        "rel_polar",
    };

    std::vector<std::string> rsa_data{
        std::string{freesasa_node_name(freesasa_node_parent(residue))[0]},
        std::string{freesasa_node_residue_number(residue)},
        std::string{freesasa_node_name(residue)},
        correct_inf_nan_values(abs->total),
        correct_inf_nan_values(rel.total),
        correct_inf_nan_values(abs->side_chain),
        correct_inf_nan_values(rel.side_chain),
        correct_inf_nan_values(abs->main_chain),
        correct_inf_nan_values(rel.main_chain),
        correct_inf_nan_values(abs->apolar),
        correct_inf_nan_values(rel.apolar),
        correct_inf_nan_values(abs->polar),
        correct_inf_nan_values(rel.polar),
    };

    gemmi::cif::Loop *rsa_loop;
    if (block.find(rsa_prefix, rsa_tags).ok()) {
        rsa_loop = block.find(rsa_prefix, rsa_tags).get_loop();
    } else {
        gemmi::cif::Loop &temp_loop = block.init_loop(rsa_prefix, rsa_tags);
        rsa_loop = &temp_loop;
    }

    rsa_loop->add_row(rsa_data);
}

static int
append_freesasa_result_summary_to_block(gemmi::cif::Block &block, freesasa_node *result)
{
    assert(freesasa_node_type(result) == FREESASA_NODE_RESULT);

    freesasa_node *structure = NULL, *chain = NULL;
    const freesasa_nodearea *area = NULL;

    structure = freesasa_node_children(result);
    if (structure == NULL) {
        return freesasa_fail("Result node has no structure nodes.");
    }

    area = freesasa_node_area(structure);
    if (area == NULL) {
        return freesasa_fail("Structure node has no area.");
    }

    std::string results_prefix{"_freeSASA_results."};
    std::vector<std::string> result_tags{"model", "chains", "atoms", "type", "surface_area"};
    std::vector<std::string> template_data(result_tags.size());
    std::vector<std::vector<std::string>> result_data;

    template_data[0] = std::to_string(freesasa_node_structure_model(structure));
    template_data[1] = std::string{freesasa_node_structure_chain_labels(structure)};
    template_data[2] = std::to_string(freesasa_node_structure_n_atoms(structure));

    template_data[3] = "Total";
    template_data[4] = std::to_string(area->total);
    result_data.push_back(template_data);

    template_data[3] = "Apolar";
    template_data[4] = std::to_string(area->apolar);
    result_data.push_back(template_data);

    template_data[3] = "Polar";
    template_data[4] = std::to_string(area->polar);
    result_data.push_back(template_data);

    if (area->unknown > 0) {
        template_data[3] = "Unknown";
        template_data[4] = std::to_string(area->unknown);
        result_data.push_back(template_data);
    }

    chain = freesasa_node_children(structure);
    while (chain) {
        area = freesasa_node_area(chain);
        if (area == NULL) {
            return freesasa_fail("Chain has no stored area");
        }

        template_data[3] = gemmi::cif::quote(std::string{"CHAIN "} + std::string{freesasa_node_name(chain)});
        template_data[4] = std::to_string(area->total);
        result_data.push_back(template_data);

        chain = freesasa_node_next(chain);
    }

    gemmi::cif::Loop *result_loop;
    if (block.find(results_prefix, result_tags).ok()) {
        result_loop = block.find(results_prefix, result_tags).get_loop();
    } else {
        gemmi::cif::Loop &temp_loop = block.init_loop(results_prefix, result_tags);
        result_loop = &temp_loop;
    }

    for (auto &row : result_data) {
        result_loop->add_row(row);
    }

    return FREESASA_SUCCESS;
}

static void
populate_freesasa_result_vectors(gemmi::cif::Table &table, freesasa_node *result,
                                 std::vector<std::string> &sasa_vals,
                                 std::vector<std::string> &sasa_radii)
{
    assert(freesasa_node_type(result) == FREESASA_NODE_RESULT);

    assert(table.ok());

    freesasa_node *structure, *chain, *residue, *atom;
    int rowNum{0}, model{0};

    structure = freesasa_node_children(result);
    while (structure) {
        rowNum = 0;
        model = freesasa_node_structure_model(structure);
        chain = freesasa_node_children(structure);
        while (chain) {
            residue = freesasa_node_children(chain);
            while (residue) {
                append_freesasa_rsa_residue_to_block(table.bloc, residue);
                atom = freesasa_node_children(residue);
                while (atom) {
                    auto cName = std::string(1, freesasa_node_atom_chain(atom));
                    // TODO figure out why this returns decimal string sometimes
                    auto rNum = std::to_string(std::atoi(freesasa_node_atom_residue_number(atom)));
                    auto rName = freesasa_node_atom_residue_name(atom);
                    auto aName = freesasa_node_name(atom);
                    auto area = freesasa_node_area(atom);
                    auto radius = freesasa_node_atom_radius(atom);

                    rowNum = freesasa_MCRA{model, cName, rNum, rName, aName}.find_row(table, rowNum);
                    if (rowNum == FREESASA_FAIL)
                        freesasa_fail(
                            "In %s(), unable to find freesasa_node atom (%d, %s, %s, %s, %s) in cif %s",
                            __func__, model, cName.c_str(), rNum.c_str(), rName, aName, table.bloc.name.c_str());

                    sasa_vals[rowNum] = std::to_string(area->total);
                    sasa_radii[rowNum] = std::to_string(radius);

                    atom = freesasa_node_next(atom);
                }
                residue = freesasa_node_next(residue);
            }
            chain = freesasa_node_next(chain);
        }
        structure = freesasa_node_next(structure);
    }
}

static int
replace_sasa_columns(gemmi::cif::Table &table,
                     const std::vector<std::string> &sasa_vals,
                     const std::vector<std::string> &sasa_radii,
                     const std::vector<std::string> &sasa_tags)
{
    std::vector<std::vector<std::string>> sasa_columns = {sasa_vals, sasa_radii};

    assert(sasa_columns.size() == sasa_tags.size());

    for (unsigned i = 0; i != sasa_tags.size(); ++i) {
        auto column = table.bloc.find_loop(sasa_tags[i]);
        std::copy(sasa_columns[i].begin(), sasa_columns[i].end(), column.begin());
    }

    return FREESASA_SUCCESS;
}

static int
append_sasa_columns(gemmi::cif::Table &table,
                    const std::vector<std::string> &sasa_vals,
                    const std::vector<std::string> &sasa_radii,
                    const std::vector<std::string> &sasa_tags)
{
    auto &loop = *table.get_loop();

    unsigned long orig_tag_size = loop.tags.size();
    unsigned long new_tag_size = orig_tag_size + 2;

    // Creates a new table full of empty strings with the correct number of dimensions
    // Outside vector size is the # of columns, inside vector size is the # of rows.
    std::vector<std::vector<std::string>> new_columns(new_tag_size, {loop.length(), {"Empty"}});

    // Copies data from original columns to their respecitve column in the new table filled with empty strings.
    // Leaving only the new appended columns as empty strings
    for (unsigned i = 0; i != orig_tag_size; ++i) {
        auto column = table.bloc.find_loop(loop.tags[i]);
        std::copy(column.begin(), column.end(), new_columns[i].begin());
    }

    new_columns[new_tag_size - 2] = std::move(sasa_vals);
    new_columns[new_tag_size - 1] = std::move(sasa_radii);

    for (const auto &tag : sasa_tags)
        loop.tags.push_back(tag);

    loop.set_all_values(new_columns);

    return FREESASA_SUCCESS;
}

static int
rewrite_atom_site(gemmi::cif::Table &table,
                  std::vector<std::string> &sasa_vals,
                  std::vector<std::string> &sasa_radii)
{
    std::vector<std::string> sasa_tags{
        "_atom_site.FreeSASA_value",
        "_atom_site.FreeSASA_radius"};

    auto &loop = *table.get_loop();

    if (loop.has_tag(sasa_tags[0]) && loop.has_tag(sasa_tags[1])) {
        return replace_sasa_columns(table, sasa_vals, sasa_radii, sasa_tags);
    } else {
        return append_sasa_columns(table, sasa_vals, sasa_radii, sasa_tags);
    }
}

static int
reset_freesasa_tables(gemmi::cif::Block &block)
{
    std::vector<gemmi::cif::Table> result_tables{
        block.find_mmcif_category("_freeSASA_results."),
        block.find_mmcif_category("_freeSASA_rsa."),
        block.find_mmcif_category("_freeSASA_parameters.")};

    for (auto &table : result_tables) {
        if (table.ok()) {
            // Table exist. Making sure its a loop.
            if (!table.loop_item) {
                // Table is a pair so turning it into a loop
                table.convert_pair_to_loop();
            }
            table.get_loop()->clear();
        }
    }

    return FREESASA_SUCCESS;
}

static int
write_result(std::ostream &out, freesasa_node *root)
{
    freesasa_node *result{freesasa_node_children(root)};

    size_t cif_ref, prev_cif_ref = 0;
    bool write = false;
    std::vector<std::string> sasa_vals, sasa_radii;

    while (result) {
        cif_ref = freesasa_node_structure_cif_ref(freesasa_node_children(result));
        if (docs.find(cif_ref) == docs.end()) {
            return freesasa_fail("In %s(), unable to find gemmi doc for result node: %s."
                                 "This can happen when using the --separate-chains option",
                                 __func__, freesasa_node_name(result));
        }

        auto &block = docs[cif_ref].sole_block();

        if (cif_ref != prev_cif_ref) {
            reset_freesasa_tables(block);
        }

        auto table = block.find("_atom_site.", atom_site_columns);
        if (prev_cif_ref != cif_ref) {
            sasa_vals = std::vector<std::string>{table.length(), "?"};
            sasa_radii = std::vector<std::string>{table.length(), "?"};
        }

        populate_freesasa_result_vectors(table, result, sasa_vals, sasa_radii);

        int added_data = rewrite_atom_site(table, sasa_vals, sasa_radii);
        int added_summary = append_freesasa_result_summary_to_block(block, result);

        if (added_data == FREESASA_FAIL || added_summary == FREESASA_FAIL) {
            return freesasa_fail("Unable to build CIF output");
        }

        append_freesasa_params_to_block(block, result);

        prev_cif_ref = cif_ref;
        result = freesasa_node_next(result);

        if (!result) {
            // There is no next result so write out current file.
            write = true;
        } else if (freesasa_node_structure_cif_ref(freesasa_node_children(result)) != prev_cif_ref) {
            // Next result node is from a new file so write out current file.
            write = true;
        } else {
            // Next result node is from the same doc so do not write current file yet.
            write = false;
        }

        if (write) gemmi::cif::write_cif_block_to_stream(out, block);
    }
    return FREESASA_SUCCESS;
}

int freesasa_export_tree_to_cif(const char *filename,
                                freesasa_node *root)
{
    assert(root);
    assert(freesasa_node_type(root) == FREESASA_NODE_ROOT);
    int ret;

    std::streambuf *buf;
    std::ofstream fout;

    if (filename != NULL) {
        fout.open(filename);
        buf = fout.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }

    std::ostream out(buf);

    try {
        ret = write_result(out, root);
    } catch (...) {
        ret = FREESASA_FAIL;
    }

    if (filename != NULL) {
        fout.close();
    }

    if (ret == FREESASA_FAIL) {
        freesasa_fail("Unable to output CIF file");
    }

    return ret;
}
