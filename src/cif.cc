#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
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
        }
    }
    return structure;
}

static gemmi::cif::Document &
generate_gemmi_doc(std::FILE *input)
{
    docs.emplace_back(gemmi::cif::read_cstream(input, 8192, "cif-input"));
    auto &doc = docs.back();

    gemmi::Structure gemmi_struct = gemmi::make_structure_from_block(doc.blocks[0]);

    if (gemmi_struct.name.find(".cif") != std::string::npos) {
        doc.source = gemmi_struct.name;
    } else {
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

    const auto models = get_models(doc);

    std::unique_ptr<const ModelSetDiscriminator> discriminator;
    if (structure_options & FREESASA_JOIN_MODELS) {
        discriminator = std::make_unique<const ModelSetDiscriminator>(std::move(*models));
    } else {
        auto firstModel = models->begin();
        auto singleModel = std::set<int>{*firstModel};
        discriminator = std::make_unique<const ModelSetDiscriminator>(singleModel);
    }
    return structure_from_pred(doc, *discriminator, classifier, structure_options);
}

static freesasa_structure *
structure_from_model(const gemmi::cif::Document &doc,
                     const std::string &model_name,
                     const freesasa_classifier *classifier,
                     int structure_options)
{
    const ModelDiscriminator discriminator(model_name);
    return structure_from_pred(doc, discriminator, classifier, structure_options);
    freesasa_structure *structure = freesasa_structure_new();
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
        }
        *n = n_models;
    }
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

static std::string
get_cif_filename(std::string filename)
{
    // Remove :<model number> if present in the cif filename
    if (filename.find(".cif") != std::string::npos)
        return filename.substr(0, filename.find(".cif") + 4);
    return filename;
}

static int
find_doc_idx(std::string filename)
{
    filename = get_cif_filename(filename);

    if (filename.find("stdin") != std::string::npos) {
        if (docs.size() == 1)
            return 0;
        else
            freesasa_fail(
                "In %s(), CIF input is from stdin but there are more than 1 gemmi documents to choose from. "
                "Unable to select correct doc. exiting...",
                __func__);
    }

    std::transform(filename.begin(), filename.end(), filename.begin(), tolower);
    for (int i = 0; i != docs.size(); ++i) {
        // only compare the first 4 characters of the doc.source value
        // this allows different file suffix and extensions
        if (filename.find(docs[i].source.substr(0,4)) != std::string::npos) {
            return i;
        }
    }
    return FREESASA_WARN;
}

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

static void
write_cif_block(std::ostream &out, gemmi::cif::Table &table,
                std::vector<std::string> &sasa_vals,
                std::vector<std::string> &sasa_radii)
{
    auto &loop = *table.get_loop();

    unsigned long orig_tag_size = loop.tags.size();
    unsigned long new_tag_size = orig_tag_size + 2;

    // Creates a new table full of empty strings with the correct number of dimensions
    // Outside vector size is the # of columns, inside vector size is the # of rows.
    std::vector<std::vector<std::string>> newCols(new_tag_size, {loop.length(), {"Empty"}});

    // Copies data from original columns to their respecitve column in the new table filled with empty strings.
    // Leaving only the new appended columns as empty strings
    for (unsigned int i = 0; i != orig_tag_size; ++i) {
        auto iCol = table.bloc.find_loop(loop.tags[i]);
        std::copy(iCol.begin(), iCol.end(), newCols[i].begin());
    }

    newCols[new_tag_size - 2] = std::move(sasa_vals);
    newCols[new_tag_size - 1] = std::move(sasa_radii);

    std::vector<std::string> new_tags{
        "_atom_site.FreeSASA_value",
        "_atom_site.FreeSASA_radius"};
    for (auto tag : new_tags)
        loop.tags.push_back(tag);

    loop.set_all_values(newCols);

    gemmi::cif::write_cif_block_to_stream(out, table.bloc);
}

static int
write_result(std::ostream &out, freesasa_node *root)
{
    freesasa_node *result{freesasa_node_children(root)};

    int prev_doc_idx = -1, doc_idx = 0;
    bool write = false;
    std::vector<std::string> sasa_vals, sasa_radii;

    while (result) {
        doc_idx = find_doc_idx(freesasa_node_name(result));
        if (doc_idx == FREESASA_WARN) {
            freesasa_warn("In %s(), unable to find gemmi doc for result node: %s. Skipping...",
                          __func__, freesasa_node_name(result));
            result = freesasa_node_next(result);
            continue;
        }

        auto &block = docs[doc_idx].sole_block();

        append_freesasa_params_to_block(block, result);

        if (append_freesasa_result_summary_to_block(block, result) == FREESASA_FAIL) {
            return freesasa_fail("Unable to build CIF output");
        }

        auto table = block.find("_atom_site.", atom_site_columns);
        if (prev_doc_idx != doc_idx) {
            sasa_vals = std::vector<std::string>{table.length(), "?"};
            sasa_radii = std::vector<std::string>{table.length(), "?"};
        }

        populate_freesasa_result_vectors(table, result, sasa_vals, sasa_radii);

        prev_doc_idx = doc_idx;
        result = freesasa_node_next(result);

        if (!result) {
            // There is no next result so write out current file.
            write = true;
        } else if (find_doc_idx(freesasa_node_name(result)) != prev_doc_idx) {
            // Next result node is from a new file so write out current file.
            write = true;
        } else {
            // Next result node is from the same doc so do not write current file yet.
            write = false;
        }

        if (write) write_cif_block(out, table, sasa_vals, sasa_radii);
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
