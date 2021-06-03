import os

import pandas
import numpy

from lxml import etree
import statsmodels.formula.api as smf
from patsy import ModelDesc


def read_xml(fn_in):
    tree = etree.parse(fn_in)
    return tree


def validation_results(fn_in):
    data = pandas.read_csv(fn_in)
    data = data.loc[~numpy.isnan(data['Error'])]
    return data


def write_xml(tree, fn_out):
    if os.path.exists(fn_out):
        raise ValueError("{0} already exists. Not overwriting!".format(fn_out))
    with open(fn_out, "wb") as fid:
        fid.write(etree.tostring(tree, encoding="UTF-8", xml_declaration=True))


def update_connectivity_recipe(tree, set_key, sm_result, min_val, max_val):
    set_fac = sm_result.params["data"]
    root = tree.getroot()
    for c in root.getchildren():
        if c.get(set_key) is not None:
            v_out = set_fac * float(c.get(set_key))
            v_out = numpy.maximum(min_val,
                                  numpy.minimum(max_val, v_out))
            c.set(set_key, str(v_out))

    return tree


def linear_fit_model(data):
    sm_model = ModelDesc.from_formula("ref ~ data - 1")  # TODO: model with offset?
    sm_result = smf.ols(sm_model, {"ref": data["Mean (ref.)"], "data": data["Mean (data)"]}).fit()
    print(sm_result.params)
    return sm_result


def main():
    import sys
    import getopt
    opts, args = getopt.getopt(sys.argv[1:], "i:o:b:s:")
    opts = dict(opts)

    if "-i" not in opts:
        print("""Need to specify an input s2f-recipe file:
        {0} [...] -i builderConnectivityRecipeAllPathways.xml""".format(__file__))
        sys.exit(2)
    tree = read_xml(opts["-i"])

    if "-b" not in opts and "-s" not in opts:
        print("""Neither results on bouton density nor synapses per connection provided.
        No adjustments will be performed!""")

    if "-o" not in opts:
        import os
        fn, ext = os.path.splitext(opts["-i"])
        opts["-o"] = fn + "_adjusted" + ext

    if "-b" in opts:
        data = validation_results(opts["-b"])
        recipe_set_key = "bouton_reduction_factor"
        min_val = 0.0
        max_val = 1.0
        print("Adjusting {0} based on results in {1}".format(recipe_set_key, opts["-b"]))
        mdl_fit = linear_fit_model(data)
        tree = update_connectivity_recipe(tree, recipe_set_key, mdl_fit, min_val, max_val)

    if "-s" in opts:
        data = validation_results(opts["-s"])
        recipe_set_key = "mean_syns_connection"
        min_val = 1.0
        max_val = 1E20
        print("Adjusting {0} based on results in {1}".format(recipe_set_key, opts["-s"]))
        mdl_fit = linear_fit_model(data)
        tree = update_connectivity_recipe(tree, recipe_set_key, mdl_fit, min_val, max_val)

    write_xml(tree, opts["-o"])


if __name__ == "__main__":
    main()
