#!/usr/bin/env python3
"""Build an SBML Level 3 Version 2 document from a JSON spec produced by
dMod's `exportSbml()`. Mirrors the import helper `sbmlImport.py` but in
the opposite direction. Uses `libsbml.parseL3Formula` to lift the text-form
kinetic laws into MathML ASTs.
"""

import sys
import json
import libsbml


def _check(rc, ctx):
    if rc != libsbml.LIBSBML_OPERATION_SUCCESS:
        raise RuntimeError("libsbml failure at %s: %d" % (ctx, rc))


def build_sbml(spec):
    ns = libsbml.SBMLNamespaces(3, 2)
    document = libsbml.SBMLDocument(ns)
    model = document.createModel()
    _check(model.setId(spec.get("modelId", "dMod_export")), "setId")

    for c in spec.get("compartments", []):
        comp = model.createCompartment()
        _check(comp.setId(c["id"]), "compartment.setId")
        comp.setSize(c.get("size", 1.0))
        comp.setSpatialDimensions(c.get("spatialDimensions", 3))
        comp.setConstant(True)
        # A symbolic volume is the compartment's size, not a factor hidden in
        # the kinetic laws: SBML divides by the size, so it has to carry it.
        formula = c.get("sizeAssignment")
        if formula is not None:
            ia = model.createInitialAssignment()
            _check(ia.setSymbol(c["id"]), "initialAssignment.setSymbol")
            ast = libsbml.parseL3Formula(formula)
            if ast is None:
                raise ValueError(
                    "Could not parse compartment size for %r: %s -- got %r"
                    % (c["id"], libsbml.getLastParseL3Error(), formula))
            ia.setMath(ast)

    for s in spec.get("species", []):
        sp = model.createSpecies()
        _check(sp.setId(s["id"]), "species.setId")
        _check(sp.setCompartment(s["compartment"]), "species.setCompartment")
        # Symbolic initials become InitialAssignments (created below); numeric
        # initials use initialConcentration. SBML lets both coexist -- the
        # InitialAssignment wins at sim time -- but we keep them mutually
        # exclusive for cleanliness on roundtrip.
        amount = bool(s.get("hasOnlySubstanceUnits", False))
        if "initialAssignment" not in s:
            if amount:
                sp.setInitialAmount(s.get("initialAmount", 0.0))
            else:
                sp.setInitialConcentration(s.get("initialConcentration", 0.0))
        sp.setHasOnlySubstanceUnits(amount)
        sp.setBoundaryCondition(False)
        sp.setConstant(False)

    for s in spec.get("species", []):
        formula = s.get("initialAssignment")
        if formula is None:
            continue
        ia = model.createInitialAssignment()
        _check(ia.setSymbol(s["id"]), "initialAssignment.setSymbol")
        ast = libsbml.parseL3Formula(formula)
        if ast is None:
            raise ValueError(
                "Could not parse initialAssignment for %r: %s -- got %r"
                % (s["id"], libsbml.getLastParseL3Error(), formula))
        ia.setMath(ast)

    for p in spec.get("parameters", []):
        par = model.createParameter()
        _check(par.setId(p["id"]), "parameter.setId")
        par.setValue(p.get("value", 0.0))
        par.setConstant(True)

    for r in spec.get("reactions", []):
        rxn = model.createReaction()
        _check(rxn.setId(r["id"]), "reaction.setId")
        rxn.setReversible(False)
        for reactant in r.get("reactants", []):
            sref = rxn.createReactant()
            sref.setSpecies(reactant["species"])
            sref.setStoichiometry(float(reactant["stoich"]))
            sref.setConstant(True)
        for product in r.get("products", []):
            sref = rxn.createProduct()
            sref.setSpecies(product["species"])
            sref.setStoichiometry(float(product["stoich"]))
            sref.setConstant(True)
        kl = rxn.createKineticLaw()
        ast = libsbml.parseL3Formula(r["kineticLaw"])
        if ast is None:
            raise ValueError(
                "Could not parse kinetic law %r: %s"
                % (r["kineticLaw"], libsbml.getLastParseL3Error())
            )
        kl.setMath(ast)

    return document


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: %s SPEC.json" % __file__)
        sys.exit(1)

    with open(sys.argv[1]) as fh:
        spec = json.load(fh)

    doc = build_sbml(spec)
    libsbml.writeSBMLToFile(doc, spec["outfile"])
