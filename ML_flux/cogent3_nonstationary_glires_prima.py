from cogent3 import load_tree, load_aligned_seqs, get_model, make_tree, available_models
from cogent3.evolve.substitution_model import EmpiricalProteinMatrix
from cogent3.parse.paml_matrix import PamlMatrixParser

aln = load_aligned_seqs("glires_primatomorpha_alignment.fa", moltype="protein")
matrix_file = open("Q.mammal")
empirical_matrix , empirical_frequencies = PamlMatrixParser(matrix_file)
sm = EmpiricalProteinMatrix(empirical_matrix, empirical_frequencies)
tree = load_tree("glires_primatomorpha.tre")
tree = tree.get_sub_tree(aln.names)
primate_edges = (tree.get_edge_names("CHLOR_SAB","GALEO_VAR", outgroup_name="OCHOT_PRI", clade=True, stem=True))
rodent_edges = (tree.get_edge_names("MICRO_OCH","OCHOT_PRI", outgroup_name="GALEO_VAR", clade=True, stem=True))

opt_settings = dict(max_restarts=5, tolerance=1e-9, show_progress=False)
lf = sm.make_likelihood_function(tree)
lf.set_alignment(aln)
lf.set_param_rule("mprobs",edges=primate_edges, clade=True, stem=True)
lf.set_param_rule("mprobs",edges=rodent_edges, clade=True, stem=True)
lf.set_param_rule("mprobs",edge="root")
print(tree.ascii_art())
lf.optimise(**opt_settings)
print(lf)