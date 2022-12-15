# for posterity, my failed biotite attempts

# Generate example alignment
# (the same as in the bacterial luciferase example)
query =   entrez.SimpleQuery("luxA", "Gene Name") \
        & entrez.SimpleQuery("srcdb_swiss-prot", "Properties")
uids = entrez.search(query, db_name="protein")
fasta_file = fasta.FastaFile.read(entrez.fetch_single_file(
    uids, None, db_name="protein", ret_type="fasta"
))
sequences = [seq.ProteinSequence(seq_str) for seq_str in fasta_file.values()]
matrix = align.SubstitutionMatrix.std_protein_matrix()
alignment, order, _, _ = align.align_multiple(sequences, matrix)
# Order alignment according to the guide tree
alignment = alignment[:, order]
alignment = alignment[220:300]

alphabet = seq.ProteinSequence.alphabet

# define own version of plotting function here > this works but I can't plot the colours sep

def plot_alignment_occupancy(axes, alignment, symbols_per_line=50,
                              show_numbers=False, number_size=None,
                              number_functions=None,
                              labels=None, label_size=None,
                              show_line_position=False,
                              spacing=1,
                              color_scheme=None, color_symbols=False,
                              symbol_size=None, symbol_param=None,
                              symbol_spacing=None):
    alphabet = alignment.sequences[0].get_alphabet()
    symbol_plotter = biotite.sequence.graphics.LetterTypePlotter(
        axes, alphabet, font_size=symbol_size, font_param=symbol_param,
        color_symbols=color_symbols, color_scheme=None )
        #color_scheme
    biotite.sequence.graphics.plot_alignment(
        axes=axes, alignment=alignment, symbol_plotter=symbol_plotter,
        symbols_per_line=symbols_per_line,
        show_numbers=show_numbers, number_size=number_size,
        number_functions=number_functions,
        labels=labels, label_size=label_size,
        show_line_position=show_line_position,
        spacing=spacing) 

fig = plt.figure(figsize=(8.0, 2))
gridspec = GridSpec(2, count)
ax = fig.add_subplot()
ax.set_ylabel(name)
print(type(alignment))
print(alignment)
alignment_part = alignment[:40]
plot_alignment_occupancy( ax, alignment_part, symbols_per_line=len(alignment_part) )
fig.tight_layout()
fig.subplots_adjust(wspace=0)
plt.show()
